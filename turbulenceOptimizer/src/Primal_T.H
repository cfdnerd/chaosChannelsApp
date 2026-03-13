auto clampTemperature =
[&]() -> label
{
  if (!clipTemperatureField)
  {
    return 0;
  }
  label clippedCount(0);
  scalarField& tInternal = T.primitiveFieldRef();
  forAll(tInternal, cellI)
  {
    const scalar originalValue = tInternal[cellI];
    const scalar boundedValue =
      Foam::min(Foam::max(originalValue, temperatureLowerBound), temperatureUpperBound);
    if (boundedValue != originalValue)
    {
      tInternal[cellI] = boundedValue;
      clippedCount++;
    }
  }
  T.correctBoundaryConditions();
  return clippedCount;
};

auto evaluateThermalField =
[&](scalar& tMinGlobal, scalar& tMaxGlobal, label& invalidCountGlobal)
{
  scalar tMinLocal(GREAT);
  scalar tMaxLocal(-GREAT);
  label invalidCountLocal(0);
  const scalarField& tInternal = T.primitiveField();
  forAll(tInternal, cellI)
  {
    const scalar value = tInternal[cellI];
    if (!std::isfinite(static_cast<double>(value)))
    {
      invalidCountLocal++;
      continue;
    }
    tMinLocal = Foam::min(tMinLocal, value);
    tMaxLocal = Foam::max(tMaxLocal, value);
  }
  reduce(tMinLocal, minOp<scalar>());
  reduce(tMaxLocal, maxOp<scalar>());
  reduce(invalidCountLocal, sumOp<label>());
  tMinGlobal = tMinLocal;
  tMaxGlobal = tMaxLocal;
  invalidCountGlobal = invalidCountLocal;
};

T = T2;
DTeff = DT + turbulence->nut()/thermalTurbulentPrandtl;
DTeff.correctBoundaryConditions();
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(DTeff, T)
);
if(DDHS)
{
 TEqn -= thermalSource*(1-xh);
}
else
{
 TEqn -= thermalSource;
}
const auto TLinearSolvePerformancePrimary = TEqn.solve();
recordLinearSolve
(
  TLinearSolveStats,
  TLinearSolvePerformancePrimary.initialResidual(),
  TLinearSolvePerformancePrimary.finalResidual(),
  TLinearSolvePerformancePrimary.nIterations()
);
thermalFallbackLastPrimaryFinalResidual = TLinearSolvePerformancePrimary.finalResidual();
thermalFallbackLastPrimaryIterations = TLinearSolvePerformancePrimary.nIterations();

if (thermalUnderRelaxation < 1.0)
{
  T = thermalUnderRelaxation*T + (1.0 - thermalUnderRelaxation)*T2;
  T.correctBoundaryConditions();
}
label tClippedCount = clampTemperature();

scalar tMinValue(0.0);
scalar tMaxValue(0.0);
label tInvalidCount(0);
evaluateThermalField(tMinValue, tMaxValue, tInvalidCount);
const bool thermalFieldDivergedPrimary =
  (tInvalidCount > 0)
  ||
  (!std::isfinite(static_cast<double>(tMinValue)))
  ||
  (!std::isfinite(static_cast<double>(tMaxValue)))
  ||
  (Foam::mag(tMinValue) > thermalFieldAbsMax)
  ||
  (Foam::mag(tMaxValue) > thermalFieldAbsMax);

const bool thermalFallbackTrigger =
  enableThermalFallback
  &&
  (
    TLinearSolvePerformancePrimary.finalResidual() > thermalFallbackResidualThreshold
    ||
    TLinearSolvePerformancePrimary.nIterations() >= thermalFallbackIterationThreshold
    ||
    thermalFieldDivergedPrimary
  );

bool thermalSolveHealthy = !thermalFieldDivergedPrimary;
if (thermalFallbackTrigger)
{
  thermalFallbackTriggerCount++;
  const dictionary& solversDict = mesh.solutionDict().subDict("solvers");
  if (solversDict.found("TFallback"))
  {
    thermalSolveHealthy = false;
    for (label fallbackAttempt = 0; fallbackAttempt < thermalFallbackMaxAttempts; ++fallbackAttempt)
    {
      T = thermalFallbackBlend*T2 + (1.0 - thermalFallbackBlend)*T;
      T.correctBoundaryConditions();

      volScalarField DTFallback(DTeff);
      if (thermalFallbackDiffusivityFactor > 1.0)
      {
        DTFallback *= thermalFallbackDiffusivityFactor;
      }
      fvScalarMatrix TEqnFallback
      (
          fvm::div(phi,T)
         -fvm::laplacian(DTFallback, T)
      );
      if(DDHS)
      {
        TEqnFallback -= thermalSource*(1-xh);
      }
      else
      {
        TEqnFallback -= thermalSource;
      }
      const auto TLinearSolvePerformanceFallback = TEqnFallback.solve(solversDict.subDict("TFallback"));
      recordLinearSolve
      (
        TFallbackLinearSolveStats,
        TLinearSolvePerformanceFallback.initialResidual(),
        TLinearSolvePerformanceFallback.finalResidual(),
        TLinearSolvePerformanceFallback.nIterations()
      );
      thermalFallbackLastFallbackFinalResidual = TLinearSolvePerformanceFallback.finalResidual();
      thermalFallbackLastFallbackIterations = TLinearSolvePerformanceFallback.nIterations();

      if (thermalUnderRelaxation < 1.0)
      {
        T = thermalUnderRelaxation*T + (1.0 - thermalUnderRelaxation)*T2;
        T.correctBoundaryConditions();
      }
      tClippedCount += clampTemperature();
      evaluateThermalField(tMinValue, tMaxValue, tInvalidCount);
      const bool thermalFieldDivergedFallback =
        (tInvalidCount > 0)
        ||
        (!std::isfinite(static_cast<double>(tMinValue)))
        ||
        (!std::isfinite(static_cast<double>(tMaxValue)))
        ||
        (Foam::mag(tMinValue) > thermalFieldAbsMax)
        ||
        (Foam::mag(tMaxValue) > thermalFieldAbsMax);
      if
      (
        !thermalFieldDivergedFallback
        &&
        TLinearSolvePerformanceFallback.finalResidual() <= thermalFallbackResidualThreshold
      )
      {
        thermalFallbackSuccessCount++;
        thermalSolveHealthy = true;
        break;
      }
    }
    if (!thermalSolveHealthy)
    {
      thermalFallbackFailureCount++;
    }
  }
  else
  {
    thermalFallbackMissingSolverCount++;
    thermalSolveHealthy = false;
    if (Pstream::master())
    {
      WarningInFunction
        << "Thermal fallback triggered but solver entry 'TFallback' is missing in fvSolution."
        << endl;
    }
  }
}

thermalLastMin = tMinValue;
thermalLastMax = tMaxValue;

if (enableThermalDivergenceGuard && !thermalSolveHealthy)
{
  stopOptimization = true;
  optimizerStopReason = "thermal_divergence";
  T = T2;
  T.correctBoundaryConditions();
  if (Pstream::master())
  {
    WarningInFunction
      << "Thermal solve failed stability guard at optimizer iteration " << opt
      << ". Clipped cells=" << tClippedCount
      << ", min/max(T)=" << tMinValue << "/" << tMaxValue
      << ". Optimization marked for stop."
      << endl;
  }
}
else
{
  T2 = T;
}

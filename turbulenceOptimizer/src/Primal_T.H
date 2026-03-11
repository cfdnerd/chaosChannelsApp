T=T2;
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(DT, T)
);
if(DDHS)
{
 TEqn-=thermalSource*(1-xh);
}
else
{
 TEqn-=thermalSource;
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

const bool thermalFallbackTrigger =
  enableThermalFallback
  &&
  (
    TLinearSolvePerformancePrimary.finalResidual() > thermalFallbackResidualThreshold
    ||
    TLinearSolvePerformancePrimary.nIterations() >= thermalFallbackIterationThreshold
  );

if (thermalFallbackTrigger)
{
  thermalFallbackTriggerCount++;
  const dictionary& solversDict = mesh.solutionDict().subDict("solvers");
  if (solversDict.found("TFallback"))
  {
    // Damp the thermal state with the previous iterate before fallback solve.
    T = thermalFallbackBlend*T2 + (1.0 - thermalFallbackBlend)*T;
    T.correctBoundaryConditions();

    fvScalarMatrix TEqnFallback
    (
        fvm::div(phi,T)
       -fvm::laplacian(DT, T)
    );
    if(DDHS)
    {
      TEqnFallback-=thermalSource*(1-xh);
    }
    else
    {
      TEqnFallback-=thermalSource;
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
    if (TLinearSolvePerformanceFallback.finalResidual() <= TLinearSolvePerformancePrimary.finalResidual())
    {
      thermalFallbackSuccessCount++;
    }
  }
  else
  {
    thermalFallbackMissingSolverCount++;
    if (Pstream::master())
    {
      WarningInFunction
        << "Thermal fallback triggered but solver entry 'TFallback' is missing in fvSolution."
        << endl;
    }
  }
}

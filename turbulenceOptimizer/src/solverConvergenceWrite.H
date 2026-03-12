const scalar expectedPrimalCalls =
  static_cast<scalar>(primalFlowIterations);
const scalar expectedAdjointCalls =
  static_cast<scalar>(adjointFlowIterations);
const scalar paCorrectionsPerUa =
  UaLinearSolveStats.callCount > 0
    ? static_cast<scalar>(paLinearSolveStats.callCount)
     /static_cast<scalar>(UaLinearSolveStats.callCount)
    : 0.0;
const scalar pbCorrectionsPerUb =
  UbLinearSolveStats.callCount > 0
    ? static_cast<scalar>(pbLinearSolveStats.callCount)
     /static_cast<scalar>(UbLinearSolveStats.callCount)
    : 0.0;

scalar UaMeanMagConv(0.0), UaMaxMagConv(0.0), UbMeanMagConv(0.0), UbMaxMagConv(0.0);
scalar paMeanConv(0.0), paMinConv(0.0), paMaxConv(0.0), pbMeanConv(0.0), pbMinConv(0.0), pbMaxConv(0.0);
scalar TbMeanConv(0.0), TbMinConv(0.0), TbMaxConv(0.0);
{
  auto vectorMagnitudeStats =
  [&](const volVectorField& field, scalar& meanMagnitude, scalar& maxMagnitude)
  {
    const vectorField& internal = field.primitiveField();
    scalar sumLocal(0.0);
    scalar maxLocal(0.0);
    scalar countGlobal = internal.size();
    forAll(internal, cellI)
    {
      const scalar magValue = Foam::mag(internal[cellI]);
      sumLocal += magValue;
      maxLocal = Foam::max(maxLocal, magValue);
    }
    reduce(sumLocal, sumOp<scalar>());
    reduce(maxLocal, maxOp<scalar>());
    reduce(countGlobal, sumOp<scalar>());
    meanMagnitude = sumLocal/Foam::max(countGlobal, scalar(SMALL));
    maxMagnitude = maxLocal;
  };

  auto scalarFieldStats =
  [&](const volScalarField& field, scalar& meanValue, scalar& minValue, scalar& maxValue)
  {
    const scalarField& internal = field.primitiveField();
    scalar sumLocal(0.0);
    scalar minLocal(GREAT);
    scalar maxLocal(-GREAT);
    scalar countGlobal = internal.size();
    forAll(internal, cellI)
    {
      const scalar value = internal[cellI];
      sumLocal += value;
      minLocal = Foam::min(minLocal, value);
      maxLocal = Foam::max(maxLocal, value);
    }
    reduce(sumLocal, sumOp<scalar>());
    reduce(minLocal, minOp<scalar>());
    reduce(maxLocal, maxOp<scalar>());
    reduce(countGlobal, sumOp<scalar>());
    meanValue = sumLocal/Foam::max(countGlobal, scalar(SMALL));
    minValue = minLocal;
    maxValue = maxLocal;
  };

  vectorMagnitudeStats(Ua, UaMeanMagConv, UaMaxMagConv);
  vectorMagnitudeStats(Ub, UbMeanMagConv, UbMaxMagConv);
  scalarFieldStats(pa, paMeanConv, paMinConv, paMaxConv);
  scalarFieldStats(pb, pbMeanConv, pbMinConv, pbMaxConv);
  scalarFieldStats(Tb, TbMeanConv, TbMinConv, TbMaxConv);
}

const bool hasThermalFallbackSolver =
  mesh.solutionDict().subDict("solvers").found("TFallback");
label thermalFallbackTriggerCountConv = thermalFallbackTriggerCount;
label thermalFallbackSuccessCountConv = thermalFallbackSuccessCount;
label thermalFallbackMissingSolverCountConv = thermalFallbackMissingSolverCount;
label thermalFallbackFailureCountConv = thermalFallbackFailureCount;
label thermalFallbackLastPrimaryIterationsConv = thermalFallbackLastPrimaryIterations;
label thermalFallbackLastFallbackIterationsConv = thermalFallbackLastFallbackIterations;
scalar thermalFallbackLastPrimaryFinalResidualConv = thermalFallbackLastPrimaryFinalResidual;
scalar thermalFallbackLastFallbackFinalResidualConv = thermalFallbackLastFallbackFinalResidual;
scalar thermalLastMinConv = thermalLastMin;
scalar thermalLastMaxConv = thermalLastMax;
label turbulenceKClipCountConv = turbulenceKClipCount;
label turbulenceEpsilonClipCountConv = turbulenceEpsilonClipCount;
label turbulenceNutClipCountConv = turbulenceNutClipCount;
reduce(thermalFallbackTriggerCountConv, maxOp<label>());
reduce(thermalFallbackSuccessCountConv, maxOp<label>());
reduce(thermalFallbackMissingSolverCountConv, maxOp<label>());
reduce(thermalFallbackFailureCountConv, maxOp<label>());
reduce(thermalFallbackLastPrimaryIterationsConv, maxOp<label>());
reduce(thermalFallbackLastFallbackIterationsConv, maxOp<label>());
reduce(thermalFallbackLastPrimaryFinalResidualConv, maxOp<scalar>());
reduce(thermalFallbackLastFallbackFinalResidualConv, maxOp<scalar>());
reduce(thermalLastMinConv, minOp<scalar>());
reduce(thermalLastMaxConv, maxOp<scalar>());
reduce(turbulenceKClipCountConv, sumOp<label>());
reduce(turbulenceEpsilonClipCountConv, sumOp<label>());
reduce(turbulenceNutClipCountConv, sumOp<label>());

if (Pstream::master())
{
  std::ofstream solverConvergenceLogFile(solverConvergenceFile, std::ios::app);

  auto writeTableHeader =
  [&]()
  {
    solverConvergenceLogFile
      << "  "
      << std::left
      << std::setw(22) << "Solve"
      << std::setw(6) << "Stat"
      << std::right
      << std::setw(7) << "Calls"
      << std::setw(14) << "AvgInit"
      << std::setw(14) << "MaxInit"
      << std::setw(14) << "AvgFinal"
      << std::setw(14) << "MaxFinal"
      << std::setw(10) << "AvgIter"
      << std::setw(10) << "MaxIter"
      << std::setw(14) << "LastFinal"
      << std::setw(10) << "LastItr"
      << "\n"
      << "  "
      << std::string(22 + 6 + 7 + 14 + 14 + 14 + 14 + 10 + 10 + 14 + 10, '-')
      << "\n";
  };

  auto writeSolveStats =
  [&](const char* solveName, const LinearSolveStats& stats)
  {
    const scalar avgInitialResidual =
      stats.callCount > 0 ? stats.initialResidualSum/static_cast<scalar>(stats.callCount) : 0;
    const scalar avgFinalResidual =
      stats.callCount > 0 ? stats.finalResidualSum/static_cast<scalar>(stats.callCount) : 0;
    const scalar avgIterations =
      stats.callCount > 0 ? static_cast<scalar>(stats.iterationSum)/static_cast<scalar>(stats.callCount) : 0;
    const char* status = classifyLinearSolveStatus(stats);

    solverConvergenceLogFile
      << "  "
      << std::left << std::setw(22) << solveName
      << std::setw(6) << status
      << std::right
      << std::setw(7) << stats.callCount
      << std::scientific << std::setprecision(3)
      << std::setw(14) << avgInitialResidual
      << std::setw(14) << stats.initialResidualMax
      << std::setw(14) << avgFinalResidual
      << std::setw(14) << stats.finalResidualMax
      << std::fixed << std::setprecision(2)
      << std::setw(10) << avgIterations
      << std::setw(10) << stats.iterationMax
      << std::scientific << std::setprecision(3)
      << std::setw(14) << stats.lastFinalResidual
      << std::fixed << std::setprecision(0)
      << std::setw(10) << stats.lastIterations
      << "\n";
  };

  auto writeCorrectionTableHeader =
  [&]()
  {
    solverConvergenceLogFile
      << "  "
      << std::left
      << std::setw(22) << "Field"
      << std::setw(6) << "Stat"
      << std::right
      << std::setw(7) << "Calls"
      << std::setw(16) << "AvgMean|d|"
      << std::setw(16) << "MaxMean|d|"
      << std::setw(16) << "Max|d|"
      << std::setw(16) << "AvgRelL2"
      << std::setw(16) << "MaxRelL2"
      << std::setw(16) << "LastRelL2"
      << "\n"
      << "  "
      << std::string(22 + 6 + 7 + 16 + 16 + 16 + 16 + 16 + 16, '-')
      << "\n";
  };

  auto writeCorrectionStats =
  [&](const char* fieldName, const FieldCorrectionStats& stats)
  {
    const scalar avgMeanAbsDelta =
      stats.callCount > 0 ? stats.meanAbsDeltaSum/static_cast<scalar>(stats.callCount) : 0.0;
    const scalar avgRelL2Delta =
      stats.callCount > 0 ? stats.relL2DeltaSum/static_cast<scalar>(stats.callCount) : 0.0;
    const char* status = stats.callCount > 0 ? "OK" : "SKIP";

    solverConvergenceLogFile
      << "  "
      << std::left << std::setw(22) << fieldName
      << std::setw(6) << status
      << std::right
      << std::setw(7) << stats.callCount
      << std::scientific << std::setprecision(3)
      << std::setw(16) << avgMeanAbsDelta
      << std::setw(16) << stats.meanAbsDeltaMax
      << std::setw(16) << stats.maxAbsDeltaMax
      << std::setw(16) << avgRelL2Delta
      << std::setw(16) << stats.relL2DeltaMax
      << std::setw(16) << stats.lastRelL2Delta
      << "\n";
  };

  solverConvergenceLogFile
    << "\n"
    << "=========================================================================================================\n"
    << " Solver Convergence Report\n"
    << "---------------------------------------------------------------------------------------------------------\n"
    << "  Optimization Iteration : " << opt << "\n"
    << "  Status Legend         : OK=healthy  OBS=watch  WARN=likely issue  SKIP=not called\n"
    << "=========================================================================================================\n"
    << "[Primary Equations]\n";
  writeTableHeader();
  writeSolveStats("U momentum", ULinearSolveStats);
  writeSolveStats("p pressure", pLinearSolveStats);
  writeSolveStats("T thermal", TLinearSolveStats);
  writeSolveStats("T thermal fallback", TFallbackLinearSolveStats);
  writeSolveStats("Tb adjoint thermal", TbLinearSolveStats);
  writeSolveStats("Ub adjoint momentum", UbLinearSolveStats);
  writeSolveStats("pb adjoint pressure", pbLinearSolveStats);
  writeSolveStats("Ua adjoint momentum", UaLinearSolveStats);
  writeSolveStats("pa adjoint pressure", paLinearSolveStats);

  solverConvergenceLogFile
    << "\n"
    << "[Filter And Projection Equations]\n";
  writeTableHeader();
  writeSolveStats("xp filter", xpFilterLinearSolveStats);
  writeSolveStats("fsens filter", fsensFilterLinearSolveStats);
  writeSolveStats("gsensPower filter", gsensPowerFilterLinearSolveStats);
  writeSolveStats("gsensVol filter", gsensVolFilterLinearSolveStats);

  solverConvergenceLogFile
    << "\n"
    << "[Turbulence-Primal Corrections]\n";
  writeCorrectionTableHeader();
  writeCorrectionStats("k correction", kTurbulenceCorrectionStats);
  writeCorrectionStats("epsilon correction", epsilonTurbulenceCorrectionStats);
  writeCorrectionStats("nut correction", nutTurbulenceCorrectionStats);

  solverConvergenceLogFile
    << "\n"
    << "[Turbulence And Adjoint Diagnostics]\n"
    << "  Turbulence model          : " << turbulence->type() << "\n"
    << "  turbulenceProMode         : " << (turbulenceProMode ? "true" : "false") << "\n"
    << "  Turbulence corrector intv : " << turbulenceCorrectorInterval << "\n"
    << "  Primal flow solves        : " << ULinearSolveStats.callCount
    << " / expected ~ " << expectedPrimalCalls << "\n"
    << "  Adjoint momentum solves   : Ua=" << UaLinearSolveStats.callCount
    << ", Ub=" << UbLinearSolveStats.callCount
    << " / expected each ~ " << expectedAdjointCalls << "\n"
    << "  Thermal fallback enabled  : " << (enableThermalFallback ? "true" : "false")
    << ", solver configured: " << (hasThermalFallbackSolver ? "true" : "false") << "\n"
    << "  Thermal fallback triggers : " << thermalFallbackTriggerCountConv
    << ", successful improves: " << thermalFallbackSuccessCountConv
    << ", missing solver: " << thermalFallbackMissingSolverCountConv
    << ", failed: " << thermalFallbackFailureCountConv << "\n"
    << "  T primary/fallback        : finalRes="
    << thermalFallbackLastPrimaryFinalResidualConv << " / "
    << thermalFallbackLastFallbackFinalResidualConv
    << ", iters="
    << thermalFallbackLastPrimaryIterationsConv << " / "
    << thermalFallbackLastFallbackIterationsConv << "\n"
    << "  T guarded min/max         : " << thermalLastMinConv << " / " << thermalLastMaxConv << "\n"
    << "  Turbulence clips (k/e/n)  : "
    << turbulenceKClipCountConv << " / "
    << turbulenceEpsilonClipCountConv << " / "
    << turbulenceNutClipCountConv << "\n"
    << "  Adjoint pressure loops    : pa/Ua=" << paCorrectionsPerUa
    << ", pb/Ub=" << pbCorrectionsPerUb << "\n"
    << "  Frozen turbulence adjoint : " << (frozenTurbulenceAdjoint ? "true" : "false") << "\n"
    << "  |Ua| mean/max             : " << UaMeanMagConv << " / " << UaMaxMagConv << "\n"
    << "  |Ub| mean/max             : " << UbMeanMagConv << " / " << UbMaxMagConv << "\n"
    << "  pa mean/min/max           : " << paMeanConv << " / " << paMinConv << " / " << paMaxConv << "\n"
    << "  pb mean/min/max           : " << pbMeanConv << " / " << pbMinConv << " / " << pbMaxConv << "\n"
    << "  Tb mean/min/max           : " << TbMeanConv << " / " << TbMinConv << " / " << TbMaxConv << "\n";
}

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
}

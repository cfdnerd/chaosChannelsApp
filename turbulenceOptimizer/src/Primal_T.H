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

if (thermalUnderRelaxation < 1.0)
{
  T = thermalUnderRelaxation*T + (1.0 - thermalUnderRelaxation)*T2;
  T.correctBoundaryConditions();
}
T2 = T;

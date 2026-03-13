T=T2;
volScalarField DTeff("DTeff", DT + (xh*(1+qu)/(qu+xh))*turbulence->nut()/Prt);
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(DTeff, T)
);
if(DDHS)
{
 TEqn-=thermalSource*(1-xh);
}
else
{
 TEqn-=thermalSource;
}
const auto TLinearSolvePerformance = TEqn.solve();
recordLinearSolve
(
  TLinearSolveStats,
  TLinearSolvePerformance.initialResidual(),
  TLinearSolvePerformance.finalResidual(),
  TLinearSolvePerformance.nIterations()
);

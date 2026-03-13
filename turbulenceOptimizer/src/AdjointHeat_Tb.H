Tb = Tb2;
sumT=0;
maxT=gMax(T);
if (objFunction==2 && maxT > SMALL)
{
  for(i=0;i<n;i++)
  {
     sumT=sumT+std::pow(T[i]/maxT,Pnorm);
  }
  reduce(sumT, sumOp<scalar>());
  sumT=std::pow(sumT,1.0/Pnorm-1);
}
DTeff.correctBoundaryConditions();

fvScalarMatrix TEqnb
(
   -fvm::div(-phi,Tb) 
   +fvm::laplacian(DTeff, Tb)
);
if (objFunction==1)
{
   TEqnb -= unitThermalSource;
}
if (objFunction==2)
{
   if (maxT > SMALL)
   {
      TEqnb -= b2*sumT*Foam::pow(T/maxT,Pnorm-1);
   }
}
if (enableAdjointThermalMatrixRelaxation)
{
  TEqnb.relax();
}
const auto TbLinearSolvePerformance = TEqnb.solve();
recordLinearSolve
(
  TbLinearSolveStats,
  TbLinearSolvePerformance.initialResidual(),
  TbLinearSolvePerformance.finalResidual(),
  TbLinearSolvePerformance.nIterations()
);
if (adjointThermalUnderRelaxation < 1.0)
{
  Tb = adjointThermalUnderRelaxation*Tb + (1.0 - adjointThermalUnderRelaxation)*Tb2;
  Tb.correctBoundaryConditions();
}
Tb2 = Tb;

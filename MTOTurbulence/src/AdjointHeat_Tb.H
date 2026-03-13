//Tb=T2;
sumT=0;
maxT=gMax(T);
scalar safeMaxT=Foam::max(maxT, 1e-10);
for(i=0;i<n;i++)
{
   sumT=sumT+std::pow(Foam::max(T[i]/safeMaxT, 0.0),Pnorm);
}      
reduce(sumT, sumOp<scalar>());
sumT=std::pow(Foam::max(sumT, 1e-15),1.0/Pnorm-1);

fvScalarMatrix TEqnb
(
   -fvm::div(-phi,Tb) 
   +fvm::laplacian(DT + (xh*(1+qu)/(qu+xh))*turbulence->nut()/Prt, Tb)
);
if (objFunction==1)
{
   TEqnb -= unitThermalSource;
}
if (objFunction==2)
{
   TEqnb -= b2*sumT*Foam::pow(Foam::max(T/safeMaxT, scalar(0.0)),Pnorm-1);
}
const auto TbLinearSolvePerformance = TEqnb.solve();
recordLinearSolve
(
  TbLinearSolveStats,
  TbLinearSolvePerformance.initialResidual(),
  TbLinearSolvePerformance.finalResidual(),
  TbLinearSolvePerformance.nIterations()
);

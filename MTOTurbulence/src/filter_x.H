
//***********************************//PDE filter
const auto xpFilterLinearSolvePerformance =
solve(fvm::laplacian(xp)-fvm::Sp(b,xp)+x*b);
recordLinearSolve
(
  xpFilterLinearSolveStats,
  xpFilterLinearSolvePerformance.initialResidual(),
  xpFilterLinearSolvePerformance.finalResidual(),
  xpFilterLinearSolvePerformance.nIterations()
);
//***********************************//Heaviside_filter
del=Foam::min(0.2*opt,100);
eta0=0.0001;
eta1=1.0;
y0=diff(xp,mesh.V(),del,eta0,n);
reduce(y0, sumOp<scalar>());
do
{
  eta5=(eta0+eta1)/2.0;
  y5=diff(xp,mesh.V(),del,eta5,n);
  reduce(y5, sumOp<scalar>());
  if(y0*y5<0)
  { 
     eta1=eta5;
  }
  else
  { 
     eta0=eta5;
     y0=y5;
  }
} while ((eta1-eta0)>0.0001);

for(i=0;i<n;i++)
{
   if(xp[i]<=eta5)
   {
     xh[i]=eta5*(Foam::exp(-del*(1-xp[i]/eta5))-(1-xp[i]/eta5)*Foam::exp(-del)); 
   }
   else
   {
     xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+(xp[i]-eta5)*Foam::exp(-del)/(1-eta5));
   }   
}

if (enableMaxFeatureConstraint)
{
  const scalar etaEroded = maxFeatureEtaEroded; // User-controlled from tuneOptParameters
  VEroded = 0.0;
  for(i=0;i<n;i++)
  {
     scalar xhErodedI;
     if(xp[i]<=etaEroded)
     {
       xhErodedI=etaEroded*(Foam::exp(-del*(1-xp[i]/etaEroded))-(1-xp[i]/etaEroded)*Foam::exp(-del)); 
     }
     else
     {
       xhErodedI=etaEroded+(1-etaEroded)*(1-Foam::exp(-del*(xp[i]-etaEroded)/(1-etaEroded))+(xp[i]-etaEroded)*Foam::exp(-del)/(1-etaEroded));
     }
     
     // Sensitivity chain rule drho_eroded / dxp
     scalar drhoEroded;
     if(xp[i]<=etaEroded)
     {
       drhoEroded=del*Foam::exp(-del*(1-xp[i]/etaEroded))+Foam::exp(-del); 
     }
     else
     {
       drhoEroded=del*Foam::exp(-del*(xp[i]-etaEroded)/(1-etaEroded))+Foam::exp(-del); 
     }
     
     gsenshMaxFeature[i] = drhoEroded;
     VEroded += xhErodedI * mesh.V()[i];
  }
  
  // Aggregate eroded volume globally across MPI ranks and normalize.
  reduce(VEroded, sumOp<scalar>());
  VEroded = VEroded / area;
  
  // Express as a constraint violation: positive = too-thick features present,
  // negative = solid features are already thin enough (no penalty needed).
  VEroded = VEroded - maxFeatureVolFracLimit;
}



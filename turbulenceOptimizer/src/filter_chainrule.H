fsenshMeanT=fsenshMeanT/gMax(mag(fsenshMeanT.primitiveField()));
gsenshPowerDiss=gsenshPowerDiss/gMax(mag(gsenshPowerDiss.primitiveField()));
gsenshVol.primitiveFieldRef() = 1.0;
//***********************************//Heaviside_filter
if(opt>1)
{
  for(i=0;i<n;i++)
  {
     if(xp[i]<=eta5)
     {
       drho[i]=del*Foam::exp(-del*(1-xp[i]/eta5))+Foam::exp(-del); 
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(xp[i]-eta5)/(1-eta5))+Foam::exp(-del); 
     } 
  }
  fsenshMeanT=fsenshMeanT*drho;
  gsenshPowerDiss=gsenshPowerDiss*drho;
  gsenshVol=gsenshVol*drho;
}
//***********************************//PDE filter
const auto fsensFilterLinearSolvePerformance =
solve(fvm::laplacian(fsensMeanT)-fvm::Sp(b,fsensMeanT)+fsenshMeanT*b);
recordLinearSolve
(
  fsensFilterLinearSolveStats,
  fsensFilterLinearSolvePerformance.initialResidual(),
  fsensFilterLinearSolvePerformance.finalResidual(),
  fsensFilterLinearSolvePerformance.nIterations()
);

const auto gsensPowerFilterLinearSolvePerformance =
solve(fvm::laplacian(gsensPowerDiss)-fvm::Sp(b,gsensPowerDiss)+gsenshPowerDiss*b);
recordLinearSolve
(
  gsensPowerFilterLinearSolveStats,
  gsensPowerFilterLinearSolvePerformance.initialResidual(),
  gsensPowerFilterLinearSolvePerformance.finalResidual(),
  gsensPowerFilterLinearSolvePerformance.nIterations()
);

const auto gsensVolFilterLinearSolvePerformance =
solve(fvm::laplacian(gsensVol)-fvm::Sp(b,gsensVol)+gsenshVol*b);
recordLinearSolve
(
  gsensVolFilterLinearSolveStats,
  gsensVolFilterLinearSolvePerformance.initialResidual(),
  gsensVolFilterLinearSolvePerformance.finalResidual(),
  gsensVolFilterLinearSolvePerformance.nIterations()
);

fsensMeanT.primitiveFieldRef()=fsensMeanT.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensPowerDiss.primitiveFieldRef()=gsensPowerDiss.primitiveFieldRef()*mesh.V()/gMax(mesh.V());
gsensVol.primitiveFieldRef()=gsensVol.primitiveFieldRef()*mesh.V()/gMax(mesh.V());

gsenshVol.primitiveFieldRef() = 1.0;
//***********************************//Heaviside_filter
if(opt>1)
{
  const scalar etaSafe = Foam::min(Foam::max(eta5, projectionEtaEps), scalar(1.0)-projectionEtaEps);
  const scalar etaDenominator = Foam::max(etaSafe, projectionEtaEps);
  const scalar oneMinusEtaDenominator = Foam::max(scalar(1.0)-etaSafe, projectionEtaEps);
  for(i=0;i<n;i++)
  {
     const scalar xpClamped = Foam::min(Foam::max(xp[i], scalar(0.0)), scalar(1.0));
     if(xpClamped<=etaSafe)
     {
       drho[i]=del*Foam::exp(-del*(1-xpClamped/etaDenominator))+Foam::exp(-del);
     }
     else
     {
       drho[i]=del*Foam::exp(-del*(xpClamped-etaSafe)/oneMinusEtaDenominator)+Foam::exp(-del);
     }
     if (!std::isfinite(static_cast<double>(drho[i])))
     {
       FatalErrorInFunction
         << "Non-finite projection derivative at cell " << i
         << " (xp=" << xp[i]
         << ", xpClamped=" << xpClamped
         << ", eta=" << etaSafe
         << ", del=" << del << ")"
         << exit(FatalError);
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

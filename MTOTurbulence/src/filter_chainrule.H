// Normalize sensitivities — guard against div-by-zero on iteration 1
const scalar maxFsens = Foam::max(gMax(mag(fsenshMeanT.primitiveField())), SMALL);
const scalar maxGsPower = Foam::max(gMax(mag(gsenshPowerDiss.primitiveField())), SMALL);
fsenshMeanT   = fsenshMeanT   / maxFsens;
gsenshPowerDiss = gsenshPowerDiss / maxGsPower;
gsenshVol.primitiveFieldRef() = 1.0;

// Always zero-initialize max-feature sensitivity so the Helmholtz solve
// receives a valid field even when the Heaviside is not yet active (opt == 1).
if (enableMaxFeatureConstraint)
{
    gsenshMaxFeature.primitiveFieldRef() = 0.0;
}

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
  } // End of for loop
  fsenshMeanT=fsenshMeanT*drho;
  gsenshPowerDiss=gsenshPowerDiss*drho;
  gsenshVol=gsenshVol*drho;

  if (enableMaxFeatureConstraint)
  {
      gsenshMaxFeature = gsenshMaxFeature * drho;
  }
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

if (enableMaxFeatureConstraint)
{
    const auto gsensMaxFeatureFilterSolvePerformance =
    solve(fvm::laplacian(gsensMaxFeature)-fvm::Sp(b,gsensMaxFeature)+gsenshMaxFeature*b);
    recordLinearSolve
    (
      gsensMaxFeatureFilterLinearSolveStats,
      gsensMaxFeatureFilterSolvePerformance.initialResidual(),
      gsensMaxFeatureFilterSolvePerformance.finalResidual(),
      gsensMaxFeatureFilterSolvePerformance.nIterations()
    );
}

// Scale sensitivities by local cell volume (normalized by max cell volume)
const scalar maxVol = Foam::max(gMax(mesh.V()), SMALL);
fsensMeanT.primitiveFieldRef()    = fsensMeanT.primitiveFieldRef()    * mesh.V() / maxVol;
gsensPowerDiss.primitiveFieldRef()= gsensPowerDiss.primitiveFieldRef()* mesh.V() / maxVol;
gsensVol.primitiveFieldRef()      = gsensVol.primitiveFieldRef()      * mesh.V() / maxVol;
if (enableMaxFeatureConstraint)
{
    gsensMaxFeature.primitiveFieldRef() = gsensMaxFeature.primitiveFieldRef() * mesh.V() / maxVol;
}



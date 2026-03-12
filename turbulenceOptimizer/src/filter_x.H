
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
{
  scalar maxAbsXStepForProjectionLocal(0.0);
  for(i=0;i<n;i++)
  {
    maxAbsXStepForProjectionLocal =
      Foam::max(maxAbsXStepForProjectionLocal, Foam::mag(x[i] - xOld[i]));
  }
  reduce(maxAbsXStepForProjectionLocal, maxOp<scalar>());

  const scalar powerConstraintLimitForProjection =
    Foam::max(Foam::max(PowerDissRelax-opt*0.2, PowerDissMax), scalar(SMALL));
  const scalar powerConstraintViolationForProjection =
    Foam::max(powerConstraintNormalized/powerConstraintLimitForProjection - 1.0, scalar(0.0));
  scalar objectiveRelativeDeltaForProjection(GREAT);
  if (hasPreviousObjectiveValue)
  {
    objectiveRelativeDeltaForProjection =
      Foam::mag(objectiveNormalized-previousObjectiveValue)
     /Foam::max(Foam::mag(objectiveNormalized), scalar(SMALL));
  }

  const bool projectionUpdateIteration =
    (opt == 1) || ((opt % projectionUpdateInterval) == 0);
  bool projectionUpdateStable =
    objectiveRelativeDeltaForProjection <= projectionObjectiveTolForBeta
    &&
    maxAbsXStepForProjectionLocal <= projectionStepTolForBeta;
  if (projectionRequireConstraintForBeta)
  {
    projectionUpdateStable =
      projectionUpdateStable
      &&
      powerConstraintViolationForProjection <= projectionConstraintTolForBeta;
  }

  if (opt == 1)
  {
    del = projectionBetaInitial;
  }
  else if (projectionUpdateIteration && projectionUpdateStable)
  {
    del = Foam::min(del + projectionBetaRate, projectionBetaMax);
  }
  else
  {
    del = Foam::min(Foam::max(del, projectionBetaInitial), projectionBetaMax);
  }
}
eta0=projectionEtaLower;
eta1=projectionEtaUpper;
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
} while ((eta1-eta0)>projectionEtaTolerance);

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
   xh[i]=Foam::min(Foam::max(xh[i],rhoMin),scalar(1.0));
}

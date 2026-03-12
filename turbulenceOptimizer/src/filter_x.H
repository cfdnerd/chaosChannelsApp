
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
if (adaptiveHeavisideProjection)
{
  scalar maxAbsXStepForProjectionLocal(0.0);
  for(i=0;i<n;i++)
  {
    maxAbsXStepForProjectionLocal =
      Foam::max(maxAbsXStepForProjectionLocal, Foam::mag(x[i] - xOld[i]));
  }
  reduce(maxAbsXStepForProjectionLocal, maxOp<scalar>());

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
else
{
  // Conservative deterministic continuation: multiplicative beta updates.
  const label betaUpdateCount = Foam::max((opt-1)/projectionUpdateInterval, label(0));
  del =
    Foam::min
    (
      projectionBetaInitial*Foam::pow(projectionBetaGrowthFactor, betaUpdateCount),
      projectionBetaMax
    );
}
eta0=projectionEtaLower;
eta1=projectionEtaUpper;
y0=diff(xp,mesh.V(),del,eta0,n);
reduce(y0, sumOp<scalar>());
scalar y1=diff(xp,mesh.V(),del,eta1,n);
reduce(y1, sumOp<scalar>());
const bool etaBracketed = (y0*y1 <= 0.0);
const label projectionEtaMaxIterations(100);
label projectionEtaIterations(0);
if (!etaBracketed)
{
  eta5=0.5*(eta0+eta1);
  y5=y0;
}
else
{
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
    projectionEtaIterations++;
  } while ((eta1-eta0)>projectionEtaTolerance && projectionEtaIterations < projectionEtaMaxIterations);
}

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

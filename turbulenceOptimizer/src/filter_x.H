
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
eta0=Foam::min(Foam::max(projectionEtaLower, projectionEtaEps), scalar(1.0)-projectionEtaEps);
eta1=Foam::min(Foam::max(projectionEtaUpper, projectionEtaEps), scalar(1.0)-projectionEtaEps);
if (eta1 <= eta0)
{
  eta0 = projectionEtaEps;
  eta1 = scalar(1.0)-projectionEtaEps;
}
y0=diff(xp,mesh.V(),del,eta0,n);
reduce(y0, sumOp<scalar>());
scalar y1=diff(xp,mesh.V(),del,eta1,n);
reduce(y1, sumOp<scalar>());
if
(
  !std::isfinite(static_cast<double>(del))
 || !std::isfinite(static_cast<double>(y0))
 || !std::isfinite(static_cast<double>(y1))
)
{
  FatalErrorInFunction
    << "Non-finite projection state detected:"
    << " del=" << del
    << " y0=" << y0
    << " y1=" << y1
    << " eta0=" << eta0
    << " eta1=" << eta1
    << exit(FatalError);
}
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
    if (!std::isfinite(static_cast<double>(y5)))
    {
      FatalErrorInFunction
        << "Non-finite projection root function value at eta5=" << eta5
        << ", del=" << del
        << exit(FatalError);
    }
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

eta5 = Foam::min(Foam::max(eta5, projectionEtaEps), scalar(1.0)-projectionEtaEps);
const scalar etaDenominator = Foam::max(eta5, projectionEtaEps);
const scalar oneMinusEtaDenominator = Foam::max(scalar(1.0)-eta5, projectionEtaEps);
for(i=0;i<n;i++)
{
   const scalar xpClamped = Foam::min(Foam::max(xp[i], scalar(0.0)), scalar(1.0));
   if(xpClamped<=eta5)
   {
     const scalar ratio = xpClamped/etaDenominator;
     xh[i]=eta5*(Foam::exp(-del*(1-ratio))-(1-ratio)*Foam::exp(-del));
   }
   else
   {
     const scalar ratio = (xpClamped-eta5)/oneMinusEtaDenominator;
     xh[i]=eta5+(1-eta5)*(1-Foam::exp(-del*ratio)+ratio*Foam::exp(-del));
   }
   if (!std::isfinite(static_cast<double>(xh[i])))
   {
     FatalErrorInFunction
       << "Non-finite xh in projection at cell " << i
       << " (xp=" << xp[i]
       << ", xpClamped=" << xpClamped
       << ", eta5=" << eta5
       << ", del=" << del << ")"
       << exit(FatalError);
   }
   xh[i]=Foam::min(Foam::max(xh[i],rhoMin),scalar(1.0));
}

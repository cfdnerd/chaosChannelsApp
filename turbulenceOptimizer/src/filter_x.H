
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
const bool projectionContinuationActive = (opt >= projectionStartIteration);
const bool projectionUpdateIteration =
  projectionContinuationActive
  &&
  (
    (opt == projectionStartIteration)
    || (((opt - projectionStartIteration) % projectionUpdateInterval) == 0)
  );
const scalar projectionConstraintViolation =
  Foam::max(weightedPowerConstraint, scalar(0.0));
const bool allowProjectionContinuation =
  projectionConstraintViolation <= projectionConstraintViolationTolForBeta;
const scalar projectionSoftGrowthFactor =
  scalar(1.0) + scalar(0.15)*(projectionBetaGrowthFactor - scalar(1.0));
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

  bool projectionUpdateStable =
    objectiveRelativeDeltaForProjection <= projectionObjectiveTolForBeta
    &&
    maxAbsXStepForProjectionLocal <= projectionStepTolForBeta;

  if (opt == 1)
  {
    del = projectionBetaInitial;
  }
  else if (projectionUpdateIteration && projectionUpdateStable && allowProjectionContinuation)
  {
    del = Foam::min(del + projectionBetaRate, projectionBetaMax);
  }
  else if (projectionUpdateIteration && !allowProjectionContinuation)
  {
    // Under constraint stress, grow slowly instead of freezing at beta=1.
    del = Foam::min(Foam::max(del*projectionSoftGrowthFactor, projectionBetaInitial), projectionBetaMax);
  }
  else
  {
    del = Foam::min(Foam::max(del, projectionBetaInitial), projectionBetaMax);
  }
}
else
{
  // Deterministic continuation: keep sharpening active in non-adaptive mode.
  if (opt == 1)
  {
    del = projectionBetaInitial;
  }
  else if (projectionUpdateIteration && allowProjectionContinuation)
  {
    del = Foam::min(del*projectionBetaGrowthFactor, projectionBetaMax);
  }
  else if (projectionUpdateIteration && !allowProjectionContinuation)
  {
    // Hold projection sharpness when power-drop constraint is stressed.
    del = Foam::min(Foam::max(del, projectionBetaInitial), projectionBetaMax);
  }
  else
  {
    del = Foam::min(Foam::max(del, projectionBetaInitial), projectionBetaMax);
  }
}
if (projectionNumericalSafeguards)
{
  eta0=Foam::min(Foam::max(projectionEtaLower, projectionEtaEps), scalar(1.0)-projectionEtaEps);
  eta1=Foam::min(Foam::max(projectionEtaUpper, projectionEtaEps), scalar(1.0)-projectionEtaEps);
  if (eta1 <= eta0)
  {
    eta0 = projectionEtaEps;
    eta1 = scalar(1.0)-projectionEtaEps;
  }
}
else
{
  eta0 = projectionEtaLower;
  eta1 = projectionEtaUpper;
}
y0=diff(xp,mesh.V(),del,eta0,n,projectionNumericalSafeguards);
reduce(y0, sumOp<scalar>());
scalar y1=diff(xp,mesh.V(),del,eta1,n,projectionNumericalSafeguards);
reduce(y1, sumOp<scalar>());
if
(
  projectionNumericalSafeguards
  &&
  (
    !std::isfinite(static_cast<double>(del))
   || !std::isfinite(static_cast<double>(y0))
   || !std::isfinite(static_cast<double>(y1))
  )
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
    y5=diff(xp,mesh.V(),del,eta5,n,projectionNumericalSafeguards);
    reduce(y5, sumOp<scalar>());
    if (projectionNumericalSafeguards && !std::isfinite(static_cast<double>(y5)))
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

const scalar etaProjection =
  projectionNumericalSafeguards
  ? Foam::min(Foam::max(eta5, projectionEtaEps), scalar(1.0)-projectionEtaEps)
  : eta5;
const scalar etaDenominator =
  projectionNumericalSafeguards ? Foam::max(etaProjection, projectionEtaEps) : etaProjection;
const scalar oneMinusEtaDenominator =
  projectionNumericalSafeguards
  ? Foam::max(scalar(1.0)-etaProjection, projectionEtaEps)
  : (scalar(1.0)-etaProjection);
for(i=0;i<n;i++)
{
   const scalar xpProjection =
     projectionNumericalSafeguards
     ? Foam::min(Foam::max(xp[i], scalar(0.0)), scalar(1.0))
     : xp[i];
   if(xpProjection<=etaProjection)
   {
     const scalar ratio = xpProjection/etaDenominator;
     xh[i]=etaProjection*(Foam::exp(-del*(1-ratio))-(1-ratio)*Foam::exp(-del));
   }
   else
   {
     const scalar ratio = (xpProjection-etaProjection)/oneMinusEtaDenominator;
     xh[i]=etaProjection+(1-etaProjection)*(1-Foam::exp(-del*ratio)+ratio*Foam::exp(-del));
   }
   if (projectionNumericalSafeguards && !std::isfinite(static_cast<double>(xh[i])))
   {
     FatalErrorInFunction
       << "Non-finite xh in projection at cell " << i
       << " (xp=" << xp[i]
       << ", xpProjection=" << xpProjection
       << ", etaProjection=" << etaProjection
       << ", del=" << del << ")"
       << exit(FatalError);
   }
   xh[i]=Foam::min(Foam::max(xh[i],rhoMin),scalar(1.0));
}
eta5 = etaProjection;

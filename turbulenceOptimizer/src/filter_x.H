
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
scalar grayVolumeForProjectionLocal(0.0);
for(i=0;i<n;i++)
{
  const scalar xhValue = xhOld[i];
  if (xhValue > 0.05 && xhValue < 0.95)
  {
    grayVolumeForProjectionLocal += mesh.V()[i];
  }
}
reduce(grayVolumeForProjectionLocal, sumOp<scalar>());
projectionGrayFraction =
  grayVolumeForProjectionLocal/Foam::max(area, scalar(SMALL));
const scalar projectionConstraintViolation =
  Foam::max(weightedPowerConstraint, scalar(0.0));
const bool allowProjectionContinuation =
  projectionConstraintViolation <= projectionConstraintViolationTolForBeta;
scalar maxAbsXStepForProjectionLocal(0.0);
scalar l2XStepForProjectionLocal(0.0);
for(i=0;i<n;i++)
{
  const scalar xStep = x[i] - xOld[i];
  maxAbsXStepForProjectionLocal =
    Foam::max(maxAbsXStepForProjectionLocal, Foam::mag(xStep));
  l2XStepForProjectionLocal += xStep*xStep;
}
reduce(maxAbsXStepForProjectionLocal, maxOp<scalar>());
reduce(l2XStepForProjectionLocal, sumOp<scalar>());
const scalar normalizedL2XStepForProjection =
  std::sqrt(l2XStepForProjectionLocal)
 /std::sqrt(Foam::max(activeDesignDOF, scalar(1)));

scalar objectiveRelativeDeltaForProjection(GREAT);
if (hasPreviousObjectiveValue)
{
  objectiveRelativeDeltaForProjection =
    Foam::mag(objectiveNormalized-previousObjectiveValue)
   /Foam::max(Foam::mag(objectiveNormalized), scalar(SMALL));
}

const bool projectionUpdateStable =
  objectiveRelativeDeltaForProjection <= projectionObjectiveTolForBeta
  &&
  maxAbsXStepForProjectionLocal <= projectionStepTolForBeta;
const bool projectionGrayStallStepSmall =
  maxAbsXStepForProjectionLocal <= grayStallStepTol
  ||
  normalizedL2XStepForProjection <= grayStallStepTol;
const bool projectionGrayStallObserved =
  projectionContinuationActive
  &&
  projectionGrayFraction >= grayStallThreshold
  &&
  objectiveRelativeDeltaForProjection <= grayStallObjectiveTol
  &&
  projectionGrayStallStepSmall;
projectionGrayStallCounter =
  projectionGrayStallObserved ? projectionGrayStallCounter + 1 : 0;
projectionGrayStallRecovery =
  projectionGrayStallCounter >= grayStallWindow;
projectionLowGrayFreezeActive =
  projectionContinuationActive
  &&
  opt >= projectionLowGrayFreezeStartIteration
  &&
  del >= projectionLowGrayFreezeMinBeta
  &&
  projectionGrayFraction <= projectionLowGrayFreezeThreshold
  &&
  objectiveRelativeDeltaForProjection <= projectionLowGrayFreezeObjectiveTol
  &&
  normalizedL2XStepForProjection <= projectionLowGrayFreezeStepTol;
const scalar projectionConstraintViolationBand =
  Foam::max
  (
    projectionConstraintViolationSoftTol - projectionConstraintViolationTolForBeta,
    scalar(SMALL)
  );
const scalar projectionConstraintStress =
  allowProjectionContinuation
  ? scalar(0.0)
  : Foam::min
    (
      Foam::max
      (
        projectionConstraintViolation - projectionConstraintViolationTolForBeta,
        scalar(0.0)
      )/projectionConstraintViolationBand,
      scalar(1.0)
    );
const scalar projectionGrowthFactorUnderStress =
  projectionBetaGrowthFactor
  - projectionConstraintStress
   *(projectionBetaGrowthFactor - projectionBetaGrowthFactorMin);
const scalar projectionStressStepScale =
  projectionBetaGrowthFactor > scalar(1.0)
  ? Foam::max
    (
      (projectionGrowthFactorUnderStress - scalar(1.0))
     /(projectionBetaGrowthFactor - scalar(1.0)),
      scalar(0.0)
    )
  : scalar(0.0);
const scalar delPrevious(del);
if (adaptiveHeavisideProjection)
{
  if (opt == 1)
  {
    del = projectionBetaInitial;
  }
  else if (projectionLowGrayFreezeActive)
  {
    del = Foam::min(Foam::max(del, projectionBetaInitial), projectionBetaMax);
  }
  else if (projectionUpdateIteration)
  {
    scalar betaRateStep(projectionBetaRate);
    if (!allowProjectionContinuation)
    {
      betaRateStep *= Foam::max(projectionStressStepScale, scalar(0.1));
    }
    if (!projectionUpdateStable)
    {
      betaRateStep *= projectionUnstableGrowthScale;
    }
    del = Foam::min(del + betaRateStep, projectionBetaMax);
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
  else if (projectionLowGrayFreezeActive)
  {
    del = Foam::min(Foam::max(del, projectionBetaInitial), projectionBetaMax);
  }
  else if (projectionUpdateIteration)
  {
    scalar betaGrowthFactor =
      allowProjectionContinuation
      ? projectionBetaGrowthFactor
      : projectionGrowthFactorUnderStress;
    if (!projectionUpdateStable)
    {
      betaGrowthFactor =
        scalar(1.0)
       +projectionUnstableGrowthScale*(betaGrowthFactor - scalar(1.0));
    }
    del = Foam::min(del*betaGrowthFactor, projectionBetaMax);
  }
  else
  {
    del = Foam::min(Foam::max(del, projectionBetaInitial), projectionBetaMax);
  }
}
if (projectionUpdateIteration && projectionGrayStallRecovery)
{
  del =
    Foam::min
    (
      Foam::max(del*grayStallBetaGrowthFactor, del + projectionBetaRate),
      projectionBetaMax
    );
}
if (projectionLowGrayFreezeActive)
{
  del = delPrevious;
}
projectionGrowthFactorApplied =
  delPrevious > SMALL ? del/delPrevious : scalar(1.0);
eta0 = projectionEtaLower;
eta1 = projectionEtaUpper;
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

const scalar etaProjection = eta5;
const scalar etaDenominator = etaProjection;
const scalar oneMinusEtaDenominator = scalar(1.0)-etaProjection;
for(i=0;i<n;i++)
{
   const scalar xpProjection = xp[i];
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
   xh[i]=Foam::min(Foam::max(xh[i],rhoMin),scalar(1.0));
}
eta5 = etaProjection;

//***********************************//update mu and lambda
opt++;
forAll(x, i)
{
   x[i]=Foam::min(Foam::max(x[i],rhoMin),scalar(1.0));
}
if(solid_area)
{
   forAll(cells_solid, i)
   {
      x[cells_solid[i]]=rhoMin;
      xh[cells_solid[i]]=rhoMin;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      x[cells_fluid[i]]=1.0;
      xh[cells_fluid[i]]=1.0;
   }
}
if(test_area)
{
   forAll(cells_test, i)
   {
      x[cells_test[i]]=1.0;
      xh[cells_test[i]]=1.0;
   }
}
forAll(xh, i)
{
   xh[i]=Foam::min(Foam::max(xh[i],rhoMin),scalar(1.0));
}
xh.correctBoundaryConditions();

const scalar previousConstraintViolation =
  Foam::max(static_cast<scalar>(gx[0]), scalar(0.0));
const bool continuationFreeze =
  previousConstraintViolation > continuationConstraintViolationFreeze;
const scalar continuationBlend =
  continuationFreeze ? continuationBlendSlow : continuationBlendFast;

scalar alphaTargetValue(alphaMax.value());
if(opt<=alphaRampIterations)
{
  alphaTargetValue = alphamax.value()*(scalar(opt)/alphaRampDivisor+1.0);
}
if(opt>alphaGrowthStartIteration)
{
  alphaTargetValue *= alphaGrowthFactor;
}
alphaTargetValue =
  Foam::min
  (
    alphaTargetValue,
    Foam::min(alphaMAX.value(), alphamax.value()*alphaContinuationCapMultiplier)
  );
alphaMax.value() =
  alphaMax.value() + continuationBlend*(alphaTargetValue-alphaMax.value());

const scalar quPrevious(qu);
const scalar quContinuationBlend = continuationBlendFast;
scalar quTarget(quStage1Value);
if(opt<quStage1Iteration)
{
  quTarget = quStage1Value;
}
else if(opt<quStage2Iteration)
{
  const scalar t =
    scalar(opt-quStage1Iteration)/scalar(quStage2Iteration-quStage1Iteration);
  quTarget = quStage1Value + (quStage2Value-quStage1Value)*t;
}
else if(opt<quStage3Iteration)
{
  const scalar t =
    scalar(opt-quStage2Iteration)/scalar(quStage3Iteration-quStage2Iteration);
  quTarget = quStage2Value + (quStage3Value-quStage2Value)*t;
}
else
{
  const label quStage4RampEnd =
    quStage3Iteration + Foam::max(quStage3Iteration-quStage2Iteration, label(1));
  if (opt<quStage4RampEnd)
  {
    const scalar t =
      scalar(opt-quStage3Iteration)/scalar(quStage4RampEnd-quStage3Iteration);
    quTarget = quStage3Value + (quStage4Value-quStage3Value)*t;
  }
  else
  {
    quTarget = quStage4Value;
  }
}
qu = quPrevious + quContinuationBlend*(quTarget-quPrevious);
qu = Foam::min(qu,quMax);

scalar filterSchedule(1.0);
if (filterRContinuationEnd > filterRContinuationStart)
{
  filterSchedule =
    Foam::min
    (
      Foam::max
      (
        scalar(opt-filterRContinuationStart)
       /scalar(filterRContinuationEnd-filterRContinuationStart),
        scalar(0.0)
      ),
      scalar(1.0)
    );
}
const scalar currentFilterR = filterRStart + (filterREnd-filterRStart)*filterSchedule;
b.value()=1.0/(currentFilterR*len/3.464)/(currentFilterR*len/3.464);

alpha=alphaMax*qu*(1-xh)/(qu+xh);
DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoCp;

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
      xh[cells_solid[i]]=rhoMin;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      xh[cells_fluid[i]]=1.0;
   }
}
forAll(xh, i)
{
   xh[i]=Foam::min(Foam::max(xh[i],rhoMin),scalar(1.0));
}
xh.correctBoundaryConditions();

if(opt<=alphaRampIterations)
{
  alphaMax=alphamax*(opt/alphaRampDivisor+1.0);
}
if(opt>alphaGrowthStartIteration)
{
  alphaMax=alphaMax*alphaGrowthFactor;
}
alphaMax.value()=
  Foam::min
  (
    alphaMax.value(),
    Foam::min(alphaMAX.value(), alphamax.value()*alphaContinuationCapMultiplier)
  );
if(opt<quStage1Iteration)
{
  qu=quStage1Value;
}
else if(opt<quStage2Iteration)
{
  qu=quStage2Value;
}
else if(opt<quStage3Iteration)
{
  qu=quStage3Value;
}
else
{
  qu=quStage4Value;
}
qu=Foam::min(qu,quMax);

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

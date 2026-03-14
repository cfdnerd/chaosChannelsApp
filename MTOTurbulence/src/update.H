//***********************************//update mu and lambda
opt++;
if(solid_area)
{
   forAll(cells_solid, i)
   {
      xh[cells_solid[i]]=0.0;
      x[cells_solid[i]]=0.0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, i)
   {
      xh[cells_fluid[i]]=1.0;
      x[cells_fluid[i]]=1.0;
   }
}
if(test_area)
{
   forAll(cells_test, i)
   {
      xh[cells_test[i]]=1.0;
      x[cells_test[i]]=1.0;
   }
}
xh.correctBoundaryConditions();

if(opt<=63)
{
  alphaMax=alphamax*(opt/7.0+1.0);
}
if(opt>100)
{
  alphaMax=alphaMax*1.05;
}
alphaMax.value()=Foam::min(alphaMax.value(),alphaMAX.value());
if(opt>=80)
{
  qu=0.005+(opt-80)*1e-4;
}
qu=Foam::min(qu,0.01);
alpha=alphaMax*qu*(1-xh)/(qu+xh);
DT=(ks+(kf-ks)*xh*(1+qu)/(qu+xh))/rhoCp;

// Sanitize design fields immediately after update to prevent NaN/Inf propagation
forAll(xh, i)
{
    if (!std::isfinite(xh[i]) || xh[i] > 1.0) xh[i] = 1.0;
    if (xh[i] < 0.0) xh[i] = 0.0;
}
forAll(xp, i)
{
    if (!std::isfinite(xp[i]) || xp[i] > 1.0) xp[i] = 1.0;
    if (xp[i] < 0.0) xp[i] = 0.0;
}
xh.correctBoundaryConditions();
xp.correctBoundaryConditions();

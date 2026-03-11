
Info<<"sensitivity analysis"<<endl;

std::vector<scalar> xOld(n);
std::vector<scalar> xhOld(n);
for(i=0;i<n;i++)
{
 xOld[i]=x[i];
 xhOld[i]=xh[i];
}

fsenshMeanT=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ub)+b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tb))/rhoCp;
if(DDHS)//design dependent heat source
{
 fsenshMeanT += b1*Tb*thermalSource;
}
gsenshPowerDiss=-alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua);


#include "filter_chainrule.H"

const scalar weightedPowerConstraint =
  optPowerConstraintWeight*(PowerDiss/Foam::max(PowerDissRelax-opt*0.2,PowerDissMax)-1.0);
const scalar weightedVolumeConstraint = optVolumeConstraintWeight*V;
gx[0]=weightedPowerConstraint;
gx[1]=weightedVolumeConstraint;

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=fsensMeanT[i]/N;//sensitivities of objective function
 dgdx[0][i]=optPowerConstraintWeight*(gsensPowerDiss[i]/N);//sensitivities of constraint function
 dgdx[1][i]=optVolumeConstraintWeight*(gsensVol[i]/N);//sensitivities of constraint function
}


Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=xmma[i];
}

#include "filter_x.H"
#include "debugOptimizer.H"
#include "solverConvergenceWrite.H"

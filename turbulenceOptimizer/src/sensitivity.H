Info<<"sensitivity analysis"<<endl;

volScalarField dAlphaDx
(
  -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*alphaResistanceScale
);
volScalarField fsenshCurrent
(
  dAlphaDx*(U&Ub)
 +b1*(kf-ks)*(1+qu)*qu/((qu+xh)*(qu+xh))*(fvc::grad(T)&fvc::grad(Tb))/rhoCp
);
if(DDHS)//design dependent heat source
{
 fsenshCurrent += b1*Tb*thermalSource;
}
volScalarField gsenshPowerCurrent(dAlphaDx*(U&Ua));

std::vector<scalar> xOld(n);
std::vector<scalar> xhOld(n);
for(i=0;i<n;i++)
{
 xOld[i]=x[i];
 xhOld[i]=xh[i];
}

fsenshMeanT = fsenshCurrent;
gsenshPowerDiss = gsenshPowerCurrent;
#include "filter_chainrule.H"

const scalar tighteningIterations =
  scalar(Foam::max(opt - powerConstraintTighteningStartIteration, label(0)));
const scalar powerConstraintDenominator =
  Foam::max(PowerDissRelax - powerConstraintTighteningRate*tighteningIterations, PowerDissMax);
const scalar weightedPowerConstraint =
  optPowerConstraintWeight*(powerConstraintNormalized/powerConstraintDenominator-1.0);
const scalar weightedVolumeConstraint = optVolumeConstraintWeight*V;
gx[0]=weightedPowerConstraint;
gx[1]=weightedVolumeConstraint;

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 dfdx[i]=fsensMeanT[i]/N;//sensitivities of objective function
 dgdx[0][i]=optPowerConstraintWeight*(gsensPowerDiss[i]/N)/powerConstraintDenominator;//sensitivities of constraint function
 dgdx[1][i]=optVolumeConstraintWeight*(gsensVol[i]/N);//sensitivities of constraint function
}

Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
for(i=0;i<n;i++)
{
   x[i]=Foam::min(Foam::max(xmma[i],rhoMin),scalar(1.0));
}

#include "filter_x.H"
if (Pstream::master())
{
  Info<< "post-MMA: entering debug diagnostics" << endl;
}
#include "debugOptimizer.H"
if (Pstream::master())
{
  Info<< "post-debug: entering solver convergence logging" << endl;
}
#include "solverConvergenceWrite.H"
if (Pstream::master())
{
  Info<< "post-logging: optimizer iteration diagnostics complete" << endl;
}

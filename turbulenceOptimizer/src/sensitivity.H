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
scalar powerConstraintReference(1.0);
if (normalizeScenarioObjectives)
{
 const scalar objectiveReference =
   objectiveIterationOneReferenceCaptured
   ? objectiveIterationOneReference
   : Foam::max(Foam::mag(MeanTPhysical), scalar(SMALL));
 powerConstraintReference =
   powerIterationOneReferenceCaptured
   ? powerIterationOneReference
   : Foam::max(Foam::mag(PowerDissPhysical), scalar(SMALL));
 fsenshMeanT = fsenshMeanT/Foam::max(objectiveReference, scalar(SMALL));
 gsenshPowerDiss = gsenshPowerDiss/Foam::max(powerConstraintReference, scalar(SMALL));
}
#include "filter_chainrule.H"

const scalar tighteningIterations =
  scalar(Foam::max(opt - powerConstraintTighteningStartIteration, label(0)));
scalar powerConstraintRelaxLimit(PowerDissRelax);
scalar powerConstraintHardLimit(PowerDissMax);
scalar powerConstraintTighteningRateLocal(powerConstraintTighteningRate);
if (normalizeScenarioObjectives)
{
  const scalar normalizationRef = Foam::max(powerConstraintReference, scalar(SMALL));
  powerConstraintRelaxLimit /= normalizationRef;
  powerConstraintHardLimit /= normalizationRef;
  powerConstraintTighteningRateLocal /= normalizationRef;
}
const scalar powerConstraintDenominatorRaw =
  Foam::max
  (
    powerConstraintRelaxLimit - powerConstraintTighteningRateLocal*tighteningIterations,
    powerConstraintHardLimit
  );
const scalar powerConstraintDenominator =
  Foam::max(powerConstraintDenominatorRaw, scalar(SMALL));
const scalar weightedPowerConstraint =
  optPowerConstraintWeight*(powerConstraintNormalized/powerConstraintDenominator-1.0);
const scalar weightedVolumeConstraint = optVolumeConstraintWeight*V;
gx[0]=weightedPowerConstraint;
gx[1]=weightedVolumeConstraint;

for(i=0;i<n;i++)
{
 xmma[i] = x[i];
 if (passiveCellMask[i])
 {
   dfdx[i]=0.0;
   dgdx[0][i]=0.0;
   dgdx[1][i]=0.0;
 }
 else
 {
   dfdx[i]=fsensMeanT[i]/activeDesignDOF;//sensitivities of objective function
   dgdx[0][i]=optPowerConstraintWeight*(gsensPowerDiss[i]/activeDesignDOF)/powerConstraintDenominator;//sensitivities of constraint function
   dgdx[1][i]=optVolumeConstraintWeight*(gsensVol[i]/activeDesignDOF);//sensitivities of constraint function
 }
}

Info<< "\n run MMA \n" << endl;
mma.MMAsolver(xmma, dfdx, gx, dgdx);
const scalar mmaPowerConstraintViolation =
  Foam::max(weightedPowerConstraint, scalar(0.0));
const scalar mmaDampingViolation =
  Foam::max(mmaPowerConstraintViolation - mmaConstraintDampingActivation, scalar(0.0));
const scalar mmaStepDamping =
  Foam::max
  (
    mmaConstraintDampingMin,
    scalar(1.0)
   /(
      scalar(1.0) + mmaConstraintDampingGain*mmaDampingViolation
    )
  );
for(i=0;i<n;i++)
{
   const scalar xCandidate =
     Foam::min(Foam::max(static_cast<scalar>(xmma[i]),rhoMin),scalar(1.0));
   x[i]=xOld[i]+mmaStepDamping*(xCandidate-xOld[i]);
   x[i]=Foam::min(Foam::max(x[i],rhoMin),scalar(1.0));
}
if(solid_area)
{
   forAll(cells_solid, iCell)
   {
      x[cells_solid[iCell]]=rhoMin;
      xh[cells_solid[iCell]]=rhoMin;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, iCell)
   {
      x[cells_fluid[iCell]]=1.0;
      xh[cells_fluid[iCell]]=1.0;
   }
}
if(test_area)
{
   forAll(cells_test, iCell)
   {
      x[cells_test[iCell]]=1.0;
      xh[cells_test[iCell]]=1.0;
   }
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

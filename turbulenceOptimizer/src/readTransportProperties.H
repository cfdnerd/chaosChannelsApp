Info<< "Reading TransportProperties\n" << endl;
IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);
//fluid properties
dimensionedScalar nu(transportProperties.lookup("nu"));//Viscosity of fluid
dimensionedScalar alphaMax(transportProperties.lookup("alphaMax"));
dimensionedScalar alphaMAX(transportProperties.lookup("alphaMAX"));
dimensionedScalar alphamax(transportProperties.lookup("alphaMax"));
scalar qu(quInitial);
volScalarField alpha(alphaMax*qu*(1-xh)/(qu+xh));
#include "createMRF.H"
Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);
#include "createPhi.H"
label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());
singlePhaseTransportModel laminarTransport(U, phi);
autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);
turbulence->validate();

label primalFlowIterations(60);
label adjointFlowIterations(40);
label turbulenceCorrectorInterval(1);
Switch useTurbulentResistance(true);
Switch frozenTurbulenceAdjoint(true);
Switch normalizeScenarioObjectives(true);
scalar turbulenceResistanceFactor(1.0);
scalar turbulenceResistanceCap(100.0);
scalar turbulenceStatePenaltyFactor(1.0);
scalar turbulenceIntensity(0.05);
scalar turbulenceLengthScale(1e-3);
scalar kMinBound(1e-10);
scalar epsilonMinBound(1e-10);

label inletPatchID(mesh.boundaryMesh().findPatchID("inlet"));
if (inletPatchID < 0)
{
  FatalErrorInFunction
    << "Missing required patch 'inlet' for multi-scenario setup."
    << exit(FatalError);
}
vector inletVelocityDirection(0, -1, 0);
scalar defaultInletSpeed(1.0);
if (U.boundaryField()[inletPatchID].size() > 0)
{
  const vector inletVelocitySample(U.boundaryField()[inletPatchID][0]);
  const scalar inletSpeedSample(mag(inletVelocitySample));
  if (inletSpeedSample > SMALL)
  {
    defaultInletSpeed = inletSpeedSample;
    inletVelocityDirection = inletVelocitySample/inletSpeedSample;
  }
}

scalarList scenarioVelocities(1, defaultInletSpeed);
scalarList scenarioWeights(1, 1.0);
if (tuneOptParameters.found("turbulenceOptimization"))
{
  const dictionary& turbulenceOptimization = tuneOptParameters.subDict("turbulenceOptimization");
  primalFlowIterations =
    turbulenceOptimization.lookupOrDefault<label>("primalFlowIterations", primalFlowIterations);
  adjointFlowIterations =
    turbulenceOptimization.lookupOrDefault<label>("adjointFlowIterations", adjointFlowIterations);
  turbulenceCorrectorInterval =
    turbulenceOptimization.lookupOrDefault<label>("turbulenceCorrectorInterval", turbulenceCorrectorInterval);
  useTurbulentResistance =
    turbulenceOptimization.lookupOrDefault<Switch>("useTurbulentResistance", useTurbulentResistance);
  frozenTurbulenceAdjoint =
    turbulenceOptimization.lookupOrDefault<Switch>("frozenTurbulenceAdjoint", frozenTurbulenceAdjoint);
  normalizeScenarioObjectives =
    turbulenceOptimization.lookupOrDefault<Switch>("normalizeScenarioObjectives", normalizeScenarioObjectives);
  turbulenceResistanceFactor =
    turbulenceOptimization.lookupOrDefault<scalar>("turbulenceResistanceFactor", turbulenceResistanceFactor);
  turbulenceResistanceCap =
    turbulenceOptimization.lookupOrDefault<scalar>("turbulenceResistanceCap", turbulenceResistanceCap);
  turbulenceStatePenaltyFactor =
    turbulenceOptimization.lookupOrDefault<scalar>("turbulenceStatePenaltyFactor", turbulenceStatePenaltyFactor);
  turbulenceIntensity =
    turbulenceOptimization.lookupOrDefault<scalar>("turbulenceIntensity", turbulenceIntensity);
  turbulenceLengthScale =
    turbulenceOptimization.lookupOrDefault<scalar>("turbulenceLengthScale", turbulenceLengthScale);
  kMinBound = turbulenceOptimization.lookupOrDefault<scalar>("kMinBound", kMinBound);
  epsilonMinBound = turbulenceOptimization.lookupOrDefault<scalar>("epsilonMinBound", epsilonMinBound);
  if (turbulenceOptimization.found("scenarioVelocities"))
  {
    turbulenceOptimization.lookup("scenarioVelocities") >> scenarioVelocities;
  }
  if (turbulenceOptimization.found("scenarioWeights"))
  {
    turbulenceOptimization.lookup("scenarioWeights") >> scenarioWeights;
  }
}

primalFlowIterations = Foam::max(primalFlowIterations, label(1));
adjointFlowIterations = Foam::max(adjointFlowIterations, label(1));
turbulenceCorrectorInterval = Foam::max(turbulenceCorrectorInterval, label(1));
turbulenceResistanceFactor = Foam::max(turbulenceResistanceFactor, scalar(0));
turbulenceResistanceCap = Foam::max(turbulenceResistanceCap, scalar(0));
turbulenceStatePenaltyFactor = Foam::max(turbulenceStatePenaltyFactor, scalar(0));
turbulenceIntensity = Foam::max(turbulenceIntensity, scalar(1e-6));
turbulenceLengthScale = Foam::max(turbulenceLengthScale, scalar(1e-8));
kMinBound = Foam::max(kMinBound, scalar(1e-14));
epsilonMinBound = Foam::max(epsilonMinBound, scalar(1e-14));

if (scenarioVelocities.size() == 0)
{
  scenarioVelocities.setSize(1);
  scenarioVelocities[0] = defaultInletSpeed;
}
forAll(scenarioVelocities, sI)
{
  scenarioVelocities[sI] = Foam::max(scenarioVelocities[sI], scalar(SMALL));
}
if (scenarioWeights.size() != scenarioVelocities.size())
{
  if (scenarioWeights.size() == 1)
  {
    scenarioWeights.setSize(scenarioVelocities.size());
    forAll(scenarioWeights, sI)
    {
      scenarioWeights[sI] = 1.0;
    }
  }
  else
  {
    FatalErrorInFunction
      << "scenarioWeights size (" << scenarioWeights.size()
      << ") must match scenarioVelocities size (" << scenarioVelocities.size()
      << ") or be omitted."
      << exit(FatalError);
  }
}
scalar scenarioWeightSum(0.0);
forAll(scenarioWeights, sI)
{
  scenarioWeights[sI] = Foam::max(scenarioWeights[sI], scalar(0.0));
  scenarioWeightSum += scenarioWeights[sI];
}
if (scenarioWeightSum <= SMALL)
{
  scenarioWeightSum = scenarioWeights.size();
  forAll(scenarioWeights, sI)
  {
    scenarioWeights[sI] = 1.0;
  }
}
forAll(scenarioWeights, sI)
{
  scenarioWeights[sI] /= scenarioWeightSum;
}
const label nVelocityScenarios(scenarioVelocities.size());
scalarList scenarioThermalReference(nVelocityScenarios, -1.0);
scalarList scenarioPowerReference(nVelocityScenarios, -1.0);
scalarList scenarioMeanTLast(nVelocityScenarios, 0.0);
scalarList scenarioPowerLast(nVelocityScenarios, 0.0);
scalarList scenarioMaxTLast(nVelocityScenarios, 0.0);
Info<< "Velocity scenarios [m/s]: " << scenarioVelocities
    << ", normalized weights: " << scenarioWeights << endl;

volScalarField alphaResistanceScale
(
    IOobject
    (
        "alphaResistanceScale",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("alphaResistanceScale", dimless, 1.0),
    zeroGradientFvPatchScalarField::typeName
);

//****************************************************************//adj-fluid-a
Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Info<< "Reading/calculating face flux field phia\n" << endl;
surfaceScalarField phia
(
    IOobject
    (
        "phia",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(Ua) & mesh.Sf()
);
label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.setFluxRequired(pa.name());
//****************************************************************//adj-fluid-b
Info<< "Reading field pb\n" << endl;
volScalarField pb
(
    IOobject
    (
        "pb",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Info<< "Reading field Ub\n" << endl;
volVectorField Ub
(
    IOobject
    (
        "Ub",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);
Info<< "Reading/calculating face flux field phib\n" << endl;
surfaceScalarField phib
(
    IOobject
    (
        "phib",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::NO_WRITE
    ),
    linearInterpolate(Ub) & mesh.Sf()
);
label pbRefCell = 0;
scalar pbRefValue = 0.0;
setRefCell
(
    pb,
    simple.dict(),
    pbRefCell,
    pbRefValue
);
mesh.setFluxRequired(pb.name());

//*****************************************************//
volScalarField gsensPowerDiss
(
    IOobject
    (
        "gsensPowerDiss",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    -alphaMax*(1+qu)*qu/((qu+xh)*(qu+xh))*(U&Ua),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshPowerDiss(gsensPowerDiss);

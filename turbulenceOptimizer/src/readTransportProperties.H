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
label freezeTurbulenceUntilIteration(0);
Switch turbulenceProMode(false);
Switch useTurbulentResistance(true);
Switch frozenTurbulenceAdjoint(true);
Switch normalizeScenarioObjectives(true);
Switch enableThermalFallback(true);
Switch enforceTurbulenceBounds(true);
Switch updateTurbulenceInletFromVelocity(true);
Switch enableThermalDivergenceGuard(true);
Switch clipTemperatureField(true);
Switch enableAdjointHeatSourceClipping(true);
scalar turbulenceResistanceFactor(1.0);
scalar turbulenceResistanceCap(100.0);
scalar turbulenceStatePenaltyFactor(1.0);
scalar turbulenceIntensity(0.05);
scalar turbulenceLengthScale(1e-3);
scalar kMinBound(1e-10);
scalar epsilonMinBound(1e-10);
scalar kMaxBound(1e4);
scalar epsilonMaxBound(1e8);
scalar nutMaxBound(1e3);
scalar thermalFallbackResidualThreshold(1e-7);
label thermalFallbackIterationThreshold(15000);
scalar thermalFallbackBlend(0.3);
label thermalFallbackMaxAttempts(2);
scalar thermalFallbackDiffusivityFactor(3.0);
scalar thermalUnderRelaxation(0.7);
scalar thermalFieldAbsMax(1e6);
scalar temperatureLowerBound(-1000.0);
scalar temperatureUpperBound(1000.0);
scalar adjointHeatSourceClip(1e4);
scalar turbulenceModelCmu(0.09);

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
  turbulenceProMode =
    turbulenceOptimization.lookupOrDefault<Switch>("turbulenceProMode", turbulenceProMode);
  freezeTurbulenceUntilIteration =
    turbulenceOptimization.lookupOrDefault<label>
    (
      "freezeTurbulenceUntilIteration",
      freezeTurbulenceUntilIteration
    );
  useTurbulentResistance =
    turbulenceOptimization.lookupOrDefault<Switch>("useTurbulentResistance", useTurbulentResistance);
  frozenTurbulenceAdjoint =
    turbulenceOptimization.lookupOrDefault<Switch>("frozenTurbulenceAdjoint", frozenTurbulenceAdjoint);
  normalizeScenarioObjectives =
    turbulenceOptimization.lookupOrDefault<Switch>("normalizeScenarioObjectives", normalizeScenarioObjectives);
  enableThermalFallback =
    turbulenceOptimization.lookupOrDefault<Switch>("enableThermalFallback", enableThermalFallback);
  enforceTurbulenceBounds =
    turbulenceOptimization.lookupOrDefault<Switch>("enforceTurbulenceBounds", enforceTurbulenceBounds);
  updateTurbulenceInletFromVelocity =
    turbulenceOptimization.lookupOrDefault<Switch>
    (
      "updateTurbulenceInletFromVelocity",
      updateTurbulenceInletFromVelocity
    );
  enableThermalDivergenceGuard =
    turbulenceOptimization.lookupOrDefault<Switch>
    (
      "enableThermalDivergenceGuard",
      enableThermalDivergenceGuard
    );
  clipTemperatureField =
    turbulenceOptimization.lookupOrDefault<Switch>("clipTemperatureField", clipTemperatureField);
  enableAdjointHeatSourceClipping =
    turbulenceOptimization.lookupOrDefault<Switch>
    (
      "enableAdjointHeatSourceClipping",
      enableAdjointHeatSourceClipping
    );
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
  kMaxBound = turbulenceOptimization.lookupOrDefault<scalar>("kMaxBound", kMaxBound);
  epsilonMaxBound = turbulenceOptimization.lookupOrDefault<scalar>("epsilonMaxBound", epsilonMaxBound);
  nutMaxBound = turbulenceOptimization.lookupOrDefault<scalar>("nutMaxBound", nutMaxBound);
  thermalFallbackResidualThreshold =
    turbulenceOptimization.lookupOrDefault<scalar>
    (
      "thermalFallbackResidualThreshold",
      thermalFallbackResidualThreshold
    );
  thermalFallbackIterationThreshold =
    turbulenceOptimization.lookupOrDefault<label>
    (
      "thermalFallbackIterationThreshold",
      thermalFallbackIterationThreshold
    );
  thermalFallbackBlend =
    turbulenceOptimization.lookupOrDefault<scalar>("thermalFallbackBlend", thermalFallbackBlend);
  thermalFallbackMaxAttempts =
    turbulenceOptimization.lookupOrDefault<label>("thermalFallbackMaxAttempts", thermalFallbackMaxAttempts);
  thermalFallbackDiffusivityFactor =
    turbulenceOptimization.lookupOrDefault<scalar>
    (
      "thermalFallbackDiffusivityFactor",
      thermalFallbackDiffusivityFactor
    );
  thermalUnderRelaxation =
    turbulenceOptimization.lookupOrDefault<scalar>("thermalUnderRelaxation", thermalUnderRelaxation);
  thermalFieldAbsMax =
    turbulenceOptimization.lookupOrDefault<scalar>("thermalFieldAbsMax", thermalFieldAbsMax);
  temperatureLowerBound =
    turbulenceOptimization.lookupOrDefault<scalar>("temperatureLowerBound", temperatureLowerBound);
  temperatureUpperBound =
    turbulenceOptimization.lookupOrDefault<scalar>("temperatureUpperBound", temperatureUpperBound);
  adjointHeatSourceClip =
    turbulenceOptimization.lookupOrDefault<scalar>("adjointHeatSourceClip", adjointHeatSourceClip);
  turbulenceModelCmu =
    turbulenceOptimization.lookupOrDefault<scalar>("turbulenceModelCmu", turbulenceModelCmu);
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
freezeTurbulenceUntilIteration = Foam::max(freezeTurbulenceUntilIteration, label(0));
turbulenceResistanceFactor = Foam::max(turbulenceResistanceFactor, scalar(0));
turbulenceResistanceCap = Foam::max(turbulenceResistanceCap, scalar(0));
turbulenceStatePenaltyFactor = Foam::max(turbulenceStatePenaltyFactor, scalar(0));
turbulenceIntensity = Foam::max(turbulenceIntensity, scalar(1e-6));
turbulenceLengthScale = Foam::max(turbulenceLengthScale, scalar(1e-8));
kMinBound = Foam::max(kMinBound, scalar(1e-14));
epsilonMinBound = Foam::max(epsilonMinBound, scalar(1e-14));
kMaxBound = Foam::max(kMaxBound, kMinBound + scalar(SMALL));
epsilonMaxBound = Foam::max(epsilonMaxBound, epsilonMinBound + scalar(SMALL));
nutMaxBound = Foam::max(nutMaxBound, scalar(SMALL));
thermalFallbackResidualThreshold = Foam::max(thermalFallbackResidualThreshold, scalar(1e-14));
thermalFallbackIterationThreshold = Foam::max(thermalFallbackIterationThreshold, label(1));
thermalFallbackBlend = Foam::min(Foam::max(thermalFallbackBlend, scalar(0.0)), scalar(1.0));
thermalFallbackMaxAttempts = Foam::max(thermalFallbackMaxAttempts, label(1));
thermalFallbackDiffusivityFactor = Foam::max(thermalFallbackDiffusivityFactor, scalar(1.0));
thermalUnderRelaxation = Foam::min(Foam::max(thermalUnderRelaxation, scalar(0.05)), scalar(1.0));
thermalFieldAbsMax = Foam::max(thermalFieldAbsMax, scalar(1.0));
if (temperatureUpperBound <= temperatureLowerBound)
{
  const scalar tCenter = 0.5*(temperatureUpperBound + temperatureLowerBound);
  temperatureLowerBound = tCenter - 1.0;
  temperatureUpperBound = tCenter + 1.0;
}
adjointHeatSourceClip = Foam::max(adjointHeatSourceClip, scalar(SMALL));
turbulenceModelCmu = Foam::max(turbulenceModelCmu, scalar(1e-6));

if (!turbulenceProMode)
{
  // Non-publication safety heuristics are centrally gated by turbulenceProMode.
  enableThermalFallback = false;
  enforceTurbulenceBounds = false;
  enableThermalDivergenceGuard = false;
  clipTemperatureField = false;
  enableAdjointHeatSourceClipping = false;
  thermalUnderRelaxation = 1.0;
}

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
Info<< "turbulenceProMode: " << (turbulenceProMode ? "on" : "off") << endl;
Info<< "Heuristic stabilizers: "
    << (turbulenceProMode ? "enabled" : "disabled")
    << " (thermal fallback / clipping / guarding)" << endl;
Info<< "Stability controls: freezeTurbulenceUntilIteration="
    << freezeTurbulenceUntilIteration
    << ", turbulence bounds=[" << kMinBound << "," << kMaxBound
    << "] epsilon=[" << epsilonMinBound << "," << epsilonMaxBound
    << "], thermal fallback attempts=" << thermalFallbackMaxAttempts
    << ", thermalUnderRelaxation=" << thermalUnderRelaxation
    << ", thermalFieldAbsMax=" << thermalFieldAbsMax
    << endl;

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

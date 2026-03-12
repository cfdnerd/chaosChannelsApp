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
Switch useTurbulentResistance(true);
Switch frozenTurbulenceAdjoint(true);
Switch enableThermalFallback(true);
Switch enforceTurbulenceBounds(true);
Switch updateTurbulenceInletFromVelocity(true);
Switch enableTbBcSmoothingMode(false);
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
    << "Missing required patch 'inlet' for single-case setup."
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
scalar inletSpeed(Foam::max(defaultInletSpeed, scalar(SMALL)));
if (tuneOptParameters.found("turbulenceOptimization"))
{
  const dictionary& turbulenceOptimization = tuneOptParameters.subDict("turbulenceOptimization");
  // Keep only the pro-mode switch configurable from tuneOptParameters.
  turbulenceProMode =
    turbulenceOptimization.lookupOrDefault<Switch>("turbulenceProMode", turbulenceProMode);
}
enableTbBcSmoothingMode =
  tuneOptParameters.lookupOrDefault<Switch>
  (
    "enableTbBcSmoothingMode",
    enableTbBcSmoothingMode
  );

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

volScalarField tbBcSmoothingMode
(
    IOobject
    (
        "tbBcSmoothingMode",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("tbBcSmoothingMode", dimless, enableTbBcSmoothingMode ? 1.0 : 0.0),
    zeroGradientFvPatchScalarField::typeName
);
tbBcSmoothingMode.correctBoundaryConditions();
Info<< "Inlet speed [m/s]: " << inletSpeed
    << " with direction " << inletVelocityDirection << endl;
Info<< "turbulenceProMode: " << (turbulenceProMode ? "on" : "off") << endl;
Info<< "Heuristic stabilizers: always enabled"
    << " (turbulence bounds / clipping / divergence guarding)" << endl;
Info<< "tbBcSmoothingMode: "
    << (enableTbBcSmoothingMode ? "enabled" : "disabled")
    << " (independent toggle)" << endl;
Info<< "Thermal fallback path: "
    << (enableThermalFallback ? "enabled" : "disabled")
    << " (general, independent of turbulenceProMode)" << endl;
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

Info<< "Reading optimization parameters\n" << endl;
IOdictionary optProperties
(
   IOobject
   (
       "optProperties",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE
   )
);
IOdictionary tuneOptParameters
(
   IOobject
   (
       "tuneOptParameters",
       runTime.constant(),
       mesh,
       IOobject::MUST_READ_IF_MODIFIED,
       IOobject::NO_WRITE
   )
);
//optimization parameters
scalar raa0(readScalar(optProperties.lookup("raa0")));
scalar mma_init(readScalar(optProperties.lookup("mma_init")));
scalar mma_dec(readScalar(optProperties.lookup("mma_dec")));
scalar mma_inc(readScalar(optProperties.lookup("mma_inc")));
scalar movlim(readScalar(optProperties.lookup("movlim")));
scalar voluse(readScalar(optProperties.lookup("voluse")));
scalar PowerDiss0(readScalar(optProperties.lookup("PowerDiss0")));
scalar PowerDissMax(readScalar(optProperties.lookup("PowerDissMax")));
scalar PowerDissRelax(readScalar(optProperties.lookup("PowerDissRelax")));
scalar filterR(readScalar(optProperties.lookup("filterR")));
dimensionedScalar b("b", dimensionSet(0,-2,0,0,0,0,0),1.0);//filter radius of the PDE filter
dimensionedScalar b1("b1", dimensionSet(0,2,-2,-2,0,0,0),1.0);
Switch solid_area(optProperties.lookup("solid_area"));
Switch fluid_area(optProperties.lookup("fluid_area"));
Switch test_area(optProperties.lookup("test_area"));
scalar Pnorm(readScalar(optProperties.lookup("Pnorm")));
scalar objFunction(readScalar(optProperties.lookup("objFunction")));
Switch turbulenceProMode(false);
if (tuneOptParameters.found("turbulenceOptimization"))
{
  const dictionary& turbulenceOptimization = tuneOptParameters.subDict("turbulenceOptimization");
  turbulenceProMode =
    turbulenceOptimization.lookupOrDefault<Switch>("turbulenceProMode", turbulenceProMode);
}
const dictionary& convergenceControl = tuneOptParameters.subDict("convergenceControl");
label optMinIterations(convergenceControl.lookupOrDefault<label>("optMinIterations", 50));
label optStallWindow(convergenceControl.lookupOrDefault<label>("optStallWindow", 5));
scalar optObjectiveTol(convergenceControl.lookupOrDefault<scalar>("optObjectiveTol", 1e-3));
scalar optConstraintTol(convergenceControl.lookupOrDefault<scalar>("optConstraintTol", 1e-3));
scalar optStepTol(convergenceControl.lookupOrDefault<scalar>("optStepTol", 5e-3));
scalar optPowerConstraintWeight
(
  convergenceControl.lookupOrDefault<scalar>("optPowerConstraintWeight", 1.0)
);
scalar optVolumeConstraintWeight
(
  convergenceControl.lookupOrDefault<scalar>("optVolumeConstraintWeight", 1.0)
);
scalar rhoMin(1e-3);
scalar projectionBetaRate(0.2);
scalar projectionBetaMax(100.0);
scalar projectionBetaInitial(1.0);
scalar projectionEtaLower(1e-4);
scalar projectionEtaUpper(1.0);
scalar projectionEtaTolerance(1e-4);
label projectionUpdateInterval(5);
scalar projectionObjectiveTolForBeta(5e-3);
scalar projectionStepTolForBeta(5e-3);
scalar projectionConstraintTolForBeta(5e-2);
Switch projectionRequireConstraintForBeta(true);
scalar alphaRampDivisor(7.0);
label alphaRampIterations(63);
label alphaGrowthStartIteration(100);
scalar alphaGrowthFactor(1.05);
scalar quInitial(0.005);
label quRampStartIteration(80);
scalar quGrowthRate(1e-4);
scalar quMax(0.01);
scalar filterRStart(filterR);
scalar filterREnd(filterR);
label filterRContinuationStart(0);
label filterRContinuationEnd(1);
// Continuation/projection tuning is intentionally hardcoded.
rhoMin = Foam::min(Foam::max(rhoMin, scalar(1e-6)), scalar(0.25));
projectionBetaInitial = Foam::max(projectionBetaInitial, scalar(1e-6));
projectionEtaLower = Foam::max(projectionEtaLower, scalar(1e-8));
projectionEtaUpper = Foam::min(Foam::max(projectionEtaUpper, projectionEtaLower + 1e-6), scalar(1.0));
projectionEtaTolerance = Foam::max(projectionEtaTolerance, scalar(1e-8));
projectionUpdateInterval = Foam::max(projectionUpdateInterval, label(1));
projectionObjectiveTolForBeta = Foam::max(projectionObjectiveTolForBeta, scalar(1e-8));
projectionStepTolForBeta = Foam::max(projectionStepTolForBeta, scalar(1e-8));
projectionConstraintTolForBeta = Foam::max(projectionConstraintTolForBeta, scalar(0.0));
alphaRampDivisor = Foam::max(alphaRampDivisor, scalar(1.0));
alphaGrowthFactor = Foam::max(alphaGrowthFactor, scalar(1.0));
quInitial = Foam::max(quInitial, scalar(1e-6));
quGrowthRate = Foam::max(quGrowthRate, scalar(0.0));
quMax = Foam::max(quMax, quInitial);
filterRStart = Foam::max(filterRStart, scalar(1e-6));
filterREnd = Foam::max(filterREnd, scalar(1e-6));
if (optPowerConstraintWeight <= 0.0)
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optPowerConstraintWeight must be > 0. Found "
    << optPowerConstraintWeight
    << exit(FatalError);
}
if (optVolumeConstraintWeight <= 0.0)
{
  FatalErrorInFunction
    << "tuneOptParameters.convergenceControl.optVolumeConstraintWeight must be > 0. Found "
    << optVolumeConstraintWeight
    << exit(FatalError);
}
Info<< "Reading pseudo-density field\n" << endl;
volScalarField xp
(
    IOobject
    (
        "xp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField xh
(
    IOobject
    (
        "xh",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    scalar(voluse),
    zeroGradientFvPatchScalarField::typeName
);
volScalarField x(xp);
volScalarField drho(x);

// find the cells in the non-design area
labelList cells_solid,cells_fluid,cells_test;
if(solid_area)
{
  word zoneName="solid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_solid=mesh.cellZones()[zoneID];
  forAll(cells_solid, i)
  {
     x[cells_solid[i]]=0;
  }
}
if(fluid_area)
{
  word zoneName="fluid_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_fluid=mesh.cellZones()[zoneID];
  forAll(cells_fluid, i)
  {
     x[cells_fluid[i]]=1.0;
  }
}
if(test_area)
{
  word zoneName="test_area";
  label zoneID=mesh.cellZones().findZoneID(zoneName);
  cells_test=mesh.cellZones()[zoneID];
  forAll(cells_test, i)
  {
     x[cells_test[i]]=1.0;
     xh[cells_test[i]]=1.0;
  }
}
volScalarField gsensVol
(
    IOobject
    (
        "gsensVol",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    xh,
    zeroGradientFvPatchScalarField::typeName
);
volScalarField gsenshVol(xh);

// find the cells on the boundary
dictionary conFunc = mesh.solutionDict().subDict("constraintFunctionDict_flow");
int nObjPatch =conFunc.lookupOrDefault<scalar>("numberConstraintPatches",2);
wordList conPatchNames = conFunc.lookup("constraintPatchesNames");
label conPatchList[nObjPatch];//patch list
int iLoop;
for(iLoop=0; iLoop<nObjPatch; iLoop++)
{
    conPatchList[iLoop]=mesh.boundaryMesh().findPatchID(conPatchNames[iLoop]);
}   

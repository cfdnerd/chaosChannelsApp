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
Switch adaptiveHeavisideProjection(false);
adaptiveHeavisideProjection =
  tuneOptParameters.lookupOrDefault<Switch>
  (
    "adaptiveHeavisideProjection",
    adaptiveHeavisideProjection
  );
Switch normalizeScenarioObjectives(false);
normalizeScenarioObjectives =
  tuneOptParameters.lookupOrDefault<Switch>
  (
    "normalizeScenarioObjectives",
    normalizeScenarioObjectives
  );
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
scalar projectionBetaGrowthFactor(1.25);
scalar projectionEtaLower(1e-4);
scalar projectionEtaUpper(1.0);
scalar projectionEtaTolerance(1e-4);
scalar projectionEtaEps(1e-6);
label projectionUpdateInterval(25);
label projectionStartIteration(75);
scalar projectionObjectiveTolForBeta(5e-3);
scalar projectionStepTolForBeta(5e-3);
scalar projectionConstraintViolationTolForBeta(10.0*optConstraintTol);
label powerConstraintTighteningStartIteration(0);
scalar powerConstraintTighteningRate(0.0);
scalar alphaRampDivisor(40.0);
label alphaRampIterations(80);
label alphaGrowthStartIteration(1000000);
scalar alphaGrowthFactor(1.0);
scalar alphaContinuationCapMultiplier(4.0);
// Staged q-continuation to suppress gray designs without runaway alpha growth.
scalar quInitial(0.01);
label quStage1Iteration(50);
label quStage2Iteration(100);
label quStage3Iteration(150);
scalar quStage1Value(0.01);
scalar quStage2Value(0.1);
scalar quStage3Value(1.0);
scalar quStage4Value(10.0);
scalar quMax(10.0);
scalar continuationConstraintViolationFreeze(10.0*optConstraintTol);
scalar continuationBlendFast(0.25);
scalar continuationBlendSlow(0.05);
scalar mmaConstraintDampingGain(1.5);
scalar mmaConstraintDampingMin(0.3);
scalar mmaConstraintDampingActivation(10.0*optConstraintTol);
scalar filterRStart(filterR);
scalar filterREnd(filterR);
label filterRContinuationStart(0);
label filterRContinuationEnd(1);
scalar continuationConstraintViolationSoftTol(2.0*continuationConstraintViolationFreeze);
scalar continuationProgressMin(0.2);
scalar grayStallContinuationBlendMin(0.1);
scalar projectionConstraintViolationSoftTol(2.0*projectionConstraintViolationTolForBeta);
scalar projectionBetaGrowthFactorMin(1.02);
scalar projectionUnstableGrowthScale(0.25);
scalar grayStallThreshold(0.85);
scalar grayStallObjectiveTol(5e-5);
scalar grayStallStepTol(2e-3);
label grayStallWindow(8);
scalar grayStallBetaGrowthFactor(1.10);
if (tuneOptParameters.found("stallRecoveryControl"))
{
  const dictionary& stallRecoveryControl = tuneOptParameters.subDict("stallRecoveryControl");
  continuationConstraintViolationSoftTol =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "continuationConstraintViolationSoftTol",
      continuationConstraintViolationSoftTol
    );
  continuationProgressMin =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "continuationProgressMin",
      continuationProgressMin
    );
  grayStallContinuationBlendMin =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "grayStallContinuationBlendMin",
      grayStallContinuationBlendMin
    );
  projectionConstraintViolationSoftTol =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "projectionConstraintViolationSoftTol",
      projectionConstraintViolationSoftTol
    );
  projectionBetaGrowthFactorMin =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "projectionBetaGrowthFactorMin",
      projectionBetaGrowthFactorMin
    );
  projectionUnstableGrowthScale =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "projectionUnstableGrowthScale",
      projectionUnstableGrowthScale
    );
  grayStallThreshold =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "grayStallThreshold",
      grayStallThreshold
    );
  grayStallObjectiveTol =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "grayStallObjectiveTol",
      grayStallObjectiveTol
    );
  grayStallStepTol =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "grayStallStepTol",
      grayStallStepTol
    );
  grayStallWindow =
    stallRecoveryControl.lookupOrDefault<label>
    (
      "grayStallWindow",
      grayStallWindow
    );
  grayStallBetaGrowthFactor =
    stallRecoveryControl.lookupOrDefault<scalar>
    (
      "grayStallBetaGrowthFactor",
      grayStallBetaGrowthFactor
    );
}
// Continuation/projection tuning is intentionally hardcoded.
rhoMin = Foam::min(Foam::max(rhoMin, scalar(1e-6)), scalar(0.25));
projectionBetaInitial = Foam::max(projectionBetaInitial, scalar(1e-6));
projectionBetaGrowthFactor = Foam::max(projectionBetaGrowthFactor, scalar(1.0));
projectionEtaEps = Foam::max(projectionEtaEps, scalar(1e-8));
projectionEtaLower = Foam::max(projectionEtaLower, projectionEtaEps);
projectionEtaUpper = Foam::min(projectionEtaUpper, scalar(1.0)-projectionEtaEps);
if (projectionEtaUpper <= projectionEtaLower)
{
  projectionEtaUpper = Foam::min(projectionEtaLower + 10.0*projectionEtaEps, scalar(1.0)-projectionEtaEps);
  if (projectionEtaUpper <= projectionEtaLower)
  {
    projectionEtaLower = projectionEtaEps;
    projectionEtaUpper = scalar(1.0)-projectionEtaEps;
  }
}
projectionEtaTolerance = Foam::max(projectionEtaTolerance, scalar(1e-8));
projectionUpdateInterval = Foam::max(projectionUpdateInterval, label(1));
projectionStartIteration = Foam::max(projectionStartIteration, label(1));
projectionObjectiveTolForBeta = Foam::max(projectionObjectiveTolForBeta, scalar(1e-8));
projectionStepTolForBeta = Foam::max(projectionStepTolForBeta, scalar(1e-8));
projectionConstraintViolationTolForBeta =
  Foam::max(projectionConstraintViolationTolForBeta, scalar(0.0));
projectionStepTolForBeta =
  Foam::max
  (
    projectionStepTolForBeta,
    scalar(0.02)*Foam::max(Foam::mag(movlim), scalar(1e-6))
  );
powerConstraintTighteningStartIteration =
  Foam::max(powerConstraintTighteningStartIteration, label(0));
powerConstraintTighteningRate = Foam::max(powerConstraintTighteningRate, scalar(0.0));
alphaRampDivisor = Foam::max(alphaRampDivisor, scalar(1.0));
alphaGrowthFactor = Foam::max(alphaGrowthFactor, scalar(1.0));
alphaContinuationCapMultiplier = Foam::max(alphaContinuationCapMultiplier, scalar(1.0));
quInitial = Foam::max(quInitial, scalar(1e-6));
quStage1Iteration = Foam::max(quStage1Iteration, label(1));
quStage2Iteration = Foam::max(quStage2Iteration, quStage1Iteration + 1);
quStage3Iteration = Foam::max(quStage3Iteration, quStage2Iteration + 1);
quStage1Value = Foam::max(quStage1Value, scalar(1e-6));
quStage2Value = Foam::max(quStage2Value, quStage1Value);
quStage3Value = Foam::max(quStage3Value, quStage2Value);
quStage4Value = Foam::max(quStage4Value, quStage3Value);
quMax = Foam::max(quMax, quStage4Value);
continuationConstraintViolationFreeze =
  Foam::max(continuationConstraintViolationFreeze, scalar(0.0));
continuationConstraintViolationSoftTol =
  Foam::max
  (
    continuationConstraintViolationSoftTol,
    continuationConstraintViolationFreeze + scalar(1e-8)
  );
continuationBlendFast =
  Foam::min(Foam::max(continuationBlendFast, scalar(1e-6)), scalar(1.0));
continuationBlendSlow =
  Foam::min(Foam::max(continuationBlendSlow, scalar(1e-6)), scalar(1.0));
continuationProgressMin =
  Foam::min(Foam::max(continuationProgressMin, scalar(1e-3)), scalar(1.0));
grayStallContinuationBlendMin =
  Foam::min(Foam::max(grayStallContinuationBlendMin, scalar(1e-6)), scalar(1.0));
projectionConstraintViolationSoftTol =
  Foam::max
  (
    projectionConstraintViolationSoftTol,
    projectionConstraintViolationTolForBeta + scalar(1e-8)
  );
projectionBetaGrowthFactorMin =
  Foam::min
  (
    Foam::max(projectionBetaGrowthFactorMin, scalar(1.0)),
    projectionBetaGrowthFactor
  );
projectionUnstableGrowthScale =
  Foam::min(Foam::max(projectionUnstableGrowthScale, scalar(0.0)), scalar(1.0));
grayStallThreshold =
  Foam::min(Foam::max(grayStallThreshold, scalar(0.0)), scalar(1.0));
grayStallObjectiveTol = Foam::max(grayStallObjectiveTol, scalar(1e-10));
grayStallStepTol = Foam::max(grayStallStepTol, scalar(1e-10));
grayStallWindow = Foam::max(grayStallWindow, label(1));
grayStallBetaGrowthFactor =
  Foam::max(grayStallBetaGrowthFactor, projectionBetaGrowthFactorMin);
mmaConstraintDampingGain = Foam::max(mmaConstraintDampingGain, scalar(0.0));
mmaConstraintDampingMin =
  Foam::min(Foam::max(mmaConstraintDampingMin, scalar(1e-3)), scalar(1.0));
mmaConstraintDampingActivation =
  Foam::max(mmaConstraintDampingActivation, scalar(0.0));
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
Info<< "normalizeScenarioObjectives: "
    << (normalizeScenarioObjectives ? "enabled" : "disabled")
    << " (first valid-iteration fixed MeanT/PDCon scaling)" << endl;
Info<< "Power constraint scaling: PDCon always uses PowerDiss/PowerDiss0; "
    << "normalizeScenarioObjectives only adds/removes the extra first-valid-iteration scenario scaling."
    << endl;
Info<< "Continuation settings: alphaRampDivisor=" << alphaRampDivisor
    << ", alphaRampIterations=" << alphaRampIterations
    << ", alphaCapMultiplier=" << alphaContinuationCapMultiplier
    << ", betaStartIteration=" << projectionStartIteration
    << ", betaUpdateInterval=" << projectionUpdateInterval
    << ", betaGrowthFactor=" << projectionBetaGrowthFactor
    << ", betaMax=" << projectionBetaMax
    << ", qStages=[" << quStage1Value << "," << quStage2Value
    << "," << quStage3Value << "," << quStage4Value << "]"
    << endl;
Info<< "Stability profile: projectionConstraintViolationTolForBeta="
    << projectionConstraintViolationTolForBeta
    << ", continuationConstraintViolationFreeze="
    << continuationConstraintViolationFreeze
    << ", mmaConstraintDampingGain=" << mmaConstraintDampingGain
    << ", mmaConstraintDampingMin=" << mmaConstraintDampingMin
    << ", mmaConstraintDampingActivation=" << mmaConstraintDampingActivation
    << endl;
Info<< "Stall recovery: continuationSoftTol=" << continuationConstraintViolationSoftTol
    << ", continuationProgressMin=" << continuationProgressMin
    << ", continuationBlendMin=" << grayStallContinuationBlendMin
    << ", projectionSoftTol=" << projectionConstraintViolationSoftTol
    << ", betaGrowthMin=" << projectionBetaGrowthFactorMin
    << ", unstableGrowthScale=" << projectionUnstableGrowthScale
    << ", grayThreshold=" << grayStallThreshold
    << ", grayWindow=" << grayStallWindow
    << ", grayBetaGrowth=" << grayStallBetaGrowthFactor
    << endl;
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

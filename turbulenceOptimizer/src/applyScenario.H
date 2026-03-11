const scalar scenarioVelocity = scenarioVelocities[scenarioI];
const vector scenarioInletVelocity = inletVelocityDirection*scenarioVelocity;
U.boundaryFieldRef()[inletPatchID] = scenarioInletVelocity;
Ua.boundaryFieldRef()[inletPatchID] = scenarioInletVelocity;

if (updateTurbulenceInletFromVelocity)
{
  const scalar inletSpeed = Foam::max(mag(scenarioInletVelocity), scalar(SMALL));
  const scalar turbulenceK =
    Foam::max(1.5*Foam::sqr(turbulenceIntensity*inletSpeed), kMinBound);
  const scalar turbulenceEpsilon =
    Foam::max
    (
      Foam::pow(turbulenceModelCmu, 0.75)
     *Foam::pow(turbulenceK, 1.5)
     /Foam::max(turbulenceLengthScale, scalar(SMALL)),
      epsilonMinBound
    );
  if
  (
    mesh.foundObject<volScalarField>("k")
    &&
    mesh.foundObject<volScalarField>("epsilon")
  )
  {
    volScalarField& kField =
      const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("k"));
    volScalarField& epsilonField =
      const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("epsilon"));
    kField.boundaryFieldRef()[inletPatchID] = turbulenceK;
    epsilonField.boundaryFieldRef()[inletPatchID] = turbulenceEpsilon;
    kField.correctBoundaryConditions();
    epsilonField.correctBoundaryConditions();
  }
}

U.correctBoundaryConditions();
Ua.correctBoundaryConditions();

Info<< "Running scenario " << (scenarioI + 1) << "/" << nVelocityScenarios
    << " with inlet speed " << scenarioVelocity << " m/s and weight "
    << scenarioWeights[scenarioI] << endl;

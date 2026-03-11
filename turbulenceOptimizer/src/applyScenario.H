const scalar scenarioVelocity = scenarioVelocities[scenarioI];
const vector scenarioInletVelocity = inletVelocityDirection*scenarioVelocity;
U.boundaryField()[inletPatchID] == scenarioInletVelocity;
Ua.boundaryField()[inletPatchID] == scenarioInletVelocity;
U.correctBoundaryConditions();
Ua.correctBoundaryConditions();

if (mesh.foundObject<volScalarField>("k") && mesh.foundObject<volScalarField>("epsilon"))
{
  volScalarField& kField =
    const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("k"));
  volScalarField& epsilonField =
    const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("epsilon"));

  const scalar kInlet =
    Foam::max(1.5*Foam::sqr(scenarioVelocity*turbulenceIntensity), kMinBound);
  const scalar epsilonInlet =
    Foam::max
    (
      Foam::pow(0.09, 0.75)*Foam::pow(kInlet, 1.5)/turbulenceLengthScale,
      epsilonMinBound
    );
  kField.boundaryField()[inletPatchID] == kInlet;
  epsilonField.boundaryField()[inletPatchID] == epsilonInlet;
  kField.max(dimensionedScalar("kMinScenario", kField.dimensions(), kMinBound));
  epsilonField.max
  (
    dimensionedScalar("epsilonMinScenario", epsilonField.dimensions(), epsilonMinBound)
  );
  kField.correctBoundaryConditions();
  epsilonField.correctBoundaryConditions();
}

Info<< "Running scenario " << (scenarioI + 1) << "/" << nVelocityScenarios
    << " with inlet speed " << scenarioVelocity << " m/s and weight "
    << scenarioWeights[scenarioI] << endl;

const scalar scenarioVelocity = scenarioVelocities[scenarioI];
const vector scenarioInletVelocity = inletVelocityDirection*scenarioVelocity;
U.boundaryFieldRef()[inletPatchID] = scenarioInletVelocity;
Ua.boundaryFieldRef()[inletPatchID] = scenarioInletVelocity;
U.correctBoundaryConditions();
Ua.correctBoundaryConditions();

Info<< "Running scenario " << (scenarioI + 1) << "/" << nVelocityScenarios
    << " with inlet speed " << scenarioVelocity << " m/s and weight "
    << scenarioWeights[scenarioI] << endl;

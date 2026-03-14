Info<< "\n run Reaction-Diffusion Level Set Method (RD-LSM) \n" << endl;

Info<< "LSM: Updating Lagrangian multipliers..." << endl;
const scalar lsmKpVol = 0.5;
lsmLambda = lsmLambda + lsmKpVol * V;

const scalar P_val = (PowerDiss / Foam::max(PowerDissRelax - opt * 0.2, PowerDissMax) - 1.0);
const scalar lsmKpPower = optPowerConstraintWeight;
if (P_val > 0.0) 
{
    lsmLambdaPower = lsmLambdaPower + lsmKpPower * P_val;
} 
else 
{
    lsmLambdaPower = Foam::max(0.0, lsmLambdaPower - 0.01);
}

Info<< "LSM: Calculating Sens_phi..." << endl;
// Strip dimensions for the Phase-Field Sensitivity calculation to prevent FOAM FATAL ERROR
// Dimensions of fsensMeanT are usually [0 2 -3 0 0 0 0], gsensVol is [0 0 0 0 0 0 0].
volScalarField Sens_phi
(
    IOobject("Sens_phi", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
    mesh,
    dimensionedScalar("zero", dimless, 0.0)
);

forAll(Sens_phi, i)
{
    Sens_phi[i] = -0.5 * fsensMeanT[i] 
                 - 0.5 * lsmLambda * gsensVol[i] 
                 - 0.5 * lsmLambdaPower * gsensPowerDiss[i];
}

if (enableMaxFeatureConstraint) 
{
    forAll(Sens_phi, i)
    {
        Sens_phi[i] -= 0.5 * maxFeatureConstraintWeight * VEroded * gsensMaxFeature[i];
    }
}


scalar maxSens = gMax(mag(Sens_phi.primitiveField()));
maxSens = Foam::max(maxSens, 1e-10);
Sens_phi = Sens_phi / maxSens;

Info<< "LSM: Advancing interface..." << endl;
scalar lsmStep = 0.1;
phiLSM = phiLSM - lsmStep * Sens_phi;

Info<< "LSM: Diffusion step..." << endl;
volScalarField phiLSMRaw = phiLSM;
const auto lsmDiffusionSolvePerf = solve(fvm::laplacian(phiLSM) - fvm::Sp(b, phiLSM) == -phiLSMRaw * b);
Info<< "LSM: Step complete." << endl;



// 6. Phase Field Bounding (Reaction well limits)
forAll(phiLSM, i)
{
    phiLSM[i] = Foam::max(Foam::min(phiLSM[i], 1.0), -1.0);
}

// 7. Map back to continuous density parameters for the standard solver algorithms
forAll(x, i)
{
    x[i] = 0.5 * (1.0 - phiLSM[i]); // phi=1 (Solid) -> x=0; phi=-1 (Fluid) -> x=1.
    x[i] = Foam::max(Foam::min(x[i], 1.0), 0.0);
}

// 8. Re-enforce zones after LSM update to ensure no drift in inlet/outlets
if(solid_area)
{
   forAll(cells_solid, j)
   {
      x[cells_solid[j]] = 0.0;
      phiLSM[cells_solid[j]] = 1.0;
   }
}
if(fluid_area)
{
   forAll(cells_fluid, j)
   {
      x[cells_fluid[j]] = 1.0;
      phiLSM[cells_fluid[j]] = -1.0;
   }
}

forAll(xmma, i)
{
    xmma[i] = x[i]; // Sync GUI/Paraview output arrays
}

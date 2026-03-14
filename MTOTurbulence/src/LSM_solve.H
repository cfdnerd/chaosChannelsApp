Info<< "\n run Reaction-Diffusion Level Set Method (RD-LSM) \n" << endl;

// 1. Update Augmented Lagrangian multiplier for Volume Constraint
// If V > 0, we have too much fluid. We need lsmLambda > 0 to drive Sens_phi negative (which drives phi positive/solid).
const scalar lsmKpVol = 0.5;
lsmLambda = lsmLambda + lsmKpVol * V;

// 2. Update Augmented Lagrangian multiplier for Power Dissipation Constraint
// If P_val > 0, power dissipation is too high. 
const scalar P_val = (PowerDiss / Foam::max(PowerDissRelax - opt * 0.2, PowerDissMax) - 1.0);
const scalar lsmKpPower = optPowerConstraintWeight;
if (P_val > 0.0) 
{
    lsmLambdaPower = lsmLambdaPower + lsmKpPower * P_val;
} 
else 
{
    // Slowly relax multiplier if constraint is over-satisfied
    lsmLambdaPower = Foam::max(0.0, lsmLambdaPower - 0.01);
}

// 3. Calculate Phase-Field Sensitivity (Gradient of Lagrangian w.r.t Phi)
// dL/dphi = -0.5 * dL/dxh (since xh = 0.5 * (1 - phi))
volScalarField Sens_phi = -0.5 * fsensMeanT 
                         - 0.5 * lsmLambda * gsensVol 
                         - 0.5 * lsmLambdaPower * gsensPowerDiss;

// Incorporate Maximum Feature Size penalty if enabled
if (enableMaxFeatureConstraint) 
{
    Sens_phi -= 0.5 * maxFeatureConstraintWeight * VEroded * gsensMaxFeature;
}

// Normalize the Raw Sensitivity
scalar maxSens = gMax(mag(Sens_phi.primitiveField()));
maxSens = Foam::max(maxSens, 1e-10);
Sens_phi = Sens_phi / maxSens;

// 4. Advance the Level-Set Interface (Reaction Step)
scalar lsmStep = 0.1; // Limit interface advection step strictly
phiLSM = phiLSM - lsmStep * Sens_phi;

// 5. Implicit Smoothing & Re-initialization (Diffusion Step)
// This strictly maintains the continuous interface width and enforces topological integrity.
// Solving: fvm::laplacian(phiLSM) - fvm::Sp(b, phiLSM) = -b * phiLSM_raw
volScalarField phiLSMRaw = phiLSM;
const auto lsmDiffusionSolvePerf = solve(fvm::laplacian(phiLSM) - fvm::Sp(b, phiLSM) + phiLSMRaw * b);

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
    xmma[i] = x[i]; // Sync GUI/Paraview output arrays
}

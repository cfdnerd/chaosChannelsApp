for(i=0;i<60;i++)
{
MRF.correctBoundaryVelocity(U);
tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  //- fvm::laplacian(nu,U)
  + fvm::Sp(alpha, U)
 ==
    fvOptions(U)
);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (simple.momentumPredictor())
{
    const auto ULinearSolvePerformance = solve(UEqn == -fvc::grad(p));
    recordLinearSolve
    (
        ULinearSolveStats,
        cmptMax(ULinearSolvePerformance.initialResidual()),
        cmptMax(ULinearSolvePerformance.finalResidual()),
        cmptMax(ULinearSolvePerformance.nIterations())
    );

    fvOptions.correct(U);
}
//****************************************//
volScalarField rAU(1.0/UEqn.A());
volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
surfaceScalarField phiHbyA("phiHbyA", fvc::flux(HbyA));
MRF.makeRelative(phiHbyA);
adjustPhi(phiHbyA, U, p);

tmp<volScalarField> rAtU(rAU);

if (simple.consistent())
{
    rAtU = 1.0/(1.0/rAU - UEqn.H1());
    phiHbyA +=
        fvc::interpolate(rAtU() - rAU)*fvc::snGrad(p)*mesh.magSf();
    HbyA -= (rAU - rAtU())*fvc::grad(p);
}

tUEqn.clear();

// Update the pressure BCs to ensure flux consistency
constrainPressure(p, U, phiHbyA, rAtU(), MRF);

// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn
    (
        fvm::laplacian(rAtU(), p) == fvc::div(phiHbyA)
    );

    pEqn.setReference(pRefCell, pRefValue);

    const auto pLinearSolvePerformance = pEqn.solve();
    recordLinearSolve
    (
        pLinearSolveStats,
        pLinearSolvePerformance.initialResidual(),
        pLinearSolvePerformance.finalResidual(),
        pLinearSolvePerformance.nIterations()
    );

    if (simple.finalNonOrthogonalIter())
    {
        phi = phiHbyA - pEqn.flux();
    }
}

#include "continuityErrs.H"

// Explicitly relax pressure for momentum corrector
p.relax();

// Momentum corrector
U = HbyA - rAtU()*fvc::grad(p);
U.correctBoundaryConditions();
fvOptions.correct(U);

turbulence->correct();

// --- Turbulence Brinkman Penalization Guard ---
// Prevent epsilon→0 SIGFPE in solid zones (kEpsilon: nut = Cmu*k²/epsilon).
// k and epsilon are floored to 1e-15 to unconditionally block division-by-zero
// across solid/fluid interface cells on fine meshes after aggressive MMA steps.
// This implements the turbulence suppression strategy documented in PlanTurbulenceOpt.md.
{
    volScalarField& k_rans =
        const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("k"));
    volScalarField& epsilon_rans =
        const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("epsilon"));
    forAll(k_rans, cellI)
    {
        k_rans[cellI]       = Foam::max(k_rans[cellI],       scalar(1e-15));
        epsilon_rans[cellI] = Foam::max(epsilon_rans[cellI], scalar(1e-15));
    }
    k_rans.correctBoundaryConditions();
    epsilon_rans.correctBoundaryConditions();
}
// ----------------------------------------------

U.storePrevIter();
p.storePrevIter();
phi.storePrevIter();

}

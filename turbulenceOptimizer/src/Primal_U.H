for(i=0;i<primalFlowIterations;i++)
{
MRF.correctBoundaryVelocity(U);
if (useTurbulentResistance)
{
  volScalarField nutByNu(turbulence->nut()/Foam::max(nu.value(), scalar(SMALL)));
  forAll(nutByNu, cellI)
  {
    const scalar clippedNutByNu =
      Foam::min(Foam::max(nutByNu[cellI], scalar(0.0)), turbulenceResistanceCap);
    alphaResistanceScale[cellI] = 1.0 + turbulenceResistanceFactor*clippedNutByNu;
  }
}
else
{
  alphaResistanceScale = scalar(1.0);
}
alphaResistanceScale.correctBoundaryConditions();
volScalarField alphaEff(alpha*alphaResistanceScale);

tmp<fvVectorMatrix> tUEqn
(
    fvm::div(phi, U)
  + MRF.DDt(U)
  + turbulence->divDevReff(U)
  //- fvm::laplacian(nu,U)
  + fvm::Sp(alphaEff, U)
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

if (((i+1) % turbulenceCorrectorInterval == 0) || (i+1 == primalFlowIterations))
{
  turbulence->correct();
}

U.storePrevIter();
p.storePrevIter();
phi.storePrevIter();

}

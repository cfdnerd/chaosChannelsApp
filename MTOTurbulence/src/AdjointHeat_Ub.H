for(i=0;i<40;i++)
{    
volVectorField adjointTransposeConvectionb((fvc::grad(U) & Ub));

tmp<fvVectorMatrix> tUbEqn
(
    fvm::div(-phi, Ub)
  + adjointTransposeConvectionb
  - fvm::laplacian(turbulence->nuEff(),Ub)
  + fvm::Sp(alpha, Ub)
  + b1*Tb*fvc::grad(T)
 ==
    fvOptions(Ub)
);
fvVectorMatrix& UbEqn = tUbEqn.ref();

UbEqn.relax();

fvOptions.constrain(UbEqn);

const auto UbLinearSolvePerformance = solve(UbEqn == -fvc::grad(pb));
recordLinearSolve
(
    UbLinearSolveStats,
    cmptMax(UbLinearSolvePerformance.initialResidual()),
    cmptMax(UbLinearSolvePerformance.finalResidual()),
    cmptMax(UbLinearSolvePerformance.nIterations())
);

fvOptions.correct(Ub);
//****************************************
volScalarField rAUb(1.0/UbEqn.A());
volVectorField HbyAb("HbyAb", Ub);
HbyAb = rAUb*UbEqn.H();
tUbEqn.clear();
surfaceScalarField phiHbyAb("phiHbyAb", fvc::flux(HbyAb));
adjustPhi(phiHbyAb, Ub, pb);

tmp<volScalarField> rAtUb(rAUb);

 if (simple.consistent())
{
    rAtUb = 1.0/(1.0/rAUb - UbEqn.H1());
    phiHbyAb +=
        fvc::interpolate(rAtUb() - rAUb)*fvc::snGrad(pb)*mesh.magSf();
    HbyAb -= (rAUb - rAtUb())*fvc::grad(pb);
}
tUbEqn.clear();
// Update the pressure BCs to ensure flux consistency
constrainPressure(pb, Ub, phiHbyAb, rAtUb(), MRF);  
// Non-orthogonal pressure corrector loop
while (simple.correctNonOrthogonal())
{
    fvScalarMatrix pbEqn
    (
        fvm::laplacian(rAUb, pb) == fvc::div(phiHbyAb)
    );

    pbEqn.setReference(pbRefCell, pbRefValue);
    const auto pbLinearSolvePerformance = pbEqn.solve();
    recordLinearSolve
    (
        pbLinearSolveStats,
        pbLinearSolvePerformance.initialResidual(),
        pbLinearSolvePerformance.finalResidual(),
        pbLinearSolvePerformance.nIterations()
    );

    if (simple.finalNonOrthogonalIter())
    {
        phib = phiHbyAb - pbEqn.flux();
    }
}

//#include "adjointContinuityErrs.H"

// Explicitly relax pressure for adjoint momentum corrector
pb.relax();

// Adjoint momentum corrector
Ub = HbyAb - rAUb*fvc::grad(pb);
Ub.correctBoundaryConditions();

forAll(Ub, cellI)
{
    Ub[cellI].x() = Foam::max(Foam::min(Ub[cellI].x(), 1e8), -1e8);
    Ub[cellI].y() = Foam::max(Foam::min(Ub[cellI].y(), 1e8), -1e8);
    Ub[cellI].z() = Foam::max(Foam::min(Ub[cellI].z(), 1e8), -1e8);
}
Ub.correctBoundaryConditions();

fvOptions.correct(Ub);

Ub.storePrevIter();
pb.storePrevIter();
phib.storePrevIter();
}

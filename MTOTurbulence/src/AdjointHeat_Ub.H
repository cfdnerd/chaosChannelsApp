for(i=0;i<40;i++)
{    
  volTensorField gradU = fvc::grad(U);
  forAll(gradU, cellI)
  {
      gradU[cellI].xx() = Foam::max(Foam::min(gradU[cellI].xx(), 5000.0), -5000.0);
      gradU[cellI].xy() = Foam::max(Foam::min(gradU[cellI].xy(), 5000.0), -5000.0);
      gradU[cellI].xz() = Foam::max(Foam::min(gradU[cellI].xz(), 5000.0), -5000.0);
      gradU[cellI].yx() = Foam::max(Foam::min(gradU[cellI].yx(), 5000.0), -5000.0);
      gradU[cellI].yy() = Foam::max(Foam::min(gradU[cellI].yy(), 5000.0), -5000.0);
      gradU[cellI].yz() = Foam::max(Foam::min(gradU[cellI].yz(), 5000.0), -5000.0);
      gradU[cellI].zx() = Foam::max(Foam::min(gradU[cellI].zx(), 5000.0), -5000.0);
      gradU[cellI].zy() = Foam::max(Foam::min(gradU[cellI].zy(), 5000.0), -5000.0);
      gradU[cellI].zz() = Foam::max(Foam::min(gradU[cellI].zz(), 5000.0), -5000.0);
  }
  volVectorField adjointTransposeConvectionb((gradU & Ub));

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

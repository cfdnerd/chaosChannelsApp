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
    alphaResistanceScale[cellI] =
      1.0 + turbulenceResistanceFactor*turbulenceStatePenaltyFactor*clippedNutByNu;
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

const bool doTurbulenceCorrection =
  opt > freezeTurbulenceUntilIteration
  &&
  (((i+1) % turbulenceCorrectorInterval == 0) || (i+1 == primalFlowIterations));
if (doTurbulenceCorrection)
{
  const volScalarField kBefore(turbulence->k());
  const volScalarField epsilonBefore(turbulence->epsilon());
  const volScalarField nutBefore(turbulence->nut());

  turbulence->correct();

  if
  (
    enforceTurbulenceBounds
    &&
    mesh.foundObject<volScalarField>("k")
    &&
    mesh.foundObject<volScalarField>("epsilon")
    &&
    mesh.foundObject<volScalarField>("nut")
  )
  {
    volScalarField& kField =
      const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("k"));
    volScalarField& epsilonField =
      const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("epsilon"));
    volScalarField& nutField =
      const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("nut"));
    scalarField& kInternal = kField.primitiveFieldRef();
    scalarField& epsilonInternal = epsilonField.primitiveFieldRef();
    scalarField& nutInternal = nutField.primitiveFieldRef();
    forAll(kInternal, cellI)
    {
      const scalar originalK = kInternal[cellI];
      const scalar originalEpsilon = epsilonInternal[cellI];
      const scalar originalNut = nutInternal[cellI];

      if (!std::isfinite(static_cast<double>(originalK)))
      {
        kInternal[cellI] = kMinBound;
      }
      else
      {
        kInternal[cellI] = Foam::min(Foam::max(originalK, kMinBound), kMaxBound);
      }
      if (!std::isfinite(static_cast<double>(originalEpsilon)))
      {
        epsilonInternal[cellI] = epsilonMinBound;
      }
      else
      {
        epsilonInternal[cellI] =
          Foam::min(Foam::max(originalEpsilon, epsilonMinBound), epsilonMaxBound);
      }
      if (!std::isfinite(static_cast<double>(originalNut)))
      {
        nutInternal[cellI] = 0.0;
      }
      else
      {
        nutInternal[cellI] = Foam::min(Foam::max(originalNut, scalar(0.0)), nutMaxBound);
      }

      if (kInternal[cellI] != originalK)
      {
        turbulenceKClipCount++;
      }
      if (epsilonInternal[cellI] != originalEpsilon)
      {
        turbulenceEpsilonClipCount++;
      }
      if (nutInternal[cellI] != originalNut)
      {
        turbulenceNutClipCount++;
      }
    }

    kField.correctBoundaryConditions();
    epsilonField.correctBoundaryConditions();
    nutField.correctBoundaryConditions();
  }

  const volScalarField kAfter(turbulence->k());
  const volScalarField epsilonAfter(turbulence->epsilon());
  const volScalarField nutAfter(turbulence->nut());

  auto accumulateCorrectionStats =
  [&](const volScalarField& fieldBefore, const volScalarField& fieldAfter, FieldCorrectionStats& stats)
  {
    scalar absDeltaSumLocal(0.0);
    scalar absDeltaMaxLocal(0.0);
    scalar relL2NumeratorLocal(0.0);
    scalar relL2DenominatorLocal(0.0);
    const scalarField& beforeInternal = fieldBefore.primitiveField();
    const scalarField& afterInternal = fieldAfter.primitiveField();
    forAll(afterInternal, cellI)
    {
      const scalar delta = Foam::mag(afterInternal[cellI] - beforeInternal[cellI]);
      absDeltaSumLocal += delta;
      absDeltaMaxLocal = Foam::max(absDeltaMaxLocal, delta);
      relL2NumeratorLocal += delta*delta;
      relL2DenominatorLocal += beforeInternal[cellI]*beforeInternal[cellI];
    }

    scalar globalCellCount = afterInternal.size();
    reduce(absDeltaSumLocal, sumOp<scalar>());
    reduce(absDeltaMaxLocal, maxOp<scalar>());
    reduce(relL2NumeratorLocal, sumOp<scalar>());
    reduce(relL2DenominatorLocal, sumOp<scalar>());
    reduce(globalCellCount, sumOp<scalar>());

    const scalar meanAbsDelta =
      absDeltaSumLocal/Foam::max(globalCellCount, scalar(SMALL));
    const scalar relL2Delta =
      Foam::sqrt(relL2NumeratorLocal/Foam::max(relL2DenominatorLocal, scalar(SMALL)));

    recordFieldCorrection(stats, meanAbsDelta, absDeltaMaxLocal, relL2Delta);
  };

  accumulateCorrectionStats(kBefore, kAfter, kTurbulenceCorrectionStats);
  accumulateCorrectionStats(epsilonBefore, epsilonAfter, epsilonTurbulenceCorrectionStats);
  accumulateCorrectionStats(nutBefore, nutAfter, nutTurbulenceCorrectionStats);
}

U.storePrevIter();
p.storePrevIter();
phi.storePrevIter();

}

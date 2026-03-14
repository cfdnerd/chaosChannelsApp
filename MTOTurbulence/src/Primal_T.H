T=T2;
volScalarField weight = (xh*(1+qu)/(qu+xh));

// Sanitize nut references and calculate effective diffusivity
const volScalarField& nut = turbulence->nut();
volScalarField DTeff("DTeff", DT + weight*nut/Prt);

// Guard DTeff against non-finite values before solving
forAll(DTeff, cellI)
{
    if (!std::isfinite(DTeff[cellI]) || DTeff[cellI] > 10.0)
    {
        DTeff[cellI] = 1e-6; // Reset to safe fluid-like value if unstable
    }
}
DTeff.correctBoundaryConditions();
fvScalarMatrix TEqn
(
    fvm::div(phi,T)
   -fvm::laplacian(DTeff, T)
);
if(DDHS)
{
 TEqn-=thermalSource*(1-xh);
}
else
{
 TEqn-=thermalSource;
}
const auto TLinearSolvePerformance = TEqn.solve();
recordLinearSolve
(
  TLinearSolveStats,
  TLinearSolvePerformance.initialResidual(),
  TLinearSolvePerformance.finalResidual(),
  TLinearSolvePerformance.nIterations()
);

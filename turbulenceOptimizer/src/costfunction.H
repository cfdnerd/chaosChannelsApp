Info<<"cost function"<<endl;
MeanT=0;
PowerDiss=0;
MeanTPhysical=0;
PowerDissPhysical=0;
maxT=-GREAT;
V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));

scalar meanTCase = fvc::domainIntegrate(T).value()/area;
scalar powerDissCase = 0;
for (i=0; i<nObjPatch; i++)
{
     powerDissCase =powerDissCase - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]])));
}
powerDissCase=powerDissCase/PowerDiss0;
reduce(powerDissCase, sumOp<scalar>());
scalar maxTCase = gMax(T);
reduce(maxTCase, maxOp<scalar>());

if
(
  turbulenceProMode
  &&
  (
    !std::isfinite(static_cast<double>(meanTCase))
    ||
    !std::isfinite(static_cast<double>(powerDissCase))
    ||
    !std::isfinite(static_cast<double>(maxTCase))
  )
)
{
  stopOptimization = true;
  optimizerStopReason = "nonfinite_cost_metrics";
  if (Pstream::master())
  {
    WarningInFunction
      << "Detected non-finite cost metrics at optimizer iteration " << opt
      << ". MeanT=" << meanTCase
      << ", Power=" << powerDissCase
      << ", MaxT=" << maxTCase
      << endl;
  }
}

MeanT = meanTCase;
PowerDiss = powerDissCase;
MeanTPhysical = meanTCase;
PowerDissPhysical = powerDissCase;
maxT = maxTCase;

if (normalizeScenarioObjectives)
{
  if (!objectiveIterationOneReferenceCaptured)
  {
    objectiveIterationOneReference = Foam::max(Foam::mag(meanTCase), scalar(SMALL));
    objectiveIterationOneReferenceCaptured = true;
    if (Pstream::master())
    {
      Info<< "Captured objective normalization reference at iteration " << opt
          << ": " << objectiveIterationOneReference << endl;
    }
  }
  if (!powerIterationOneReferenceCaptured)
  {
    powerIterationOneReference = Foam::max(Foam::mag(powerDissCase), scalar(SMALL));
    powerIterationOneReferenceCaptured = true;
    if (Pstream::master())
    {
      Info<< "Captured power normalization reference at iteration " << opt
          << ": " << powerIterationOneReference << endl;
    }
  }
  const scalar objectiveReference =
    objectiveIterationOneReferenceCaptured
    ? objectiveIterationOneReference
    : Foam::max(Foam::mag(meanTCase), scalar(SMALL));
  const scalar powerReference =
    powerIterationOneReferenceCaptured
    ? powerIterationOneReference
    : Foam::max(Foam::mag(powerDissCase), scalar(SMALL));
  MeanT = meanTCase/Foam::max(objectiveReference, scalar(SMALL));
  PowerDiss = powerDissCase/Foam::max(powerReference, scalar(SMALL));
}

objectiveNormalized = MeanT;
powerConstraintNormalized = PowerDiss;
PD0Sug = PowerDiss0;
if (historyWritePD0Sug)
{
  if
  (
    powerConstraintNormalized > PowerDissRelax
    &&
    PowerDissMax > SMALL
  )
  {
    PD0Sug =
      PowerDiss0
     *(
        powerConstraintNormalized
       /Foam::max(PowerDissMax, scalar(SMALL))
      );
  }
}

PowerDissActualOverBase = PowerDissPhysical;
PowerDissActual = PowerDissActualOverBase*PowerDiss0;

if (Pstream::master())
{
   std::ofstream historyFile(optimizationHistoryFile, std::ios::app);
   historyFile << std::right
     << std::setw(8) << opt
     << std::scientific << std::setprecision(8)
     << std::setw(16) << PowerDissActual
     << std::setw(16) << powerConstraintNormalized;
   if (historyWritePD0Sug)
   {
     historyFile << std::setw(16) << PD0Sug;
   }
   historyFile
     << std::setw(16) << V
     << std::setw(16) << objectiveNormalized
     << std::setw(16) << maxT
     << '\n';
   historyFile.flush();
}

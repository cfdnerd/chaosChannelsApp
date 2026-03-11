Info<<"cost function"<<endl;
const scalar meshFactor = std::pow(mesh.V()[0],1.0/3);
if (scenarioI == 0)
{
  MeanT=0;
  PowerDiss=0;
  MeanTPhysical=0;
  PowerDissPhysical=0;
  maxT=-GREAT;
  V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));
}

scalar scenarioMeanT = fvc::domainIntegrate(T).value()/area;
scalar scenarioPowerDiss = 0;
for (i=0; i<nObjPatch; i++)
{
     scenarioPowerDiss =scenarioPowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]])));
}
scenarioPowerDiss=scenarioPowerDiss/PowerDiss0;
if(GeoDim==2)
{
 scenarioPowerDiss=scenarioPowerDiss/meshFactor;
}
reduce(scenarioPowerDiss, sumOp<scalar>());
scalar scenarioMaxT = gMax(T);
reduce(scenarioMaxT, maxOp<scalar>());

if
(
  turbulenceProMode
  &&
  (
    !std::isfinite(static_cast<double>(scenarioMeanT))
    ||
    !std::isfinite(static_cast<double>(scenarioPowerDiss))
    ||
    !std::isfinite(static_cast<double>(scenarioMaxT))
  )
)
{
  stopOptimization = true;
  optimizerStopReason = "nonfinite_cost_metrics";
  if (Pstream::master())
  {
    WarningInFunction
      << "Detected non-finite cost metrics at optimizer iteration " << opt
      << " (scenario " << scenarioI + 1 << ")."
      << " MeanT=" << scenarioMeanT
      << ", Power=" << scenarioPowerDiss
      << ", MaxT=" << scenarioMaxT
      << endl;
  }
}

scenarioMeanTLast[scenarioI] = scenarioMeanT;
scenarioPowerLast[scenarioI] = scenarioPowerDiss;
scenarioMaxTLast[scenarioI] = scenarioMaxT;

if (!normalizeScenarioObjectives)
{
  scenarioThermalReference[scenarioI] = 1.0;
  scenarioPowerReference[scenarioI] = 1.0;
}
if (scenarioThermalReference[scenarioI] <= SMALL)
{
  scenarioThermalReference[scenarioI] = Foam::max(Foam::mag(scenarioMeanT), scalar(SMALL));
}
if (scenarioPowerReference[scenarioI] <= SMALL)
{
  scenarioPowerReference[scenarioI] = Foam::max(Foam::mag(scenarioPowerDiss), scalar(SMALL));
}

const scalar scenarioMeanTNormalized =
  scenarioMeanT/Foam::max(scenarioThermalReference[scenarioI], scalar(SMALL));
const scalar scenarioPowerNormalized =
  scenarioPowerDiss/Foam::max(scenarioPowerReference[scenarioI], scalar(SMALL));
scenarioObjectiveNormalizedLast[scenarioI] = scenarioMeanTNormalized;
scenarioPowerNormalizedLast[scenarioI] = scenarioPowerNormalized;

const scalar scenarioWeight = scenarioWeights[scenarioI];
MeanT += scenarioWeight*scenarioMeanTNormalized;
PowerDiss += scenarioWeight*scenarioPowerNormalized;
MeanTPhysical += scenarioWeight*scenarioMeanT;
PowerDissPhysical += scenarioWeight*scenarioPowerDiss;
maxT = Foam::max(maxT, scenarioMaxT);

if (scenarioI == nVelocityScenarios-1)
{
  objectiveNormalized = MeanT;
  powerConstraintNormalized = PowerDiss;

  PowerDissActualOverBase = PowerDissPhysical;
  if(GeoDim==2)
  {
    PowerDissActualOverBase = PowerDissActualOverBase*meshFactor;
  }
  PowerDissActual = PowerDissActualOverBase*PowerDiss0;

  if (Pstream::master())
  {
     std::ofstream historyFile(optimizationHistoryFile, std::ios::app);
     historyFile
       << std::right
       << std::setw(8) << opt
       << std::scientific << std::setprecision(8)
       << std::setw(16) << PowerDissActual
       << std::setw(16) << powerConstraintNormalized
       << std::setw(16) << V
       << std::setw(16) << objectiveNormalized
       << std::setw(16) << maxT
       << '\n';
     historyFile.flush();
  }
}

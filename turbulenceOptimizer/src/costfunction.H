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
  PDActLog = PowerDissPhysical*PowerDiss0;
  PD0Sug = PowerDiss0;
  if (powerConstraintNormalized > PowerDissRelax && PowerDissMax > SMALL)
  {
    PD0Sug = PowerDiss0*(powerConstraintNormalized/Foam::max(PowerDissMax, SMALL));
  }
  const bool pd0SugDiffObserved = Foam::mag(PD0Sug - PowerDiss0) > SMALL;
  if (pd0SugDiffObserved)
  {
    PD0SugMaxCum = Foam::max(PD0SugMaxCum, PD0Sug);
    PD0SugSumCum += PD0Sug;
    PD0SugCountCum++;
  }

  if (Pstream::master())
  {
     if (pd0SugDiffObserved && !PD0SugLogInitialized)
     {
       std::ofstream PD0SugLog(suggestedPowerDiss0File, std::ios::out);
       PD0SugLog
         << "# MTO_ThermalFluid suggested PowerDiss0 summary\n"
         << "# Cols: Iter PD0SugMax_1_to_Iter PD0SugAvg_1_to_Iter\n"
         << "# Stats exclude unchanged values where PD0Sug == PowerDiss0\n"
         << std::left
         << std::setw(8) << "Iter"
         << std::setw(20) << "PD0SugMax_1_to_Iter"
         << std::setw(20) << "PD0SugAvg_1_to_Iter"
         << '\n';
       PD0SugLogInitialized = true;
     }
     std::ofstream historyFile(optimizationHistoryFile, std::ios::app);
     if (opt > 1 && ((opt - 1) % 25 == 0))
     {
       historyFile
         << std::left
         << std::setw(8) << "Iter"
         << std::setw(16) << "PDActLog"
         << std::setw(16) << "PDCon"
         << std::setw(16) << "PD0Sug"
         << std::setw(16) << "V"
         << std::setw(16) << "MeanT"
         << std::setw(16) << "MaxT"
         << '\n';
     }
     historyFile
       << std::right
       << std::setw(8) << opt
       << std::scientific << std::setprecision(8)
       << std::setw(16) << PDActLog
       << std::setw(16) << powerConstraintNormalized
       << std::setw(16) << PD0Sug
       << std::setw(16) << V
       << std::setw(16) << objectiveNormalized
       << std::setw(16) << maxT
       << '\n';
     if (PD0SugLogInitialized && PD0SugCountCum > 0 && opt > 0 && (opt % 100 == 0))
     {
       const scalar PD0SugAvgCum =
         PD0SugSumCum/Foam::max(static_cast<scalar>(PD0SugCountCum), SMALL);
       std::ofstream PD0SugLog(suggestedPowerDiss0File, std::ios::app);
       PD0SugLog
         << std::right
         << std::setw(8) << opt
         << std::scientific << std::setprecision(8)
         << std::setw(20) << PD0SugMaxCum
         << std::setw(20) << PD0SugAvgCum
         << '\n';
     }
  }
}

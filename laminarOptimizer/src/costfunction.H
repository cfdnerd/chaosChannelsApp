Info<<"cost function"<<endl;
MeanT=fvc::domainIntegrate(T).value()/area;
V=fvc::domainIntegrate(xh).value()/area-(voluse+0*Foam::max(0,0.5-opt*0.015));  

PowerDiss=0;
for (i=0; i<nObjPatch; i++)
{
     PowerDiss =PowerDiss - sum(
     phi.boundaryField()[conPatchList[i]]*(p.boundaryField()[conPatchList[i]]
     + 0.5*magSqr(U.boundaryField()[conPatchList[i]]))); 
}
PowerDiss=PowerDiss/PowerDiss0;
if(GeoDim==2)
{
 PowerDiss=PowerDiss/std::pow(mesh.V()[0],1.0/3);
}
reduce(PowerDiss, sumOp<scalar>());

const scalar meshFactor = std::pow(mesh.V()[0],1.0/3);
PowerDissActualOverBase = PowerDiss;
if(GeoDim==2)
{
  PowerDissActualOverBase = PowerDissActualOverBase*meshFactor;
}
PowerDissActual = PowerDissActualOverBase*PowerDiss0;
const scalar PDAct = PowerDiss*PowerDiss0*meshFactor;
PDActLog = PDAct/meshFactor;
PD0Sug = PowerDiss0;
if (PowerDiss > PowerDissRelax && PowerDissMax > SMALL)
{
  PD0Sug = PowerDiss0*(PowerDiss/Foam::max(PowerDissMax, SMALL));
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
     << std::setw(16) << PowerDiss
     << std::setw(16) << PD0Sug
     << std::setw(16) << V
     << std::setw(16) << MeanT
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

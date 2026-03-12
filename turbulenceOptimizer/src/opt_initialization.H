Info<< "Initialization\n" << endl;
const char* optimizationHistoryFile = "optimization.hst";
const char* debugOptimizerJsonFile = "debugOptimizer.jsonl";
const char* debugOptimizerLogFile = "debugOptimizer.log";
const char* solverConvergenceFile = "solverConvergences.log";
if (Pstream::master())
{
  std::remove(optimizationHistoryFile);
  std::remove(debugOptimizerJsonFile);
  std::remove(debugOptimizerLogFile);
  std::remove(solverConvergenceFile);
  std::ofstream historyFile(optimizationHistoryFile, std::ios::out);
  historyFile
    << "# MTO_ThermalFluid optimization history\n"
    << "# Cols: Iter PDAct PDCon PD0Sug V MeanT MaxT\n"
    << "# PDAct=PowerDissActual physical-scale pumping power metric\n"
    << "# PDCon=power constraint value used by MMA (normalized)\n"
    << "# PD0Sug=PowerDiss0*(PDCon/PowerDissMax) when PDCon>PowerDissRelax; else PowerDiss0\n";
  if (normalizeScenarioObjectives)
  {
    historyFile
      << "# MeanT is normalized by its iteration-1 absolute value reference\n";
  }
  historyFile
    << std::left
    << std::setw(8) << "Iter"
    << std::setw(16) << "PDAct"
    << std::setw(16) << "PDCon"
    << std::setw(16) << "PD0Sug"
    << std::setw(16) << "V"
    << std::setw(16) << "MeanT"
    << std::setw(16) << "MaxT"
    << '\n';

  std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::out);
  debugLogFile
    << "# MTO_ThermalFluid optimizer debug log\n"
    << "# Structured per-iteration diagnostics for optimizer development\n"
    << "# Format: sectioned human-readable report, one block per optimization iteration\n";

  std::ofstream solverConvergenceLogFile(solverConvergenceFile, std::ios::out);
  solverConvergenceLogFile
    << "# MTO_ThermalFluid solver convergence log\n"
    << "# Linear-solver statistics grouped by optimization iteration\n";
}

struct LinearSolveStats
{
  label callCount;
  scalar initialResidualSum;
  scalar initialResidualMax;
  scalar finalResidualSum;
  scalar finalResidualMax;
  label iterationSum;
  label iterationMax;
  scalar lastInitialResidual;
  scalar lastFinalResidual;
  label lastIterations;

  LinearSolveStats()
  :
    callCount(0),
    initialResidualSum(0),
    initialResidualMax(0),
    finalResidualSum(0),
    finalResidualMax(0),
    iterationSum(0),
    iterationMax(0),
    lastInitialResidual(0),
    lastFinalResidual(0),
    lastIterations(0)
  {}
};

struct FieldCorrectionStats
{
  label callCount;
  scalar meanAbsDeltaSum;
  scalar meanAbsDeltaMax;
  scalar maxAbsDeltaMax;
  scalar relL2DeltaSum;
  scalar relL2DeltaMax;
  scalar lastMeanAbsDelta;
  scalar lastMaxAbsDelta;
  scalar lastRelL2Delta;

  FieldCorrectionStats()
  :
    callCount(0),
    meanAbsDeltaSum(0),
    meanAbsDeltaMax(0),
    maxAbsDeltaMax(0),
    relL2DeltaSum(0),
    relL2DeltaMax(0),
    lastMeanAbsDelta(0),
    lastMaxAbsDelta(0),
    lastRelL2Delta(0)
  {}
};

auto resetLinearSolveStats =
[](LinearSolveStats& stats)
{
  stats.callCount = 0;
  stats.initialResidualSum = 0;
  stats.initialResidualMax = 0;
  stats.finalResidualSum = 0;
  stats.finalResidualMax = 0;
  stats.iterationSum = 0;
  stats.iterationMax = 0;
  stats.lastInitialResidual = 0;
  stats.lastFinalResidual = 0;
  stats.lastIterations = 0;
};

auto resetFieldCorrectionStats =
[](FieldCorrectionStats& stats)
{
  stats.callCount = 0;
  stats.meanAbsDeltaSum = 0;
  stats.meanAbsDeltaMax = 0;
  stats.maxAbsDeltaMax = 0;
  stats.relL2DeltaSum = 0;
  stats.relL2DeltaMax = 0;
  stats.lastMeanAbsDelta = 0;
  stats.lastMaxAbsDelta = 0;
  stats.lastRelL2Delta = 0;
};

auto recordLinearSolve =
[](LinearSolveStats& stats, const scalar initialResidual, const scalar finalResidual, const label iterations)
{
  stats.callCount++;
  stats.initialResidualSum += initialResidual;
  stats.initialResidualMax = Foam::max(stats.initialResidualMax, initialResidual);
  stats.finalResidualSum += finalResidual;
  stats.finalResidualMax = Foam::max(stats.finalResidualMax, finalResidual);
  stats.iterationSum += iterations;
  stats.iterationMax = Foam::max(stats.iterationMax, iterations);
  stats.lastInitialResidual = initialResidual;
  stats.lastFinalResidual = finalResidual;
  stats.lastIterations = iterations;
};

auto recordFieldCorrection =
[](FieldCorrectionStats& stats, const scalar meanAbsDelta, const scalar maxAbsDelta, const scalar relL2Delta)
{
  stats.callCount++;
  stats.meanAbsDeltaSum += meanAbsDelta;
  stats.meanAbsDeltaMax = Foam::max(stats.meanAbsDeltaMax, meanAbsDelta);
  stats.maxAbsDeltaMax = Foam::max(stats.maxAbsDeltaMax, maxAbsDelta);
  stats.relL2DeltaSum += relL2Delta;
  stats.relL2DeltaMax = Foam::max(stats.relL2DeltaMax, relL2Delta);
  stats.lastMeanAbsDelta = meanAbsDelta;
  stats.lastMaxAbsDelta = maxAbsDelta;
  stats.lastRelL2Delta = relL2Delta;
};

const scalar linearSolveObsResidual = 1e-6;
const scalar linearSolveWarnResidual = 1e-4;

auto classifyLinearSolveStatus =
[&](const LinearSolveStats& stats) -> const char*
{
  if (stats.callCount == 0)
  {
    return "SKIP";
  }
  if (stats.finalResidualMax > linearSolveWarnResidual)
  {
    return "WARN";
  }
  if (stats.finalResidualMax > linearSolveObsResidual)
  {
    return "OBS";
  }
  return "OK";
};

auto linearSolveSeverity =
[](const char* status) -> label
{
  if (std::strcmp(status, "WARN") == 0)
  {
    return 2;
  }
  if (std::strcmp(status, "OBS") == 0)
  {
    return 1;
  }
  return 0;
};

auto linearSolveSeverityFromStats =
[&](const LinearSolveStats& stats) -> label
{
  return linearSolveSeverity(classifyLinearSolveStatus(stats));
};

LinearSolveStats ULinearSolveStats;
LinearSolveStats pLinearSolveStats;
LinearSolveStats TLinearSolveStats;
LinearSolveStats TFallbackLinearSolveStats;
LinearSolveStats TbLinearSolveStats;
LinearSolveStats UbLinearSolveStats;
LinearSolveStats pbLinearSolveStats;
LinearSolveStats UaLinearSolveStats;
LinearSolveStats paLinearSolveStats;
LinearSolveStats xpFilterLinearSolveStats;
LinearSolveStats fsensFilterLinearSolveStats;
LinearSolveStats gsensPowerFilterLinearSolveStats;
LinearSolveStats gsensVolFilterLinearSolveStats;
FieldCorrectionStats kTurbulenceCorrectionStats;
FieldCorrectionStats epsilonTurbulenceCorrectionStats;
FieldCorrectionStats nutTurbulenceCorrectionStats;
label thermalFallbackTriggerCount(0);
label thermalFallbackSuccessCount(0);
label thermalFallbackMissingSolverCount(0);
label thermalFallbackFailureCount(0);
label thermalFallbackLastPrimaryIterations(0);
label thermalFallbackLastFallbackIterations(-1);
scalar thermalFallbackLastPrimaryFinalResidual(0.0);
scalar thermalFallbackLastFallbackFinalResidual(-1.0);
scalar thermalLastMin(0.0);
scalar thermalLastMax(0.0);
label turbulenceKClipCount(0);
label turbulenceEpsilonClipCount(0);
label turbulenceNutClipCount(0);

int i=0,opt=0;
double del=0;
scalar n(0),N(0),eta0(0),eta5(0.5),eta1(1),y0(0),y5(0);
n=x.size();//n:length of x in the i-th core
N=x.size();//N:design freedom
reduce(N, sumOp<scalar>()); 

scalar V(0),area=gSum(mesh.V()),MeanT(0),PowerDiss(0),PowerDissActual(0),PowerDissActualOverBase(0),PD0Sug(0),maxT(0),sumT(0);
scalar MeanTPhysical(0),PowerDissPhysical(0);
scalar objectiveNormalized(0),powerConstraintNormalized(0);
scalar objectiveIterationOneReference(1.0);
bool objectiveIterationOneReferenceCaptured(false);
scalar previousObjectiveValue(0);
bool hasPreviousObjectiveValue(false);
bool stopOptimization(false);
label optimizerConvergedStreak(0);
word optimizerStopReason("running");
double len=std::pow(area/N, 1.0/3);
b.value()=1.0/(filterRStart*len/3.464)/(filterRStart*len/3.464);
//***********************************//MMA
int m=2;//number of constraint
MMA mma(n,m);
mma.asyminit = mma_init;
mma.asymdec = mma_dec;
mma.asyminc = mma_inc;
mma.raa0 = raa0;
mma.movelimit = movlim;
std::vector<double> xmma(n);
std::vector<double> dfdx(n);
std::vector<double> gx(m);
std::vector<std::vector<double>> dgdx(m);
for (int i = 0; i < m; i++)
{
  dgdx[i].resize(n);
}

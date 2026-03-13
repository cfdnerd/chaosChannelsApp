scalar xLocalMin = 1e300, xLocalMax = -1e300;
scalar xpLocalMin = 1e300, xpLocalMax = -1e300;
scalar xhLocalMin = 1e300, xhLocalMax = -1e300;
scalar alphaLocalMin = 1e300, alphaLocalMax = -1e300;
scalar DTLocalMin = 1e300, DTLocalMax = -1e300;
scalar TLocalMin = 1e300, TLocalMax = -1e300;
scalar pLocalMin = 1e300, pLocalMax = -1e300;
scalar paLocalMin = 1e300, paLocalMax = -1e300;
scalar pbLocalMin = 1e300, pbLocalMax = -1e300;
scalar xVolumeSumLocal = 0;
scalar xhVolumeSumLocal = 0;
scalar grayVolumeLocal = 0;
scalar solidVolumeLocal = 0;
scalar fluidVolumeLocal = 0;
scalar maxAbsXChangeLocal = 0;
scalar l2XChangeLocal = 0;
scalar maxAbsXhChangeLocal = 0;
scalar l2XhChangeLocal = 0;
label xLowerBoundCountLocal = 0;
label xUpperBoundCountLocal = 0;
label xMoveLimitHitCountLocal = 0;
label xNonFiniteCountLocal = 0;
scalar drhoLocalMin = 1e300, drhoLocalMax = -1e300;
scalar fsensLocalMin = 1e300, fsensLocalMax = -1e300, fsensL2Local = 0;
scalar gsensPowerLocalMin = 1e300, gsensPowerLocalMax = -1e300, gsensPowerL2Local = 0;
scalar gsensVolLocalMin = 1e300, gsensVolLocalMax = -1e300, gsensVolL2Local = 0;
label fsensNonFiniteCountLocal = 0;
label gsensPowerNonFiniteCountLocal = 0;
label gsensVolNonFiniteCountLocal = 0;
scalar UMaxMagLocal = 0, UaMaxMagLocal = 0, UbMaxMagLocal = 0;
scalar UMeanMagLocal = 0, UaMeanMagLocal = 0, UbMeanMagLocal = 0;
scalar objectivePnormSumLocal = 0;

for(i=0;i<n;i++)
{
  const scalar cellVolume = mesh.V()[i];
  const scalar xValue = x[i];
  const scalar xpValue = xp[i];
  const scalar xhValue = xh[i];
  const scalar alphaValue = alpha[i];
  const scalar DTValue = DT[i];
  const scalar TValue = T[i];
  const scalar pValue = p[i];
  const scalar paValue = pa[i];
  const scalar pbValue = pb[i];
  const scalar xChange = xValue - xOld[i];
  const scalar xhChange = xhValue - xhOld[i];
  const scalar absXChange = Foam::mag(xChange);
  const scalar absXhChange = Foam::mag(xhChange);
  const scalar UMag = Foam::mag(U[i]);
  const scalar UaMag = Foam::mag(Ua[i]);
  const scalar UbMag = Foam::mag(Ub[i]);
  const scalar drhoValue = drho[i];
  const scalar fsensValue = fsensMeanT[i];
  const scalar gsensPowerValue = gsensPowerDiss[i];
  const scalar gsensVolValue = gsensVol[i];

  xLocalMin = Foam::min(xLocalMin, xValue);
  xLocalMax = Foam::max(xLocalMax, xValue);
  xpLocalMin = Foam::min(xpLocalMin, xpValue);
  xpLocalMax = Foam::max(xpLocalMax, xpValue);
  xhLocalMin = Foam::min(xhLocalMin, xhValue);
  xhLocalMax = Foam::max(xhLocalMax, xhValue);
  alphaLocalMin = Foam::min(alphaLocalMin, alphaValue);
  alphaLocalMax = Foam::max(alphaLocalMax, alphaValue);
  DTLocalMin = Foam::min(DTLocalMin, DTValue);
  DTLocalMax = Foam::max(DTLocalMax, DTValue);
  TLocalMin = Foam::min(TLocalMin, TValue);
  TLocalMax = Foam::max(TLocalMax, TValue);
  pLocalMin = Foam::min(pLocalMin, pValue);
  pLocalMax = Foam::max(pLocalMax, pValue);
  paLocalMin = Foam::min(paLocalMin, paValue);
  paLocalMax = Foam::max(paLocalMax, paValue);
  pbLocalMin = Foam::min(pbLocalMin, pbValue);
  pbLocalMax = Foam::max(pbLocalMax, pbValue);
  UMaxMagLocal = Foam::max(UMaxMagLocal, UMag);
  UaMaxMagLocal = Foam::max(UaMaxMagLocal, UaMag);
  UbMaxMagLocal = Foam::max(UbMaxMagLocal, UbMag);

  xVolumeSumLocal += xValue*cellVolume;
  xhVolumeSumLocal += xhValue*cellVolume;
  UMeanMagLocal += UMag*cellVolume;
  UaMeanMagLocal += UaMag*cellVolume;
  UbMeanMagLocal += UbMag*cellVolume;

  if (xhValue > 0.05 && xhValue < 0.95)
  {
    grayVolumeLocal += cellVolume;
  }
  if (xhValue <= 0.05)
  {
    solidVolumeLocal += cellVolume;
  }
  if (xhValue >= 0.95)
  {
    fluidVolumeLocal += cellVolume;
  }

  maxAbsXChangeLocal = Foam::max(maxAbsXChangeLocal, absXChange);
  l2XChangeLocal += xChange*xChange;
  maxAbsXhChangeLocal = Foam::max(maxAbsXhChangeLocal, absXhChange);
  l2XhChangeLocal += xhChange*xhChange;

  if (xValue <= 1e-8)
  {
    xLowerBoundCountLocal++;
  }
  if (xValue >= 1.0 - 1e-8)
  {
    xUpperBoundCountLocal++;
  }
  if (absXChange >= 0.999*movlim)
  {
    xMoveLimitHitCountLocal++;
  }
  if (!std::isfinite(static_cast<double>(xValue)))
  {
    xNonFiniteCountLocal++;
  }

  drhoLocalMin = Foam::min(drhoLocalMin, drhoValue);
  drhoLocalMax = Foam::max(drhoLocalMax, drhoValue);

  fsensLocalMin = Foam::min(fsensLocalMin, fsensValue);
  fsensLocalMax = Foam::max(fsensLocalMax, fsensValue);
  fsensL2Local += fsensValue*fsensValue;
  if (!std::isfinite(static_cast<double>(fsensValue)))
  {
    fsensNonFiniteCountLocal++;
  }

  gsensPowerLocalMin = Foam::min(gsensPowerLocalMin, gsensPowerValue);
  gsensPowerLocalMax = Foam::max(gsensPowerLocalMax, gsensPowerValue);
  gsensPowerL2Local += gsensPowerValue*gsensPowerValue;
  if (!std::isfinite(static_cast<double>(gsensPowerValue)))
  {
    gsensPowerNonFiniteCountLocal++;
  }

  gsensVolLocalMin = Foam::min(gsensVolLocalMin, gsensVolValue);
  gsensVolLocalMax = Foam::max(gsensVolLocalMax, gsensVolValue);
  gsensVolL2Local += gsensVolValue*gsensVolValue;
  if (!std::isfinite(static_cast<double>(gsensVolValue)))
  {
    gsensVolNonFiniteCountLocal++;
  }

  if (objFunction == 2 && maxT > SMALL)
  {
    objectivePnormSumLocal += std::pow(TValue/maxT, Pnorm);
  }
}

reduce(xLocalMin, minOp<scalar>());
reduce(xLocalMax, maxOp<scalar>());
reduce(xpLocalMin, minOp<scalar>());
reduce(xpLocalMax, maxOp<scalar>());
reduce(xhLocalMin, minOp<scalar>());
reduce(xhLocalMax, maxOp<scalar>());
reduce(alphaLocalMin, minOp<scalar>());
reduce(alphaLocalMax, maxOp<scalar>());
reduce(DTLocalMin, minOp<scalar>());
reduce(DTLocalMax, maxOp<scalar>());
reduce(TLocalMin, minOp<scalar>());
reduce(TLocalMax, maxOp<scalar>());
reduce(pLocalMin, minOp<scalar>());
reduce(pLocalMax, maxOp<scalar>());
reduce(paLocalMin, minOp<scalar>());
reduce(paLocalMax, maxOp<scalar>());
reduce(pbLocalMin, minOp<scalar>());
reduce(pbLocalMax, maxOp<scalar>());
reduce(UMaxMagLocal, maxOp<scalar>());
reduce(UaMaxMagLocal, maxOp<scalar>());
reduce(UbMaxMagLocal, maxOp<scalar>());
reduce(xVolumeSumLocal, sumOp<scalar>());
reduce(xhVolumeSumLocal, sumOp<scalar>());
reduce(grayVolumeLocal, sumOp<scalar>());
reduce(solidVolumeLocal, sumOp<scalar>());
reduce(fluidVolumeLocal, sumOp<scalar>());
reduce(UMeanMagLocal, sumOp<scalar>());
reduce(UaMeanMagLocal, sumOp<scalar>());
reduce(UbMeanMagLocal, sumOp<scalar>());
reduce(maxAbsXChangeLocal, maxOp<scalar>());
reduce(l2XChangeLocal, sumOp<scalar>());
reduce(maxAbsXhChangeLocal, maxOp<scalar>());
reduce(l2XhChangeLocal, sumOp<scalar>());
reduce(xLowerBoundCountLocal, sumOp<label>());
reduce(xUpperBoundCountLocal, sumOp<label>());
reduce(xMoveLimitHitCountLocal, sumOp<label>());
reduce(xNonFiniteCountLocal, sumOp<label>());
reduce(drhoLocalMin, minOp<scalar>());
reduce(drhoLocalMax, maxOp<scalar>());
reduce(fsensLocalMin, minOp<scalar>());
reduce(fsensLocalMax, maxOp<scalar>());
reduce(fsensL2Local, sumOp<scalar>());
reduce(fsensNonFiniteCountLocal, sumOp<label>());
reduce(gsensPowerLocalMin, minOp<scalar>());
reduce(gsensPowerLocalMax, maxOp<scalar>());
reduce(gsensPowerL2Local, sumOp<scalar>());
reduce(gsensPowerNonFiniteCountLocal, sumOp<label>());
reduce(gsensVolLocalMin, minOp<scalar>());
reduce(gsensVolLocalMax, maxOp<scalar>());
reduce(gsensVolL2Local, sumOp<scalar>());
reduce(gsensVolNonFiniteCountLocal, sumOp<label>());
reduce(objectivePnormSumLocal, sumOp<scalar>());

const scalar xMean = xVolumeSumLocal/area;
const scalar xhMean = xhVolumeSumLocal/area;
const scalar grayVolumeFraction = grayVolumeLocal/area;
const scalar solidVolumeFraction = solidVolumeLocal/area;
const scalar fluidVolumeFraction = fluidVolumeLocal/area;
const scalar l2XChange = std::sqrt(l2XChangeLocal);
const scalar l2XhChange = std::sqrt(l2XhChangeLocal);
const scalar normalizedL2XChange = l2XChange/std::sqrt(Foam::max(N, scalar(1)));
const scalar normalizedL2XhChange = l2XhChange/std::sqrt(Foam::max(N, scalar(1)));
const scalar UMeanMag = UMeanMagLocal/area;
const scalar UaMeanMag = UaMeanMagLocal/area;
const scalar UbMeanMag = UbMeanMagLocal/area;
const scalar fsensL2 = std::sqrt(fsensL2Local);
const scalar gsensPowerL2 = std::sqrt(gsensPowerL2Local);
const scalar gsensVolL2 = std::sqrt(gsensVolL2Local);

scalar objectiveValue = MeanT;
word objectiveLabel("MeanT");
if (objFunction == 2)
{
  objectiveLabel = "PnormT";
  if (maxT > SMALL)
  {
    objectiveValue = maxT*std::pow(objectivePnormSumLocal, 1.0/Pnorm);
  }
  else
  {
    objectiveValue = 0;
  }
}
const scalar objectiveDelta = hasPreviousObjectiveValue ? objectiveValue - previousObjectiveValue : 0;
previousObjectiveValue = objectiveValue;
hasPreviousObjectiveValue = true;

label solverSeverity(0);
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(ULinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(pLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(TLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(TbLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(UbLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(pbLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(UaLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(paLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(xpFilterLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(fsensFilterLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(gsensPowerFilterLinearSolveStats));
solverSeverity = Foam::max(solverSeverity, linearSolveSeverityFromStats(gsensVolFilterLinearSolveStats));
reduce(solverSeverity, maxOp<label>());

const scalar objectiveRelativeDelta = Foam::mag(objectiveDelta)/Foam::max(Foam::mag(objectiveValue), SMALL);
const scalar powerConstraintWeightSafe = Foam::max(optPowerConstraintWeight, SMALL);
const scalar powerConstraintMargin = gx[0]/powerConstraintWeightSafe;
const scalar volumeConstraintWeightSafe = Foam::max(optVolumeConstraintWeight, SMALL);
const scalar volumeConstraintMargin = gx[1]/volumeConstraintWeightSafe;
const scalar maxConstraintViolation = Foam::max(Foam::max(powerConstraintMargin, volumeConstraintMargin), scalar(0));
const bool convergenceCandidate =
  opt >= optMinIterations
  &&
  solverSeverity < 2
  &&
  maxConstraintViolation <= optConstraintTol
  &&
  objectiveRelativeDelta <= optObjectiveTol
  &&
  maxAbsXChangeLocal <= optStepTol
  &&
  maxAbsXhChangeLocal <= optStepTol;
optimizerConvergedStreak = convergenceCandidate ? optimizerConvergedStreak + 1 : 0;
stopOptimization = optimizerConvergedStreak >= optStallWindow;
optimizerStopReason = stopOptimization ? word("converged") : word("running");

if (Pstream::master())
{
  const scalar clockTime = runTime.elapsedClockTime();
  const char* designHealth =
    (xNonFiniteCountLocal > 0 || xLocalMin < -1e-8 || xLocalMax > 1.0 + 1e-8 || xhLocalMin < -1e-8 || xhLocalMax > 1.0 + 1e-8)
    ? "WARN"
    : (grayVolumeFraction > 0.50 ? "OBS" : "OK");
  const char* constraintHealth = (powerConstraintMargin > 0 || volumeConstraintMargin > 0) ? "WARN" : "OK";
  const char* sensitivityHealth =
    (fsensNonFiniteCountLocal > 0 || gsensPowerNonFiniteCountLocal > 0 || gsensVolNonFiniteCountLocal > 0)
    ? "WARN"
    : "OK";
  const char* solverHealth =
    solverSeverity >= 2 ? "WARN" : (solverSeverity == 1 ? "OBS" : "OK");

  std::ofstream debugLogFile(debugOptimizerLogFile, std::ios::app);
  const int itemWidth = 18;
  const int valueWidth = 24;
  const int tableWidth = itemWidth + 4*valueWidth;
  auto fmtScalar =
  [&](const scalar value)
  {
    std::ostringstream os;
    os << std::scientific << std::setprecision(6) << value;
    return os.str();
  };
  auto fmtCount =
  [&](const label value)
  {
    std::ostringstream os;
    os << value;
    return os.str();
  };
  auto writeSectionHeader =
  [&](const char* title)
  {
    debugLogFile
      << "\n[" << title << "]\n"
      << "  "
      << std::left
      << std::setw(itemWidth) << "Item"
      << std::setw(valueWidth) << "Value 1"
      << std::setw(valueWidth) << "Value 2"
      << std::setw(valueWidth) << "Value 3"
      << std::setw(valueWidth) << "Value 4"
      << "\n"
      << "  "
      << std::string(tableWidth, '-')
      << "\n";
  };
  auto writeRow =
  [&](const char* item, const std::string& value1, const std::string& value2 = "", const std::string& value3 = "", const std::string& value4 = "")
  {
    debugLogFile
      << "  "
      << std::left
      << std::setw(itemWidth) << item
      << std::setw(valueWidth) << value1
      << std::setw(valueWidth) << value2
      << std::setw(valueWidth) << value3
      << std::setw(valueWidth) << value4
      << "\n";
  };

  debugLogFile
    << "\n"
    << "====================================================================================================\n"
    << " Optimizer Debug Dashboard\n"
    << "----------------------------------------------------------------------------------------------------\n";

  writeSectionHeader("Overview");
  writeRow("Run", "Iter="+fmtCount(opt), "Clock[s]="+fmtScalar(clockTime), "ObjMode="+std::string(objectiveLabel.c_str()), "ObjVal="+fmtScalar(objectiveValue));
  writeRow("Health", std::string("Design=")+designHealth, std::string("Constr=")+constraintHealth, std::string("Sens=")+sensitivityHealth, std::string("Solver=")+solverHealth);
  writeRow("Convergence", "objRel="+fmtScalar(objectiveRelativeDelta), "maxViol="+fmtScalar(maxConstraintViolation), "xhStep="+fmtScalar(maxAbsXhChangeLocal), "streak="+fmtCount(optimizerConvergedStreak)+"/"+fmtCount(optStallWindow));
  writeRow("Stop", "candidate="+std::string(convergenceCandidate ? "yes" : "no"), "stop="+std::string(stopOptimization ? "yes" : "no"), "reason="+std::string(optimizerStopReason.c_str()));

  writeSectionHeader("Design Diagnostics");
  writeRow("x bounds", "min="+fmtScalar(xLocalMin), "max="+fmtScalar(xLocalMax), "mean="+fmtScalar(xMean));
  writeRow("xp bounds", "min="+fmtScalar(xpLocalMin), "max="+fmtScalar(xpLocalMax));
  writeRow("xh bounds", "min="+fmtScalar(xhLocalMin), "max="+fmtScalar(xhLocalMax), "mean="+fmtScalar(xhMean));
  writeRow("xh phases", "gray="+fmtScalar(grayVolumeFraction), "solid="+fmtScalar(solidVolumeFraction), "fluid="+fmtScalar(fluidVolumeFraction));
  writeRow("x step", "max="+fmtScalar(maxAbsXChangeLocal), "l2="+fmtScalar(l2XChange), "l2/sqrtN="+fmtScalar(normalizedL2XChange));
  writeRow("xh step", "max="+fmtScalar(maxAbsXhChangeLocal), "l2="+fmtScalar(l2XhChange), "l2/sqrtN="+fmtScalar(normalizedL2XhChange));

  writeSectionHeader("Filter Diagnostics");
  writeRow("projection", "del="+fmtScalar(del), "eta5="+fmtScalar(eta5));
  writeRow("drho bounds", "min="+fmtScalar(drhoLocalMin), "max="+fmtScalar(drhoLocalMax));
  writeRow("obj sens", "min="+fmtScalar(fsensLocalMin), "max="+fmtScalar(fsensLocalMax));
  writeRow("power sens", "min="+fmtScalar(gsensPowerLocalMin), "max="+fmtScalar(gsensPowerLocalMax));
  writeRow("vol sens", "min="+fmtScalar(gsensVolLocalMin), "max="+fmtScalar(gsensVolLocalMax));

  writeSectionHeader("Objective And Constraints");
  writeRow("objective", "value="+fmtScalar(objectiveValue), "delta="+fmtScalar(objectiveDelta));
  writeRow("power diss", "actual="+fmtScalar(PowerDissActual), "act/base="+fmtScalar(PowerDissActualOverBase), "constraint="+fmtScalar(PowerDiss), "margin="+fmtScalar(powerConstraintMargin));
  writeRow("volume", "value="+fmtScalar(V), "margin="+fmtScalar(volumeConstraintMargin));
  writeRow("thermal", "meanT="+fmtScalar(MeanT), "maxT="+fmtScalar(maxT));

  writeSectionHeader("MMA Diagnostics");
  writeRow("controls", "moveLimit="+fmtScalar(movlim), "moveHits="+fmtCount(xMoveLimitHitCountLocal));
  writeRow("bound clamps", "lower="+fmtCount(xLowerBoundCountLocal), "upper="+fmtCount(xUpperBoundCountLocal), "nonFiniteX="+fmtCount(xNonFiniteCountLocal));

  writeSectionHeader("Sensitivity Diagnostics");
  writeRow("objective", "min="+fmtScalar(fsensLocalMin), "max="+fmtScalar(fsensLocalMax), "l2="+fmtScalar(fsensL2), "nonFinite="+fmtCount(fsensNonFiniteCountLocal));
  writeRow("power", "min="+fmtScalar(gsensPowerLocalMin), "max="+fmtScalar(gsensPowerLocalMax), "l2="+fmtScalar(gsensPowerL2), "nonFinite="+fmtCount(gsensPowerNonFiniteCountLocal));
  writeRow("volume", "min="+fmtScalar(gsensVolLocalMin), "max="+fmtScalar(gsensVolLocalMax), "l2="+fmtScalar(gsensVolL2), "nonFinite="+fmtCount(gsensVolNonFiniteCountLocal));

  writeSectionHeader("Solver Diagnostics");
  writeRow("config", "primal=60", "adjPower=40", "adjHeat=40");
  writeRow("primal flow", "|U|max="+fmtScalar(UMaxMagLocal), "|U|mean="+fmtScalar(UMeanMag), "pMin="+fmtScalar(pLocalMin), "pMax="+fmtScalar(pLocalMax));
  writeRow("thermal field", "TMin="+fmtScalar(TLocalMin), "TMax="+fmtScalar(TLocalMax));
  writeRow("adjoint power", "|Ua|max="+fmtScalar(UaMaxMagLocal), "|Ua|mean="+fmtScalar(UaMeanMag), "paMin="+fmtScalar(paLocalMin), "paMax="+fmtScalar(paLocalMax));
  writeRow("adjoint heat", "|Ub|max="+fmtScalar(UbMaxMagLocal), "|Ub|mean="+fmtScalar(UbMeanMag), "pbMin="+fmtScalar(pbLocalMin), "pbMax="+fmtScalar(pbLocalMax));

  writeSectionHeader("Interpolation Diagnostics");
  writeRow("interp bounds", "alphaMin="+fmtScalar(alphaLocalMin), "alphaMax="+fmtScalar(alphaLocalMax), "DTMin="+fmtScalar(DTLocalMin), "DTMax="+fmtScalar(DTLocalMax));
  writeRow("continuation", "alphaMax="+fmtScalar(alphaMax.value()), "qu="+fmtScalar(qu));
  writeRow("thermal props", "heatSource="+fmtScalar(heatSource.value()), "thermalSrc="+fmtScalar(thermalSource.value()), "rhoCp="+fmtScalar(rhoCp.value()));

  std::ofstream debugJsonFile(debugOptimizerJsonFile, std::ios::app);
  debugJsonFile
    << std::scientific << std::setprecision(8)
    << "{"
    << "\"iteration\":" << opt
    << ",\"clockTime\":" << clockTime
    << ",\"health\":{"
      << "\"design\":\"" << designHealth << "\""
      << ",\"constraints\":\"" << constraintHealth << "\""
      << ",\"sensitivity\":\"" << sensitivityHealth << "\""
      << ",\"solver\":\"" << solverHealth << "\""
    << "}"
    << ",\"convergence\":{"
      << "\"objectiveRelativeDelta\":" << objectiveRelativeDelta
      << ",\"maxConstraintViolation\":" << maxConstraintViolation
      << ",\"xStepMax\":" << maxAbsXChangeLocal
      << ",\"xhStepMax\":" << maxAbsXhChangeLocal
      << ",\"xStepL2OverSqrtN\":" << normalizedL2XChange
      << ",\"xhStepL2OverSqrtN\":" << normalizedL2XhChange
      << ",\"streak\":" << optimizerConvergedStreak
      << ",\"stallWindow\":" << optStallWindow
      << ",\"candidate\":" << (convergenceCandidate ? 1 : 0)
      << ",\"stop\":" << (stopOptimization ? 1 : 0)
      << ",\"reason\":\"" << optimizerStopReason << "\""
    << "}"
    << ",\"design\":{"
      << "\"xMin\":" << xLocalMin
      << ",\"xMax\":" << xLocalMax
      << ",\"xMean\":" << xMean
      << ",\"xpMin\":" << xpLocalMin
      << ",\"xpMax\":" << xpLocalMax
      << ",\"xhMin\":" << xhLocalMin
      << ",\"xhMax\":" << xhLocalMax
      << ",\"xhMean\":" << xhMean
      << ",\"xhGrayVolumeFraction\":" << grayVolumeFraction
      << ",\"xhSolidVolumeFraction\":" << solidVolumeFraction
      << ",\"xhFluidVolumeFraction\":" << fluidVolumeFraction
      << ",\"xStepMax\":" << maxAbsXChangeLocal
      << ",\"xStepL2\":" << l2XChange
      << ",\"xhStepMax\":" << maxAbsXhChangeLocal
      << ",\"xhStepL2\":" << l2XhChange
    << "}"
    << ",\"filter\":{"
      << "\"del\":" << del
      << ",\"eta5\":" << eta5
      << ",\"drhoMin\":" << drhoLocalMin
      << ",\"drhoMax\":" << drhoLocalMax
      << ",\"objectiveSensitivityFilteredMin\":" << fsensLocalMin
      << ",\"objectiveSensitivityFilteredMax\":" << fsensLocalMax
      << ",\"powerConstraintSensitivityFilteredMin\":" << gsensPowerLocalMin
      << ",\"powerConstraintSensitivityFilteredMax\":" << gsensPowerLocalMax
      << ",\"volumeConstraintSensitivityFilteredMin\":" << gsensVolLocalMin
      << ",\"volumeConstraintSensitivityFilteredMax\":" << gsensVolLocalMax
    << "}"
    << ",\"objective\":{"
      << "\"name\":\"" << objectiveLabel << "\""
      << ",\"value\":" << objectiveValue
      << ",\"delta\":" << objectiveDelta
      << ",\"powerDissipationActual\":" << PowerDissActual
      << ",\"powerDissipationActualOverPowerDiss0\":" << PowerDissActualOverBase
      << ",\"powerDissipationConstraintValue\":" << PowerDiss
      << ",\"powerDissipationConstraintMargin\":" << powerConstraintMargin
      << ",\"volumeConstraintValue\":" << V
      << ",\"volumeConstraintMargin\":" << volumeConstraintMargin
      << ",\"meanT\":" << MeanT
      << ",\"maxT\":" << maxT
    << "}"
    << ",\"mma\":{"
      << "\"iteration\":" << opt
      << ",\"moveLimit\":" << movlim
      << ",\"moveLimitHitCount\":" << xMoveLimitHitCountLocal
      << ",\"lowerBoundCount\":" << xLowerBoundCountLocal
      << ",\"upperBoundCount\":" << xUpperBoundCountLocal
      << ",\"designNonFiniteCount\":" << xNonFiniteCountLocal
    << "}"
    << ",\"sensitivity\":{"
      << "\"objectiveMin\":" << fsensLocalMin
      << ",\"objectiveMax\":" << fsensLocalMax
      << ",\"objectiveL2\":" << fsensL2
      << ",\"objectiveNonFiniteCount\":" << fsensNonFiniteCountLocal
      << ",\"powerConstraintMin\":" << gsensPowerLocalMin
      << ",\"powerConstraintMax\":" << gsensPowerLocalMax
      << ",\"powerConstraintL2\":" << gsensPowerL2
      << ",\"powerConstraintNonFiniteCount\":" << gsensPowerNonFiniteCountLocal
      << ",\"volumeConstraintMin\":" << gsensVolLocalMin
      << ",\"volumeConstraintMax\":" << gsensVolLocalMax
      << ",\"volumeConstraintL2\":" << gsensVolL2
      << ",\"volumeConstraintNonFiniteCount\":" << gsensVolNonFiniteCountLocal
    << "}"
    << ",\"solver\":{"
      << "\"primalFlowSweeps\":" << 60
      << ",\"adjointFlowPowerSweeps\":" << 40
      << ",\"adjointFlowHeatSweeps\":" << 40
      << ",\"uMaxMag\":" << UMaxMagLocal
      << ",\"uMeanMag\":" << UMeanMag
      << ",\"pMin\":" << pLocalMin
      << ",\"pMax\":" << pLocalMax
      << ",\"tMin\":" << TLocalMin
      << ",\"tMax\":" << TLocalMax
      << ",\"uaMaxMag\":" << UaMaxMagLocal
      << ",\"uaMeanMag\":" << UaMeanMag
      << ",\"paMin\":" << paLocalMin
      << ",\"paMax\":" << paLocalMax
      << ",\"ubMaxMag\":" << UbMaxMagLocal
      << ",\"ubMeanMag\":" << UbMeanMag
      << ",\"pbMin\":" << pbLocalMin
      << ",\"pbMax\":" << pbLocalMax
    << "}"
    << ",\"interpolation\":{"
      << "\"alphaMin\":" << alphaLocalMin
      << ",\"alphaMax\":" << alphaLocalMax
      << ",\"dtMin\":" << DTLocalMin
      << ",\"dtMax\":" << DTLocalMax
      << ",\"alphaMaxSchedule\":" << alphaMax.value()
      << ",\"qu\":" << qu
      << ",\"heatSource\":" << heatSource.value()
      << ",\"thermalSource\":" << thermalSource.value()
      << ",\"rhoCp\":" << rhoCp.value()
    << "}"
    << "}\n";
}

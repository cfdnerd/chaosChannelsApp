//Author: Yu Minghao    Updated: May 2020 
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "simpleControl.H"
#include "fvOptions.H"//
#include "MMA/MMA.h"
#include <cstring>
#include <iomanip>
#include <sstream>
#include <diff.c>

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFvOptions.H"//
    #include "createFields.H"
    #include "readTransportProperties.H" 
    #include "initContinuityErrs.H"
    #include "readThermalProperties.H" 
    #include "opt_initialization.H"
    while (simple.loop(runTime))
    {
        #include "solverConvergenceReset.H"
        #include "update.H"
        // Respect case-defined inlet BCs; derive turbulence inlet state from current primal inlet velocity.
        scalar currentInletSpeedForLog(inletSpeed);
        if (updateTurbulenceInletFromVelocity)
        {
            scalar currentInletSpeedSumLocal(0.0);
            label currentInletFaceCountLocal(0);
            const auto& inletU = U.boundaryField()[inletPatchID];
            if (inletU.size() > 0)
            {
                forAll(inletU, faceI)
                {
                    currentInletSpeedSumLocal += mag(inletU[faceI]);
                }
                currentInletFaceCountLocal = inletU.size();
            }
            scalar currentInletSpeedSumGlobal(currentInletSpeedSumLocal);
            label currentInletFaceCountGlobal(currentInletFaceCountLocal);
            reduce(currentInletSpeedSumGlobal, sumOp<scalar>());
            reduce(currentInletFaceCountGlobal, sumOp<label>());
            scalar currentInletSpeed(inletSpeed);
            if (currentInletFaceCountGlobal > 0)
            {
                currentInletSpeed =
                    currentInletSpeedSumGlobal/scalar(currentInletFaceCountGlobal);
            }
            currentInletSpeed = Foam::max(currentInletSpeed, scalar(SMALL));
            currentInletSpeedForLog = currentInletSpeed;
            const scalar turbulenceK =
                Foam::max(1.5*Foam::sqr(turbulenceIntensity*currentInletSpeed), scalar(SMALL));
            const scalar turbulenceEpsilon =
                Foam::max
                (
                    Foam::pow(turbulenceModelCmu, 0.75)
                    *Foam::pow(turbulenceK, 1.5)
                    /Foam::max(turbulenceLengthScale, scalar(SMALL)),
                    scalar(SMALL)
                );
            if
            (
                mesh.foundObject<volScalarField>("k")
                &&
                mesh.foundObject<volScalarField>("epsilon")
            )
            {
                volScalarField& kField =
                    const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("k"));
                volScalarField& epsilonField =
                    const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("epsilon"));
                kField.boundaryFieldRef()[inletPatchID] = turbulenceK;
                epsilonField.boundaryFieldRef()[inletPatchID] = turbulenceEpsilon;
            }
        }

        U.correctBoundaryConditions();

        Info<< "Running single-case optimization with inlet speed "
            << currentInletSpeedForLog << " m/s" << endl;
        #include "Primal_U.H"
        if (stopOptimization)
        {
            break;
        }
        #include "Primal_T.H"
        if (stopOptimization)
        {
            break;
        }
        #include "AdjointHeat_Tb.H"
        #include "AdjointHeat_Ub.H"
        #include "AdjointFlow_Ua.H"
        #include "costfunction.H"
        if (stopOptimization)
        {
            break;
        }
        #include "sensitivity.H"
        runTime.write();
        if (stopOptimization)
        {
            Info<< "Optimizer stop criteria met at iteration " << opt
                << " (" << optimizerStopReason << ")." << endl;
            break;
        }
    }
    #include "finalize.H"
    return 0;
}

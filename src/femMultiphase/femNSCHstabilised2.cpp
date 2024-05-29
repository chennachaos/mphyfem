
#include "femNSCHstabilised.h"
#include "util.h"
#include "ElementBase.h"
#include <chrono>
#include "KimMoinFlow.h"
#include "metis.h"
#include "QuadratureUtil.h"
#include "BasisFunctionsLagrange.h"
#include "stabilisationRoutines.h"
#include "MyTime.h"
#include "TimeFunction.h"



extern  std::vector<unique_ptr<TimeFunction> > timeFunction;
extern  MyTime                 myTime;



int femNSCHstabilised::setSolver(int slv, int *parm, bool cIO)
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n     femNSCHstabilised::setSolver()  .... STARTED ...\n\n");


    NavierStokes.setSolver(slv, parm, cIO);
    CahnHilliard.setSolver(slv, parm, cIO);

    PetscPrintf(MPI_COMM_WORLD, "\n\n     femNSCHstabilised::setSolver()  .... ENDED ...\n\n");

    return 0;
}







int  femNSCHstabilised::solveFullyImplicit()
{
    PetscPrintf(MPI_COMM_WORLD, " Solving with the Fully-Implicit Scheme \n\n\n");

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////


    int  stepsCompleted=1, ii, jj;

    double  phi, error, timerVal;
    double  mudiff = 0.5*(fluidDynamicViscosity[0]-fluidDynamicViscosity[1]), muavg = 0.5*(fluidDynamicViscosity[0]+fluidDynamicViscosity[1]);
    double  rhodiff = 0.5*(fluidDensity[0]-fluidDensity[1]), rhoavg = 0.5*(fluidDensity[0]+fluidDensity[1]);

    VectorXd  solnNSpreviter=NavierStokes.soln, solnCHpreviter=CahnHilliard.soln, solnNSerror, solnCHerror;
    
    NavierStokes.muNodal.resize(mesh->getnNodeGlobal());
    NavierStokes.muNodal.setZero();
    
    NavierStokes.rhoNodal = NavierStokes.muNodal;
    NavierStokes.phiNodal = NavierStokes.muNodal;
    NavierStokes.etaNodal = NavierStokes.muNodal;

    CahnHilliard.veloConvection.resize(mesh->getnNodeGlobal() * mesh->getNDIM());
    


    setInitialConditions();
    postProcess();
    writeOutputData();

    convergedFlagPrev = convergedFlag = false;

    computerTimeTimeLoop = MPI_Wtime();

    //Time loop
    while( (myTime.cur <= (timeFinal-EPSILON)) && (stepsCompleted <= stepsMax) )
    {
        // do a time update: reset variables and flags, and prediction step of fields
        timeUpdate();


        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n");
        PetscPrintf(MPI_COMM_WORLD, " Time step number     =  %d  \n", stepsCompleted);
        PetscPrintf(MPI_COMM_WORLD, " Time step size       =  %f  \n", myTime.dt);
        PetscPrintf(MPI_COMM_WORLD, " Current time         =  %f  \n", myTime.cur);
        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n");


        convergedFlagPrev = convergedFlag;
        convergedFlag = false;

        NSCONVERGED = false;
        CHCONVERGED = false;


        MPI_Barrier(MPI_COMM_WORLD);

        PetscPrintf(MPI_COMM_WORLD, " Outer Iteration loop \n");
        for(int outeriter=1; outeriter<=iterationsMax; outeriter++)
        {
            firstIteration = (outeriter == 1);

            // Compute the velocity and acceleration and the respective values at n+af, n+am
            NavierStokes.updateIterStep();
            CahnHilliard.updateIterStep();


            timerVal = MPI_Wtime();

            CahnHilliard.getPhi(NavierStokes.phiNodal);
            CahnHilliard.getEta(NavierStokes.etaNodal);

            for(ii=0; ii<mesh->getnNodeGlobal(); ii++)
            {
              phi = NavierStokes.phiNodal[ii];
              if(abs(phi) > 1.0)
              {
                  phi = (phi > 0) ? 1.0 : -1.0;
              }
              NavierStokes.phiNodal[ii] = phi;

              NavierStokes.muNodal[ii]  = mudiff*phi  + muavg;
              NavierStokes.rhoNodal[ii] = rhodiff*phi + rhoavg;
            }

            // solve Navier-Stokes
            PetscPrintf(MPI_COMM_WORLD, " Solving Navier-Stokes ... \n\n\n");
            NavierStokes.solveTimeStep();

            NSCONVERGED = NavierStokes.converged();

            if(!NSCONVERGED)
            {
                convergedFlag = false;
                PetscPrintf(MPI_COMM_WORLD, "\n Navier-Stokes solver didnot converge... \n");
                PetscPrintf(MPI_COMM_WORLD, "\n Exiting outer loop... \n");
                break;
            }

            NavierStokes.getVelocity(CahnHilliard.veloConvection);

            // solve Chan-Hilliard
            PetscPrintf(MPI_COMM_WORLD, " Solving Cahn-Hilliard ... \n\n\n");
            CahnHilliard.solveTimeStep();

            CHCONVERGED = CahnHilliard.converged();

            cout << " CHCONVERGED = " << CHCONVERGED << endl;
            
            if(!CHCONVERGED)
            {
                convergedFlag = false;
                PetscPrintf(MPI_COMM_WORLD, "\n Cahn-Hilliard solver didnot converge... \n");
                PetscPrintf(MPI_COMM_WORLD, "\n Exiting outer loop... \n");
                break;
            }


            timerVal = MPI_Wtime() - timerVal;
            computerTimeAssembly += timerVal;
            //PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time for matrix assembly = %f seconds \n\n", timerVal );

            solnNSerror = NavierStokes.soln - solnNSpreviter;
            solnCHerror = CahnHilliard.soln - solnCHpreviter;

            error = sqrt(solnNSerror.squaredNorm() + solnCHerror.squaredNorm());

            PetscPrintf(MPI_COMM_WORLD, "\n Outer Iteration : %2d ... %12.6E \n", outeriter, error);

            if (error < 1.0e-5)
            {
                PetscPrintf(MPI_COMM_WORLD, "\n Outer Iteration loop converged... \n");

                convergedFlag = true;

                break;
            }

            solnNSpreviter = NavierStokes.soln;
            solnCHpreviter = CahnHilliard.soln;

            MPI_Barrier(MPI_COMM_WORLD);

        } //Iteration Loop


        // if the residual is converged, then save the DOFs vectors
        if( convergedFlag )
        {
            PetscPrintf(MPI_COMM_WORLD, " Postprocessing... \n\n");

            timerVal = MPI_Wtime();

            if( stepsCompleted%outputFreq == 0 )
              postProcess();

            writeOutputData();

            computerTimePostprocess += (MPI_Wtime() - timerVal);

            //saveSolution();
            NavierStokes.saveSolution();
            CahnHilliard.saveSolution();

            myTime.stck();

            stepsCompleted++;
        }
        else
        {
            myTime.cut();

            //reset();
            NavierStokes.reset();
            CahnHilliard.reset();
        }
 
    } //Time loop

    computerTimeTimeLoop = MPI_Wtime() - computerTimeTimeLoop;

    PetscPrintf(MPI_COMM_WORLD, "\n\n Simulation reached the specified final time or maximum steps specified ... \n\n\n");

    return 0;
}



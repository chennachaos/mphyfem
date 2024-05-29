/*
! Program for solid mechanics using FEM
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 30-March-2023
! Place : Edinburgh, UK
!
!
*/


#include "headersVTK.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "femCahnHilliard.h"
#include <unistd.h>
#include "Global.h"


using namespace std;



int main(int argc, char* argv[])
{
    //Set the input file name
    string  meshfile      = "./inputs/chmesh.msh";
    string  configfile    = "./inputs/config";
    string  petscfile     = "./inputs/petsc_options.dat";

    PetscInitialize(NULL, NULL, petscfile.c_str(), NULL);
    //PetscInitialize(NULL, NULL, NULL, NULL);
    PetscMemorySetGetMaximumUsage();


    PetscPrintf(MPI_COMM_WORLD, "\n    Working directory       = %s \n ", get_current_dir_name());
    PetscPrintf(MPI_COMM_WORLD, "\n    Input mesh file         = %s \n ", meshfile.c_str());
    PetscPrintf(MPI_COMM_WORLD, "\n    Configuration file      = %s \n ", configfile.c_str());
    PetscPrintf(MPI_COMM_WORLD, "\n    PETSc options file      = %s \n ", petscfile.c_str());


    femCahnHilliard  chfem;

    chfem.readInputGMSH(meshfile);
    
    PetscPrintf(MPI_COMM_WORLD, "\n    Reading configuration file  \n ");

    chfem.readConfiguration(configfile);

    chfem.prepareInputData();

    //chfem.setInitialConditions();
    //chfem.postProcess();


    double timerStart, timerEnd;

    timerStart = MPI_Wtime();

    chfem.setSolver(1, NULL, false);

    timerEnd = MPI_Wtime();

    PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time = %f seconds \n\n", timerEnd - timerStart );

    chfem.solveFullyImplicit();
    MPI_Barrier(MPI_COMM_WORLD);

    chfem.printComputerTimes();
    MPI_Barrier(MPI_COMM_WORLD);

    PetscPrintf(MPI_COMM_WORLD, "\n\n\n Program is successful \n\n\n ");

    // this needs to done explicitly. Otherwise all the allocated PETSc objects will not be destroyed.
    chfem.deallocatePetscObjects();
    MPI_Barrier(MPI_COMM_WORLD);

    PetscFinalize();

    return 0;
}

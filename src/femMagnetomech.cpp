/*
! Program for incompressible Navier-Stokes using Stabilised Finite Element Method
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 06-May-2020
! Place : Swansea, UK
!
!
*/


#include "headersVTK.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "femMagnetomech.h"
#include "Global.h"
#include "FunctionsProgram.h"


using namespace std;




int main(int argc, char* argv[])
{
    //Set default config file name
    string  configfname = "./inputs/config";

    //Update the config file name if specified from the command line
    if(argc == 2)
    {
      configfname  = "./inputs/" + string(argv[1]);
    }

    printf("\n    Config file name ... = %s \n ", configfname.c_str());

    // PETSc options file
    string petscoptionsfile = "./inputs/petsc_options.dat";

    ifstream fin(petscoptionsfile);

    if(fin.fail())
    {
      cout << "\n\n ERROR: Could not open the petsc_options.dat file \n\n" << endl;
      exit(1);
    }
    fin.close();

    PetscInitialize(NULL, NULL, "./inputs/petsc_options.dat", NULL);


    femMagnetomech  magnfem;

    magnfem.readInput(configfname);

    magnfem.prepareInputData();

    double timerStart, timerEnd;

    timerStart = MPI_Wtime();

    magnfem.setSolver(2);

    magnfem.solveFullyImplicit();

    timerEnd = MPI_Wtime();

    printf("\n\n Elapsed time = %f seconds \n\n", timerEnd - timerStart );

    printf("\n\n\n Program is successful \n\n\n ");

    PetscFinalize(); //CHKERRQ(ierr);

    return 0;
}




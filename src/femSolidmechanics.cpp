/*
! Program for solid mechanics Finite Element Method
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 19-May-2024
! Place : Edinburgh, UK
!
!
*/


#include "headersVTK.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "femSolidmechanics.h"
#include "Global.h"


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



    femSolidmechanics  solidfem;

    solidfem.readInput(configfname);

    solidfem.prepareInputData();

    double timerStart, timerEnd;

    timerStart = MPI_Wtime();

    solidfem.setSolver(2);

    solidfem.solveFullyImplicit();

    timerEnd = MPI_Wtime();

    printf("\n\n Elapsed time = %f seconds \n\n", timerEnd - timerStart );

    //solidfem.computeElementErrors(0);
    //solidfem.printComputerTimes();

    printf("\n\n\n Program is successful \n\n\n ");

    return 0;
}


//

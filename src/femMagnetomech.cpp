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
    /*
    //Set the input file name
    //The file name is specified from the command line
    if(argc == 0)
    {
        cerr << " Error in input data " << endl;
        cerr << " You must enter name of input file" << endl;
        cerr << " Aborting..." << endl;
    }

    string  inputfile  = argv[1];

    printf("\n    Input file         = %s \n ", inputfile.c_str());
    */

    string  inputfile;

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

    magnfem.readInput(inputfile);

    magnfem.prepareInputData();

    double timerStart, timerEnd;

    timerStart = MPI_Wtime();
    magnfem.setSolver(2);
    timerEnd = MPI_Wtime();

    magnfem.solveFullyImplicit();

    //magnfem.postProcess();
    //magnfem.printComputerTimes();

    printf("\n\n\n Program is successful \n\n\n ");

    PetscFinalize(); //CHKERRQ(ierr);

    return 0;
}




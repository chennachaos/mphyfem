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
#include "GrowthFEM.h"
#include "Global.h"
#include "FunctionsProgram.h"


using namespace std;

extern Files files;




int main(int argc, char* argv[])
{
    PetscInitialize(NULL, NULL, "petsc_options.dat", NULL);

    //Set the input file name
    //The file name is specified from the command line
    if(argc == 0)
    {
        cerr << " Error in input data " << endl;
        cerr << " You must enter name of input file" << endl;
        cerr << " Aborting..." << endl;
    }

    cout << argv[1] << endl;
    cout << argv[2] << endl;

    files.projDir = argv[1];
    files.Ifile   = argv[2];

    files.Ofile.free().append(files.Ifile)[0] = 'O';
    files.Tfile.free().append(files.Ifile)[0] = 'T';
    files.Pfile.free().append(files.Ifile)[0] = 'P';


    //printf("\n    Input mesh file         = %s \n ", meshfile.c_str());
    //printf("\n    Control parameters file = %s \n ", controlfile.c_str());

    PetscPrintf(MPI_COMM_WORLD, "       project directory   : %s   \n", files.projDir.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       input file name     : %s   \n", files.Ifile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       output file name    : %s   \n", files.Ofile.asCharArray());
    PetscPrintf(MPI_COMM_WORLD, "       time plot file name : %s   \n", files.Tfile.asCharArray());


    if( prgFileExist(files.projDir,files.Ifile) )
    {
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  Directory and file exist \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
    }
    else
    {
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  No such project directory or input file! \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  Program has been aborted ... \n\n");
      PetscPrintf(MPI_COMM_WORLD, "  ================================== \n\n");

      PetscFinalize(); //CHKERRQ(ierr);

      return 1;
    }


    //MyString pathAndFile;


    GrowthFEM  growthfem;

    //string  meshfile    = argv[1];
    //string  controlfile = argv[2];


    string  inpfilename    = string(files.projDir.asCharArray()) + "/" + string(files.Ifile.asCharArray());

    growthfem.readInput(inpfilename);

    //growthfem.readControlParameters(controlfile);

    growthfem.prepareInputData();

    double timerStart, timerEnd;

    //timerStart = MPI_Wtime();

    growthfem.setSolver();
    //timerEnd = MPI_Wtime();

    printf("\n\n Elapsed time = %f seconds \n\n", timerEnd - timerStart );

    growthfem.plotGeom();

    //growthfem.postProcess();

    growthfem.solve();

    //growthfem.printComputerTimes();

    files.reset();

    printf("\n\n\n Program is successful \n\n\n ");

    PetscFinalize(); //CHKERRQ(ierr);

    return 0;
}




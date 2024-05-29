/*
! Program for 2D explicit finite element analysis of incompressible Navier-Stokes
!
!
! Author: Dr. Chennakesava Kadapa
! Date  : 17-May-2018
! Place : Swansea, UK
!
!
*/


#include "headersVTK.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "femINSmixed.h"
#include <unistd.h>
#include "Global.h"


using namespace std;



int main(int argc, char* argv[])
{
    string  meshfile       = "./inputs/fluidmesh.msh";
    string  configfile     = "./inputs/config";

    cout << " Working dir     : " << get_current_dir_name() << endl;
    cout << " meshfile        : " << meshfile << endl;
    cout << " controlfile     : " << configfile << endl;

    femINSmixed  cfdfem;

    cfdfem.readInputGMSH(meshfile);

    cfdfem.readConfiguration(configfile);

    cfdfem.prepareInputData();

    cfdfem.setSolver(0);

    //string  restfile = "Result.rst";
    //cfdfem.writeReadResult(0, restfile, 0);

    //cfdfem.solveFullyImplicit();

    cout << " Program is successful \n " << endl;

    return 1;
}

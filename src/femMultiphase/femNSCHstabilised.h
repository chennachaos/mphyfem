
#ifndef incl_femNSCHstabilised_CLASS_h
#define incl_femNSCHstabilised_CLASS_h
#define EIGEN_SUPERLU_SUPPORT

#include <string.h>
#include <vector>
#include <fstream>
#include <memory>
#include "femBase.h"
#include "SolverPetsc.h"
#include "BoundaryPatch.h"
#include "util.h"
#include "myMathFunction.h"
#include "myMesh.h"
#include "femINSstabilised.h"
#include "femCahnHilliard.h"


using namespace std;




/*
class NavierStokes
{
    public:

        double    conv_tol, spectralRadius, rhsNorm, rhsNormPrev;
        double    bodyForce[3];

        string    infilename, dirname, timeIntegrationScheme;

        bool      firstIteration, convergedFlag, convergedFlagPrev;

        PetscInt  AlgoType, SCHEME_TYPE;
        PetscErrorCode  errpetsc;
        PetscInt  n_mpi_procs, this_mpi_proc;

        PetscInt  ndof, ntotdofs_local, ntotdofs_global, bodyForceTimeFunction;

        vector<int>              assyForSoln, dofs_specified;
        vector<int>              elem_proc_id, node_proc_id, node_map_get_old, node_map_get_new;
        vector<vector<int> >     forAssyVecAll, globalDOFnumsAll;

        VectorXd  td, soln, solnPrev, solnPrev2, solnPrev3, solnPrev4, solnCur, solnInit, solnApplied;
        VectorXd  solnDot, solnDotPrev, solnDotCur, solnExtrap;

        unique_ptr<SolverPetsc>  solverPetsc;

        vector<myMathFunction>   initialConitions;

        vector<unique_ptr<BoundaryCondition> >   boundaryConitions;
};




class CahnHilliard
{
    public:

        double    conv_tol, spectralRadius, rhsNorm, rhsNormPrev;
        double    bodyForce[3];

        string    infilename, dirname, timeIntegrationScheme;

        bool      firstIteration, convergedFlag, convergedFlagPrev;

        PetscInt  AlgoType, SCHEME_TYPE;
        PetscErrorCode  errpetsc;
        PetscInt  n_mpi_procs, this_mpi_proc;

        PetscInt  ndof, ntotdofs_local, ntotdofs_global, bodyForceTimeFunction;

        vector<int>              assyForSoln, dofs_specified;
        vector<int>              elem_proc_id, node_proc_id, node_map_get_old, node_map_get_new;
        vector<vector<int> >     forAssyVecAll, globalDOFnumsAll;

        VectorXd  td, soln, solnPrev, solnPrev2, solnPrev3, solnPrev4, solnCur, solnInit, solnApplied;
        VectorXd  solnDot, solnDotPrev, solnDotCur, solnExtrap;

        unique_ptr<SolverPetsc>  solverPetsc;

        vector<myMathFunction>   initialConitions;

        vector<unique_ptr<BoundaryCondition> >   boundaryConitions;

};
*/



class femNSCHstabilised
{
    //private:

    public:

        //myMesh  mesh;
        shared_ptr<myMesh>  mesh;

        PetscErrorCode  errpetsc;
        PetscInt  n_mpi_procs, this_mpi_proc;

        PetscInt  ndim;
        PetscInt  stepsMax, iterationsMax, outputFreq, fileCount;

        double    fluidDensity[2], fluidDynamicViscosity[2];

        VectorXd  ForceVectorExternal;

        bool      firstIteration, convergedFlag, convergedFlagPrev;
        bool      NSCONVERGED, CHCONVERGED;

        double    computerTimeAssembly, computerTimeSolver, computerTimePostprocess, computerTimePattern, computerTimeTimeLoop;
        double    conv_tol, spectralRadius, timeFinal, rhsNorm, rhsNormPrev;

        string    infilename, dirname, timeIntegrationScheme;

        femINSstabilised  NavierStokes;

        femCahnHilliard   CahnHilliard;

    public:

        femNSCHstabilised();

        ~femNSCHstabilised();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  readInputGMSH(string& fname);

        int  readConfiguration(string& fname);

        int  readConfiguration(string& fname, string& fnamefluid, string& fnamephase);

        int  readInitialConditions(ifstream& infile, string& line);

        int  readBodyForce(ifstream& infile, string& line);

        int  readBoundariesData(ifstream& infile, string& line);

        int  readTimeFunctions(ifstream& infile, string& line);

        int  readSolverDetails(ifstream& infile, string& line);

        int  deallocatePetscObjects();

        int  setBoundaryConditions();

        int  prepareInputData();

        int  readInputData(string& fname);

        int  readFluidProperties(ifstream& infile, string& line);

        int  readResult(string&);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  setSolver(int, int *parm = NULL, bool cIO = false);

        int  prepareMatrixPattern();

        int  setTimeParam();

        int  timeUpdate();

        int  updateIterStep();

        int  reset();

        int  saveSolution();

        int  addExternalForces(double loadFact);

        int  computeElementErrors(int);

        int  setInitialConditions();

        int  solveFullyImplicit();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  writeOutputData();

        int  writeResult(string&);

        int  printComputerTimes();

        int  postProcess();
};






#endif







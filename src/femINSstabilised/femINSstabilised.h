
#ifndef incl_femINSstabilised_CLASS_h
#define incl_femINSstabilised_CLASS_h
#define EIGEN_SUPERLU_SUPPORT

#include <string.h>
#include <vector>
#include <fstream>
#include <memory>
#include "femBase.h"
#include "SolverPetsc.h"
#include "BoundaryPatch.h"
#include "BoundaryCondition.h"
#include "util.h"
#include "myMathFunction.h"
#include "myMesh.h"


using namespace std;



class femINSstabilised
{
    //private:

    public:

        //myMesh  mesh;
        shared_ptr<myMesh>  mesh;

        PetscErrorCode  errpetsc;
        PetscInt  AlgoType, SCHEME_TYPE;
        PetscInt  n_mpi_procs, this_mpi_proc;
        PetscInt  ndim, ndof, ntotdofs_local, ntotdofs_global, bodyForceTimeFunction, numBoundaryConditions;
        PetscInt  stepsMax, iterationsMax, outputFreq, fileCount;

        bool      firstIteration, convergedFlag, convergedFlagPrev, MULTIPHASEFLOW;

        double    computerTimeAssembly, computerTimeSolver, computerTimePostprocess, computerTimePattern, computerTimeTimeLoop;
        double    conv_tol, spectralRadius, timeFinal, rhsNorm, rhsNormPrev, timerVal;
        double    fluidProperties[20], bodyForce[3], stabilisationFactors[10];

        VectorXd  td, soln, solnPrev, solnPrev2, solnPrev3, solnPrev4, solnCur, solnInit, solnApplied;
        VectorXd  solnDot, solnDotPrev, solnDotCur, solnExtrap;
        VectorXd  muNodal, rhoNodal, phiNodal, etaNodal, ForceVectorExternal;

        ofstream  fout_convdata;

        vector<int>              assyForSoln, dofs_specified;
        vector<int>              elem_proc_id, node_proc_id, node_map_get_old, node_map_get_new;
        vector<vector<int> >     forAssyVecAll, globalDOFnumsAll;

        string    infilename, dirname, timeIntegrationScheme;

        vector<myMathFunction>   initialConitions;

        vector<unique_ptr<BoundaryCondition> >   boundaryConitions;

        unique_ptr<SolverPetsc>  solverPetsc;


    public:

        femINSstabilised();

        ~femINSstabilised();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  readInputGMSH(string& fname);

        int  readInitialConditions(ifstream& infile, string& line);

        int  readBodyForce(ifstream& infile, string& line);

        int  readBoundariesConditions(ifstream& infile, string& line);

        int  readSolverDetails(ifstream& infile, string& line);

        int  readTimeFunctions(ifstream& infile, string& line);

        int  readOutputDetailsPatch(ifstream& infile, string& line);

        int  deallocatePetscObjects();

        int  setBoundaryConditions();

        int  setSpecifiedDOFs(vector<vector<bool> >& NodeType);

        int  prepareInputData();

        int  readInputData(string& fname);

        int  readConfiguration(string& fname);

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

        bool  converged();

        int  saveSolution();

        int  addExternalForces(double loadFact);

        int  computeElementErrors(int);

        int  setInitialConditions();

        int  solveFullyImplicit();

        int  solveTimeStep();

        int  getVelocity(VectorXd& velocity);

        int  calcStiffnessAndResidual(int iter);

        int  calcStiffnessAndResidualImplicit2D(int iter);
        int  calcStiffnessAndResidualLinearised2D(int iter);
        int  calcStiffnessAndResidualExtrapolated2D(int iter);

        int  calcStiffnessAndResidualImplicit3D(int iter);
        int  calcStiffnessAndResidualLinearised3D(int iter);
        int  calcStiffnessAndResidualExtrapolated3D(int iter);

        int  calcForcesOnBoundaries();

        int  checkResult(string&);

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







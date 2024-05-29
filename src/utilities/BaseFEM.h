#ifndef incl_BaseFEM_h
#define incl_BaseFEM_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionDataSolid.h"
#include "headersVTK.h"
#include "ElementBase.h"
#include "MaterialBase.h"
#include "BoundaryPatch.h"
#include "SolverPetsc.h"


class BaseFEM
{
    //private:

    public:

        SolverLibrary  SOLVER_LIB_TYPE;
        SolverType   NONLINEAR_SOLVER_TYPE;

        int  itersMax, stepsMax, outputFrequency, tis, AlgoType, SCHEME_TYPE;
        int  ndim, ndof, lumpType, fileCount, numProc, numMaterials, nElemTypes;
        int  nNode, nElem, nElemFaces, npElem, IterNum, iterCount, nDBC, nFBC, totalDOF, dispDOF, presDOF, mpotDOF, tempDOF;
        int  nNode_Pres, nOutputFaceLoads;
        int  dispDegree,  presDegree, loadStepConverged;
        int  numBoundaryPathes;

        double  specRad, conv_tol, rNormPrev, rNorm;
        double  arclenIncr, arclenIncrPrev, arclenIncrMax, arclenIncrMin;
        double  loadFactor, loadFactorPrev, loadFactorPrev2;
        double  ctimFactSolvUpdt, ctimCalcStiffRes;
        //double  lamb, lambPrev, lambStep, lambIter, lambIncrStep, lambIncrIter, dl;
        double  dt, timeNow, timeFinal, rhsNorm, rhsNormPrev;

        bool  solverOK, firstIteration, DEBUG;
        bool  localStiffnessError, intVarFlag, IMPLICIT_SOLVER, MIXED_ELEMENT;
        bool  ARC_LENGTH, convergedFlagPrev, convergedFlag;
        std::ofstream  fout_explicit;

        PetscInt  n_mpi_procs, this_mpi_proc;
        PetscInt  bfs_start, bfs_end, bfs_local;
        PetscInt  row_start, row_end, ndofs_local;
        PetscInt  elem_start, elem_end;

        PetscInt  *colTemp;
        PetscScalar  *arrayTemp;

        PetscErrorCode ierr;

        vector<myPoint>  nodePosData, nodePosDataCur, nodePosDataTarget;

        vector<vector<bool> >  NodeType;
        vector<vector<int> >  boundaryNodes;
        vector<vector<int> >  midnodeData;

        vector<int>  assyForSoln, OutputNodes, pressure_nodes, assyForSolnPres, pressure_nodes_map_g2l;
        vector<int>  magnpote_nodes, assyForSolnMpot, magnpote_nodes_map_g2l;
        vector<int>  temperature_nodes, assyForSolnTemp, temperature_nodes_map_g2l;

        VectorXd  soln, solnInit, reac, MassVector, MassVectorPres, ForceVectorExternal;
        VectorXd  solnIncrStep, solnIncrIter, rhsVec;

        vector<BoundaryPatch>  BoundaryPatches;
        
        vector<vector<double> >  DirichletBCs, NeumannBCs, DerivativeBCs, InitialConds;
        vector<vector<double> >  DirichletBCs_Pres, DirichletBCs_Epot, DirichletBCs_Mpot, DirichletBCs_Temp;
        vector<vector<double> >  OutputData, nodeForcesData;
        vector<double>  loadfactorVec;

        vector<vector<int> >  elemConn;
        vector<vector<int> >  elemConnFaces;
        vector<vector<int> >  node_elem_conn;

        ElementBase  **elems;
        ElementBase  **elemsFaces;

        SolutionDataSolid  SolnData;
        GeomDataLagrange  GeomData;
        vector<MaterialBase*>  MatlDataList;

        SolverEigen *solverEigen;
        SolverPetsc *solverPetsc;

        vector<string>  outputlist_nodal, outputlist_elemental;

        string  infilename, dirname;
        ofstream  fout_convdata;

    public:

        BaseFEM();

        virtual ~BaseFEM();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        virtual void  printLogo()
        {  return; }

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  readInput(string& fname);

        void  readInputGMSH(string& fname);

        void  readInputMine(string& fname);

        void  readControlParameters(string& fname);

        void  readDimension(ifstream& infile, string& line);

        void  readNodes(ifstream& infile, string& line);

        void  readTargetShape(ifstream& infile, string& line);

        void  readSolidElements(ifstream& infile, string& line);

        void  readSurfaceElements(ifstream& infile, string& line);

        void  readPrescribedBCs(ifstream& infile, string& line);

        void  readNodalForces(ifstream& infile, string& line);

        void  readElementProps(ifstream& infile, string& line);

        void  readMaterialProps(ifstream& infile, string& line);

        void  readNodalDataForOutput(ifstream& infile, string& line);

        void  readOutputDetails(ifstream& infile, string& line);

        void  readSolverDetails(ifstream& infile, string& line);

        void  readTimeFunctions(ifstream& infile, string& line);

        void  deallocatePetscObjects();
/*
        virtual void readInputData(std::ifstream &, MyString &);

        void  calcForceVector();

        void  processForBernsteinElements();

        void  setSolverDataForSemiImplicit();

        void  prepareDataForMagneticPotential();

        void  prepareDataForTemperature();

        void  prepareDataForPressure();

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void setSolver();

        int  prepareMatrixPattern();

        int  solve();

        int  solveWithNewtonRaphson();

        int  solveWithArclength();

        int  solveExplicitStep(int*);

        int  solveElectricField(bool update_matrix, int niter);

        void copyElemInternalVariables();

        int  calcMassMatrixForExplicitDynamics();

        virtual int calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

        //int calcStiffnessAndResidualNR();
        //int calcStiffnessAndResidualAL();

        virtual int factoriseSolveAndUpdate();

        virtual void elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt);

        virtual bool converged();

        virtual bool diverging(double);

        virtual void setTimeParam();

        virtual void timeUpdate();

        virtual void updateIterStep();

        virtual void reset();

        virtual void addExternalForces();

        void  computeElementErrors(int);

        void  setInitialConditions();

        int  ModalAnalysis(int nn=10, bool flag=true, double fact=1.0);

        void  assignBoundaryConditions();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  writeNodalData();

        void  writeReadResult(int, MyString &);

        void  plotGeom();

        void  postProcess();

        void  projectFromElemsToNodes(bool, int, int, int);

        void  projectStrains(bool, int, int, int);

        void  projectStresses(bool, int, int, int);

        void  projectInternalVariables(bool, int, int, int);
*/
};





#endif







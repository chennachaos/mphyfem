#ifndef incl_GrowthFEM_CLASS_h
#define incl_GrowthFEM_CLASS_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "headersVTK.h"
#include "ElementBase.h"
#include "MaterialBase.h"
#include "BaseFEM.h"





class GrowthFEM : public BaseFEM
{
/*
    public:
        SolverLibrary  SOLVER_LIB_TYPE;
        SolverType   NONLINEAR_SOLVER_TYPE;

        int  max_iters, max_steps, outputfrequency;
        int  ndim, ndof, lumpType, filecount, numProc, numMaterials, nElemTypes;
        int  nNode, nNode_Pres, nElem, nElemFaces, npElem, IterNum, iterCount, nDBCs, totalDOF, dispDOF, presDOF;
        int  dispDegree,  presDegree, mpotDegree, tempDegree, loadStepConverged;

        double  tol, rNormPrev, rNorm,  totalError, tfinal;
        double  arclenIncr, arclenIncrPrev, arclenIncrMax, arclenIncrMin;
        double  loadfactor, loadfactorPrev, loadfactorPrev2;

        double  ctimFactSolvUpdt, ctimCalcStiffRes;
        //double  lamb, lambPrev, lambStep, lambIter, lambIncrStep, lambIncrIter, dl;

        bool  solverOK, firstIter;
        bool  localStiffnessError, intVarFlag, STAGGERED, IMPLICIT_SOLVER, MIXED_ELEMENT, COUPLED_PROBLEM;
        bool  ARC_LENGTH, convergedFlagPrev, convergedFlag;
        std::ofstream  fout_explicit;
        char VTKfilename[200];

        vector<myPoint>  nodePosData;

        vector<vector<int> >  ID, LM;
        vector<vector<int> >  boundaryNodes, forAssyMat;
        vector<vector<bool> >  NodeType;
        vector<vector<int> >  midnodeData;
        vector<vector<int> >  elemConn, node_elem_conn;

        vector<int>  assyForSoln, assyForSolnPres, pressure_nodes_map_g2l;

        VectorXd  soln, solnInit, reac, MassVector, MassVectorPres, ForceVectorExternal;

        vector<vector<double> >  DirichletBCs, NeumannBCs, DerivativeBCs, InitialConds, DirichletBCs_Pres;
        vector<vector<double> >  pointBCs, OutputData, nodeForcesData;
        vector<vector<double> >  ElemFaceLoadData;
        vector<double>  loadfactorVec;

        ElementBase  **elems;
        ElementBase  **elemsFaces;


        MyString   anlySolnType;

        SolutionData  SolnData;
        GeomDataLagrange  GeomData;
        vector<MaterialBase*>  MatlDataList;

        SolverEigen *solverEigen;

        SparseMatrixXd  spmtxPostproc;
        VectorXd  rhsPostproc;

*/
    public:

        GrowthFEM();

        ~GrowthFEM();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        void  printLogo();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  printComputerTime(bool reset = true, int detailFlg = 1);

        void  prepareElemProp();

        void  prepareMatlProp();

        void  prepareInputData();

        void  updateGeometry();

        void  readControlParameters(string& fname);
/*
        void  readInput(string& fname);

        void  readDimension(ifstream& infile, string& line);

        void  readNodes(ifstream& infile, string& line);

        void  readSolidElements(ifstream& infile, string& line);

        void  readSurfaceElements(ifstream& infile, string& line);

        void  readPrescribedBCs(ifstream& infile, string& line);

        void  readNodalForces(ifstream& infile, string& line);

        void  readElementProps(ifstream& infile, string& line);

        void  readMaterialProps(ifstream& infile, string& line);

        void  readOutputDetails(ifstream& infile, string& line);

        void  readSolverDetails(ifstream& infile, string& line);

        void  readTimeFunctions(ifstream& infile, string& line);
*/
        void  readInputData(std::ifstream &, MyString &);

        void  processForBernsteinElements();

        void  prepareDataForPressure();

        void  setSolverDataForSemiImplicit();

        void  setSolverDataForTIC();

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

        int  solveStepArcLengthGrowthModel();

        //int  solveStepDeflation(int solns_max, int iter_max, double tol1);
        int  solveStepDeflation();

        int  solveExplicitStep(int*);

        int  solveExplicitStepNIC(int*);

        int  solveExplicitStepFIC(int*);

        void copyElemInternalVariables();

        int  calcMassMatrixForExplicitDynamics();

        virtual int calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

        virtual int factoriseSolveAndUpdate();

        virtual void elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt);

        virtual bool converged();

        virtual bool diverging(double);

        virtual void setTimeParam();

        virtual void timeUpdate();

        virtual void updateIterStep();

        virtual void saveSolution();

        virtual void reset();

        virtual void addExternalForces();

        void  computerInternalForceDerivativeOfGrowth();

        void  applyBoundaryConditions();

        void  computeElementErrors(int);

        void  setInitialConditions();

        void  calcForceVector();

        void  assignBoundaryConditions();

        int  ModalAnalysis(int nn=10, bool flag=true, double fact=1.0);

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  plotGeom();

        void  postProcess();

        void  writeNodalData();

        void  writeNodalDataExplicitScheme();

        void  writeReadResult(int, MyString &);

        void  prepareMatrixPatternPostProcess();

        void  computeMatrixForPostProcess();

        void  projectFromElemsToNodes(bool, int, int, int);

        void  projectStrains(bool, int, int, int);

        void  projectStresses(bool, int, int, int);

        void  projectInternalVariables(bool, int, int, int);

        void  contourplot();

};





#endif







#ifndef incl_MagnetomechFEM_CLASS_h
#define incl_MagnetomechFEM_CLASS_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "headersVTK.h"
#include "ElementBase.h"
#include "MaterialBase.h"
#include "BaseFEM.h"


class MagnetomechFEM : public BaseFEM
{
    //private:

    public:
/*
        SolverLibrary  SOLVER_LIB_TYPE;
        SolverType   NONLINEAR_SOLVER_TYPE;

        int  max_iters, max_steps, outputfrequency;
        int  ndim, ndof, lumpType, filecount, numProc, numMaterials, nElemTypes;
        int  nNode, nElem, nElemFaces, npElem, IterNum, iterCount, nDBCs, totalDOF, dispDOF, mpotDOF, presDOF, tempDOF;
        int  nNode_Pres, nNode_Mpot, nNode_Temp;
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

        vector<vector<int> >  IEN, IEN2, ID, LM;
        vector<vector<int> >  boundaryNodes, forAssyMat;
        vector<vector<bool> >  NodeType;
        vector<vector<int> >  midnodeData;

        vector<int>  assyForSoln, pressure_nodes, assyForSolnPres, pressure_nodes_map_g2l;
        vector<int>  magnpote_nodes, assyForSolnMpot, magnpote_nodes_map_g2l;
        vector<int>  temperature_nodes, assyForSolnTemp, temperature_nodes_map_g2l;

        VectorXd  soln, solnInit, reac, MassVector, MassVectorPres, ForceVectorExternal;
        VectorXd  solnIncrStep, solnIncrIter, rhsVec;

        vector<vector<double> >  DirichletBCs, NeumannBCs, DerivativeBCs, InitialConds;
        vector<vector<double> >  DirichletBCs_Pres, DirichletBCs_Mpot, DirichletBCs_Temp;
        vector<vector<double> >  pointBCs, OutputData, nodeForcesData;
        vector<double>  loadfactorVec;

        vector<vector<int> >  elemConn;
        vector<vector<int> >  elemConnFaces;
        vector<vector<int> >  node_elem_conn;

        ElementBase  **elems;
        ElementBase  **elemsFaces;

        SolutionData  SolnData;
        GeomDataLagrange  GeomData;
        vector<MaterialBase*>  MatlDataList;

        SolverEigen *solverEigen;

        vector<string>  outputlist_nodal, outputlist_elemental;

        SparseMatrixXd  spmtxPostproc, spmtxElecField;
        VectorXd  rhsPostproc, rhsElecField;

        vtkSmartPointer<vtkFloatArray>           scaVTK;
        vtkSmartPointer<vtkFloatArray>           presVTK;

        SimplicialLDLT<SparseMatrix<double> > solverElecField;
        //SuperLU<SparseMatrixXd > solverElecField;
        //ConjugateGradient<SparseMatrixXd, Lower|Upper >   solverElecField;
        //ConjugateGradient<SparseMatrixXd, Lower|Upper, IncompleteLUT<double>  >   solverElecField;
        //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverElecField;
*/

    public:

        MagnetomechFEM();

        ~MagnetomechFEM();

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
        virtual  void  readInput(string& fname);

        virtual  void  readDimension(ifstream& infile, string& line);

        virtual  void  readNodes(ifstream& infile, string& line);

        virtual  void  readSolidElements(ifstream& infile, string& line);

        virtual  void  readSurfaceElements(ifstream& infile, string& line);

        virtual  void  readPrescribedBCs(ifstream& infile, string& line);

        virtual  void  readNodalForces(ifstream& infile, string& line);

        virtual  void  readElementProps(ifstream& infile, string& line);

        virtual  void  readMaterialProps(ifstream& infile, string& line);

        virtual  void  readOutputDetails(ifstream& infile, string& line);

        virtual  void  readSolverDetails(ifstream& infile, string& line);

        virtual  void  readTimeFunctions(ifstream& infile, string& line);
*/
        virtual void  readInputData(std::ifstream &, MyString &);

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

        virtual void saveSolution();

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

        void  projectStrains(bool, int, int, int, VectorXd&  output);

        void  projectStresses(bool, int, int, int, VectorXd&  output);

        void  projectInternalVariables(bool, int, int, int, VectorXd&  output);

};





#endif







#ifndef incl_ElectromechFEM_CLASS_h
#define incl_ElectromechFEM_CLASS_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "headersVTK.h"
#include "ElementBase.h"
#include "MaterialBase.h"


class PropertyItem;

#define ELECTROMECH_ELEMENT_TYPE_NAMES {                                      \
                                 "BernsteinElem2DElecMechTri11",    /*  0 - 1001 */ \
                                 "BernsteinElem2DElecMechTri22",    /*  1 - 1002 */ \
                                 "BernsteinElem2DElecMechTri22F",   /*  2 - 1003 */ \
                                 "BernsteinElem2DElecMechTri210",   /*  3 - 1004 */ \
                                 "BernsteinElem2DElecMechTri220",   /*  4 - 1005 */ \
                                 "BernsteinElem2DElecMechTri211",   /*  5 - 1006 */ \
                                 "BernsteinElem2DElecMechTri221",   /*  6 - 1007 */ \
                                 "BernsteinElem2DElecMechQua11",    /*  7 - 1008 */ \
                                 "BernsteinElem2DElecMechQua11F",   /*  8 - 1009 */ \
                                 "BernsteinElem2DElecMechQua110",   /*  9 - 1010 */ \
                                 "BernsteinElem2DElecMechQua1100",  /* 10 - 1011 */ \
                                 "BernsteinElem2DElecMechQua221",   /* 11 - 1012 */ \
                                 "BernsteinElem2DElecMechQua221P",  /* 12 - 1013 */ \
                                 "dummy",                           /* 13 - 1014 */ \
                                 "dummy",                           /* 14  */ \
                                 "dummy",                           /* 15  */ \
                                 "dummy",                           /* 16  */ \
                                 "dummy",                           /* 17  */ \
                                 "dummy",                           /* 18  */ \
                                 "dummy",                           /* 19  */ \
                                 "dummy",                           /* 20  */ \
                                 "dummy",                           /* 21  */ \
                                 "dummy",                           /* 22  */ \
                                 "dummy",                           /* 23  */ \
                                 "dummy",                           /* 24  */ \
                                 "dummy",                           /* 25  */ \
                                 "dummy",                           /* 26  */ \
                                 "dummy",                           /* 27  */ \
                                 "dummy",                           /* 28  */ \
                                 "dummy",                           /* 29  */ \
                                 "dummy",                           /* 30  */ \
                                 "dummy",                           /* 31  */ \
                                 "dummy",                           /* 32  */ \
                                 "dummy",                           /* 33  */ \
                                 "dummy",                           /* 34  */ \
                                 "dummy",                           /* 35 */ \
                                 "dummy",                           /* 36  */ \
                                 "dummy",                           /* 37  */ \
                                 "dummy",                           /* 38  */ \
                                 "dummy",                           /* 39  */ \
                                 "BernsteinElem3DElecMechTet11",    /* 40 - 1041 */ \
                                 "BernsteinElem3DElecMechTet22",    /* 41 - 1042 */ \
                                 "BernsteinElem3DElecMechTet22F",   /* 42 - 1043 */ \
                                 "BernsteinElem3DElecMechTet210",   /* 43 - 1044 */ \
                                 "BernsteinElem3DElecMechTet220",   /* 44 - 1045 */ \
                                 "BernsteinElem3DElecMechTet211",   /* 45 - 1046 */ \
                                 "BernsteinElem3DElecMechTet221",   /* 46 - 1047 */ \
                                 "BernsteinElem3DElecMechHex11",    /* 47 - 1048 */ \
                                 "BernsteinElem3DElecMechHex11F",   /* 48 - 1049 */ \
                                 "BernsteinElem3DElecMechHex110",   /* 49 - 1050 */ \
                                 "BernsteinElem3DElecMechHex1100",  /* 50 - 1051 */ \
                                 "BernsteinElem3DElecMechHex221",   /* 51 - 1052 */ \
                                 "BernsteinElem3DElecMechHex221P",  /* 52 - 1053 */ \
                                 "dummy",                           /* 53 - 1054 */ \
                                 "dummy",                           /* 54  */ \
                                 "dummy",                           /* 55  */ \
                                 "dummy",                           /* 56  */ \
                                 "dummy",                           /* 57  */ \
                                 "dummy",                           /* 58  */ \
                                 "dummy",                           /* 59  */ \
                                 "dummy",                           /* 60  */ \
                                 "dummy",                           /* 61  */ \
                                 "dummy",                           /* 62  */ \
                                 "dummy",                           /* 63  */ \
                                 "dummy",                           /* 64  */ \
                                 "dummy",                           /* 65  */ \
                                 "dummy",                           /* 66  */ \
                                 "dummy",                           /* 67  */ \
                                 "dummy",                           /* 68  */ \
                                 "dummy",                           /* 69  */ \
                                 "dummy",                           /* 70  */ \
                                 "dummy",                           /* 71  */ \
                                 "dummy",                           /* 72  */ \
                                 "dummy",                           /* 73  */ \
                                 "dummy",                           /* 74  */ \
                                 "dummy",                           /* 75  */ \
                                 "dummy",                           /* 76  */ \
                                 "dummy",                           /* 77  */ \
                                 "dummy",                           /* 78  */ \
                                 "dummy",                           /* 79  */ \
                                 "BernsteinElem2DEdge2Node",        /* 80  */ \
                                 "BernsteinElem2DEdge3Node",        /* 81  */ \
                                 "BernsteinElem3DFaceTria6Node",    /* 82  */ \
                                 NULL};




class ElectromechFEM
{
    //private:

    public:

        PhysicsType  PHYSICS_TYPE;
        SolverType   SOLVER_TYPE;

        int  ndim, ndof, tis, lumpType, filecount, numProc;
        int  nNode, nElem, nElemFaces, npElem, IterNum, iterCount, totalDOF, dispDOF, epotDOF, presDOF, tempDOF;
        int  nNode_Pres, nNode_Epot, nNode_Temp;
        int  dispDegree,  presDegree, epotDegree, tempDegree;

        double  tol, td[100];
        double  rNormPrev, rNorm,  ctimFactSolvUpdt, ctimCalcStiffRes, totalError, rhoInfty;
        double  lamb, lambPrev, lambStep, lambIter, lambIncrStep, lambIncrIter, dl;

        bool  solverOK, firstIter, localStiffnessError, intVarFlag, STAGGERED, IMPLICIT_SOLVER, MIXED_ELEMENT, ARCLENGTH;
        char  VTKfilename[500];

        vector<myPoint>  nodePosData;

        vector<vector<int> >  IEN, IEN2, ID, LM;
        vector<vector<int> >  boundaryNodes, forAssyMat;
        vector<vector<bool> >  NodeType;
        vector<vector<int> >  midnodeData;

        vector<int>  assyForSoln, pressure_nodes, assyForSolnPres, pressure_nodes_map_g2l;
        vector<int>  elecpote_nodes, assyForSolnEpot, elecpote_nodes_map_g2l;
        vector<int>  temperature_nodes, assyForSolnTemp, temperature_nodes_map_g2l;

        VectorXd  soln, solnInit, reac, MassVector, MassVectorPres, ForceVectorExternal;
        VectorXd  solnIncrStep, solnIncrIter, rhsVec;

        vector<vector<double> >  DirichletBCs, NeumannBCs, DerivativeBCs, InitialConds;
        vector<vector<double> >  DirichletBCs_Pres, DirichletBCs_Epot, DirichletBCs_Temp;
        vector<vector<double> >  pointBCs, OutputData, nodeForcesData;
        vector<vector<double> >  ElemFaceLoadData;

        vector<vector<int> >  elemConn;
        vector<vector<int> >  node_elem_conn;

        ElementBase  **elems;
        ElementBase  **elemsFaces;

        SolutionData  SolnData;
        GeomDataLagrange  GeomData;
        MaterialBase  *MatlData;

        SolverEigen *solverEigen;

        SparseMatrixXd  spmtxPostproc, spmtxElecField;
        VectorXd  rhsPostproc, rhsElecField;

        vtkSmartPointer<vtkFloatArray>           scaVTK;
        vtkSmartPointer<vtkFloatArray>           presVTK;

        SimplicialLDLT<SparseMatrix<double> > solverElecField;
        //SuperLU<SparseMatrixXd > solverElecField;
        //ConjugateGradient<SparseMatrixXd, Lower|Upper >   solverElecField;
        //ConjugateGradient<SparseMatrixXd, Lower|Upper, IncompleteLUT<double>  >   solverElecField;
        //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverElecField;


    public:

        ElectromechFEM();

        ~ElectromechFEM();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        int  solveStep(int niter);

        int  solveExplicitStep(int*);

        int  solveElectricField(bool update_matrix, int niter);

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        virtual void printComputerTime(bool reset = true, int detailFlg = 1);

        void assignBoundaryConditions();

        void  prepareElemProp();
        void  prepareMatlProp();

        virtual void prepareInputData();

        void  updateGeometry();

        virtual void readInputData(std::ifstream &, MyString &);

        void  plotGeom(int, bool, int, bool, int*);

        void  calcForceVector();

        void  writeNodalData();

        void  writeReadResult(int, MyString &);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////


        virtual void setSolver(int, int *parm = NULL, bool cIO = false);

        int  prepareMatrixPattern();

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

        int  solveWithArclength();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  projectFromElemsToNodes(bool, int, int, int);

        void  projectStrains(bool, int, int, int);

        void  projectStresses(bool, int, int, int);

        void  projectInternalVariables(bool, int, int, int);

        void  postProcess(int, int, int, bool, double, double, int*);

        void  computeTotalBodyForce(int );

        void  processForBernsteinElements();

        void  setSolverDataForSemiImplicit();

        void  prepareDataForElectricPotential();

        void  prepareDataForTemperature();

        void  prepareDataForPressure();
};





#endif







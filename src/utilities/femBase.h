#ifndef incl_femBase_h
#define incl_femBase_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "headersVTK.h"
#include "ElementBase.h"
#include "MaterialBase.h"
#include "ElementTypeData.h"
#include "Domain.h"
#include "BoundaryPatch.h"
#include "BoundaryCondition.h"
#include "TractionCondition.h"
#include "NodalForces.h"


class femBase
{
    //private:

    public:

        SolverLibrary  SOLVER_LIB_TYPE;
        SolverType   NONLINEAR_SOLVER_TYPE;

        string  inputfilename, dirname, restartfilename;
        PetscErrorCode  errpetsc;
        int  ntotdofs_local, ntotdofs_global, n_mpi_procs, this_mpi_proc, nElem_global;
        int  elem_start, elem_end, nElem_local, nNode_local;
        int  row_start, row_end, node_start, node_end, nNode_owned;
        int  iterationsMax, stepsMax, outputfrequency, outputfreq_vtk;
        int  ndim, ndof, lumpType, filecount, numProc, numMaterials, numElemTypes;
        int  numPhysicalNames, numBoundaryPatches, numDomains, numBoundaryConditions, numTractionConditions;
        int  IterNum, iterCount, dispDOF, mpotDOF, presDOF, tempDOF;
        int  nNode_global, nNode_Pres, nNode_Mpot, nNode_Temp;
        int  dispDegree,  presDegree, mpotDegree, tempDegree, loadStepConverged;

        double  conv_tol, rhsNormPrev, rhsNorm,  totalError, timeFinal;
        double  arclenIncr, arclenIncrPrev, arclenIncrMax, arclenIncrMin;
        double  loadFactor, loadFactorPrev, loadFactorPrev2;

        double  ctimFactSolvUpdt, ctimCalcStiffRes;
        //double  lamb, lambPrev, lambStep, lambIter, lambIncrStep, lambIncrIter, dl;

        bool  solverOK, firstIteration, MIXED_ELEMENT_P0;
        bool  localStiffnessError, intVarFlag, IMPLICIT_SOLVER, MIXED_ELEMENT, MIXED_STAB_ELEMENT, COUPLED_PROBLEM;
        bool  ARC_LENGTH, convergedFlagPrev, convergedFlag;
        std::ofstream  fout_explicit;

        vector<int>  dofs_specified_disp, dofs_specified_pres;
        vector<int>  elem_proc_id, node_proc_id, node_map_get_old, node_map_get_new;

        vector<vector<int> >  boundaryNodes;
        vector<vector<int> >  midnodeData;

        vector<int>  assyForSoln, pressure_nodes, assyForSolnPres, pressure_nodes_map_g2l;
        vector<int>  magnpote_nodes, assyForSolnMpot, magnpote_nodes_map_g2l;
        vector<int>  temperature_nodes, assyForSolnTemp, temperature_nodes_map_g2l;

        VectorXd  soln, solnInit, reac, MassVector, MassVectorPres, ForceVectorExternal;
        VectorXd  solnIncrStep, solnIncrIter, rhsVec;

        vector<vector<double> >  DirichletBCs, NeumannBCs, DerivativeBCs, InitialConds;
        vector<vector<double> >  DirichletBCs_Pres, DirichletBCs_Epot, DirichletBCs_Mpot, DirichletBCs_Temp;
        vector<vector<double> >  pointBCs, OutputData, nodeForcesData, NodalDataOutput;
        vector<double>  loadFactorVec;

        vector<vector<int> >  elemConn, elemConnFaces, node_elem_conn;


        ElementBase  **elems;

        SolutionDataSolid  SolnData;
        GeomDataLagrange  GeomData;

        vector<unique_ptr<BoundaryPatch> >       BoundaryPatches;
        vector<unique_ptr<BoundaryCondition> >   boundaryConitions;
        vector<unique_ptr<TractionCondition> >   tractionConitions;
        vector<unique_ptr<Domain> >              Domains;
        vector<string>                           initialConitions;
        //vector<myMathFunction>              initialConitions;

        NodalForces                              nodalForces;

        vector<MaterialBase*>      MatlDataList;
        vector<ElementTypeData*>   ElementTypeDataList;


        unique_ptr<SolverEigen>  solverEigen;
        unique_ptr<SolverPetsc>  solverPetsc;

        vector<string>  outputlist_nodal, outputlist_elemental;

        ofstream fout_nodaldata;

    public:

        femBase();

        virtual ~femBase();

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

        void  readFiles(ifstream& infile, string& line);

        void  readModelType(ifstream& infile, string& line);

        void  readDomainsData(ifstream& infile, string& line);

        void  readElementProps(ifstream& infile, string& line);

        void  readMaterialProps(ifstream& infile, string& line);

        void  readSolverDetails(ifstream& infile, string& line);

        void  readTimeFunctions(ifstream& infile, string& line);

        void  readOutputDetails(ifstream& infile, string& line);

        void  readBodyForce(ifstream& infile, string& line);

        void  readPrescribedBCs(ifstream& infile, string& line);

        void  readTractions(ifstream& infile, string& line);

        void  readControlParameters(string& fname);

        void  readNodalForces(ifstream& infile, string& line);

        void  readNodalDataOutputDetails(ifstream& infile, string& line);

        void  readInitialConditions(ifstream& infile, string& line);

        void  setSolver(int slv);


/*
        void  processForBernsteinElements();

        void  setSolverDataForSemiImplicit();

        void  prepareDataForMagneticPotential();

        void  prepareDataForTemperature();

        void  prepareDataForPressure();
*/
        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

/*
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
*/
        int  deallocatePetscObjects();

        void  writeNodalData();

        void  writeReadResult(int, string &);

        virtual  void  plotGeom()
        { cout << "   'plotGeom' is not defined for this element!\n\n"; return; }

        virtual  void  postProcess()
        { cout << "   'postProcess' is not defined for this element!\n\n"; return; }

        void  projectFromElemsToNodes(bool, int, int, int);

        void  projectStrains(bool, int, int, int);

        void  projectStresses(bool, int, int, int);

        void  projectInternalVariables(bool, int, int, int);
};





#endif







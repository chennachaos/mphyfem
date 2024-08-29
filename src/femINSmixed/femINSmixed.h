
#ifndef incl_femINSmixed_CLASS_h
#define incl_femINSmixed_CLASS_h
#define EIGEN_SUPERLU_SUPPORT

#include <string.h>
#include <vector>
#include <fstream>
#include "ImmersedSolid.h"
#include "ElementBase.h"
//#include "headersEigen.h"
#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include "ImmersedIntegrationElement.h"
#include "SolutionData.h"
#include "FluidMaterialBase.h"

#include "Domain.h"
#include "BoundaryPatch.h"
#include "ElementTypeData.h"
#include "BoundaryCondition.h"
#include "TractionCondition.h"
#include "NodalForces.h"


using std::vector;
using std::cout;
using std::string;
using std::ofstream;

class ElementBase;


#define BERNSTEIN_ELEMENT_TYPE_NAMES {\
                                "BernsteinElem2DPoissonTria6Node", \
                                "BernsteinElem2DSolidTria6Node", \
                                "BernsteinElem2DSolidBbarFbarTria6Node", \
                                "BernsteinElem2DSolidMixedTria6Node", \
                                "BernsteinElem3DPoissonTetra10Node", \
                                "BernsteinElem3DSolidTetra10Node", \
                                "BernsteinElem3DSolidBbarFbarTetra10Node", \
                                "BernsteinElem3DSolidMixedTetra10Node", \
                                "BernsteinElem2DPoissonQuad9Node", \
                                "BernsteinElem3DPoissonHex27Node", NULL};


#define BERNSTEIN_ELEMENT_TYPE_NAMES_FACE {\
                                "BernsteinElem2DEdge3Node", \
                                "BernsteinElem3DFaceTria6Node", \
                                "BernsteinElem3DFaceQuad9Node", NULL};

enum SchemeType
{
  SCHEME_TYPE_IMPLICIT = 0, SCHEME_TYPE_SEMIIMPLICIT = 1, SCHEME_TYPE_EXPLICIT = 2
};


class femINSmixed
{
    //private:

    public:

        int  ndim, ndof, numMaterials, numBoundaryPatches, numDomains, bodyForceTimeFunction, numBoundaryConditions, numTractionConditions;
        int  stepsMax, iterationsMax, outputFreq, fileCount;
        int  nElem_global, nNode_global, nElem_local, nNode_local, nNode_owned;
        int  n_mpi_procs, this_mpi_proc;
        int  node_start, node_end, elem_start, elem_end;
        int  row_start, row_end, ntotdofs_local, ntotdofs_global;

        double    computerTimeAssembly, computerTimeSolver, computerTimePostprocess, computerTimePattern, computerTimeTimeLoop;
        double    fluidProperties[10], bodyForce[3];
        double    conv_tol, timeFinal, rhsNorm, rhsNormPrev;

        VectorXd  td, reacVec;

        vector<int>              assyForSoln, dofs_specified_velo, dofs_specified_pres;
        vector<int>              elem_proc_id, node_proc_id, node_map_get_old, node_map_get_new;
        vector<vector<int> >     elemConn, node_elem_conn;                  //!< element-node connectivity array
        vector<vector<int> >     NodeDofArray, forAssyVecAll, globalDOFnumsAll;
        vector<vector<bool> >    NodeType;
        vector<vector<int> >     midnodeData;

        string    infilename, dirname, timeIntegrationScheme, ELEMENT_TYPE;

        vector<unique_ptr<BoundaryPatch> >       BoundaryPatches;
        vector<unique_ptr<BoundaryCondition> >   boundaryConitions;
        vector<unique_ptr<TractionCondition> >   tractionConitions;
        vector<unique_ptr<Domain> >              Domains;
        vector<string>                           initialConitions;

        NodalForces                              nodalForces;

        vector<ElementTypeData*>  ElementTypeDataList;

        unique_ptr<SolverEigen>  solverEigen;
        //unique_ptr<SolverPetsc>  solverPetsc;

        int  dispDOF, presDOF, lambdaDOF, totalDOF;
        int  totalDOF_Velo, totalDOF_Pres, totalDOF_Lambda, totalDOF_Solid;
        int  nNode_Velo, nNode_Pres, nDBC, nDBC_Velo, nDBC_Pres, nFBC, nOutputFaceLoads;
        int  npElemVelo, npElemPres;
        int  AlgoType, nImmersedElems, nImmersedNodes, SCHEME_TYPE;
        
        bool firstIteration, GRID_CHANGED, IB_MOVED;

        double  timeNow, loadFactor;
        double  spectralRadius, am, gamm1, gamm2, CFL, dt, amDgammaDt;
        double  elemData[50], timeData[50];

        vector<myPoint>          nodeCoords, nodeCoordsCur; //!< coordinates of the nodes
        vector<vector<int> >     midNodeData;               //!< data for processing the middle nodes
        vector<vector<int> >     outputEdges;               //!< data for computing drag/lift forces
        vector<vector<double> >  immersed_node_coords;      //!< coordinates of the nodes of immersed boundary

        vector<int>  assyForSolnVelo, assyForSolnPres, OutputNodes;
        vector<int>  pressure_nodes, pressure_nodes_map;


        vector<vector<double> >  DirichletBCsVelo;          //!< Dirichlet BCs for velocity
        vector<vector<double> >  DirichletBCsPres;          //!< Dirichlet BCs for pressure
        vector<vector<double> >  NeumannBCs;                //!< Neumann BCs
        vector<vector<double> >  InitialConds;              //!< Initial conditions
        vector<vector<double> >  OutputData;                //!< data for output
        vector<vector<double> >  nodeForcesData;
        vector<vector<double> >  ElemFaceLoadData;
        vector<vector<double> >  fluidProps;

        SolutionData  SolnData;
        vector<FluidMaterialBase*>  FluidMatlDataList;

        ElementBase  **elems;

        vector<ImmersedSolid*>  ImmersedBodyObjects;
        vector<ImmersedIntegrationElement*>  ImmersedElements;

        VectorXd  soln, solnInit, ForceVectorExternal;
        VectorXd  totalForce, totalMoment, centroid;
        VectorXd  force, forceCur, forcePrev, forcePrev2;
        VectorXd  pres, presCur, presPrev, presPrev2, presPrev3, presIncr;
        VectorXd  presDot, presDotPrev, presDotCur, presDiff;
        VectorXd  velo, veloCur, veloDiff, veloPrev, veloPrev2, veloPrev3, veloIncr;
        VectorXd  veloDot, veloDotPrev, veloDotCur;
        VectorXd  veloApplied, presApplied;
        VectorXd  globalMassVelo, globalMassPres;
        VectorXd  rhsVecVelo, rhsVecPres, rhsVec;
        VectorXd  lambdas, lambdasPrev, lambdasIncr, lambdasCur;

        VectorXd  rhsVecVeloTemp, rhsVecPresTemp;

        SparseMatrixXd  matK, matMuuInv, matKup, matKpu, matSchur;

        //SimplicialLLT<SparseMatrix<double> > solverSchur;
        //SimplicialLDLT<SparseMatrix<double> > solverSchur;
        //SparseLU<SparseMatrixXd > solverSchur;
        SuperLU<SparseMatrixXd > solverSchur;
        //ConjugateGradient<SparseMatrixXd, Lower|Upper >   solverSchur;
        //ConjugateGradient<SparseMatrixXd, Lower|Upper, IncompleteLUT<double>  >   solverSchur;
        //BiCGSTAB<SparseMatrixXd> solverSchur;
        //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverSchur;

        //SimplicialLDLT<SparseMatrix<double> > solver;
        //SparseLU<SparseMatrixXd > solver;
        SuperLU<SparseMatrixXd > solver;

        int   PT[64], IPARM[64];
        int   phase, error, SOLVER, MTYPE, MAXFCT, MNUM, NRHS, MSGLVL;

        double DPARM[64], ddum, *array;

        vector<int>  csr, col, perm;

    public:

        femINSmixed();

        ~femINSmixed();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        int  readInputGMSH(string& fname);

        int  readConfiguration(string& fname);

        int  readDomains(ifstream& infile, string& line);

        int  readFluidProperties(ifstream& infile, string& line);

        int  readInitialConditions(ifstream& infile, string& line);

        int  readBodyForce(ifstream& infile, string& line);

        int  readBoundaryConditions(ifstream& infile, string& line);

        int  readTractionConditions(ifstream& infile, string& line);

        int  readNodalForces(ifstream& infile, string& line);

        int  readSolverDetails(ifstream& infile, string& line);

        int  readTimeFunctions(ifstream& infile, string& line);

        int  readElementProps(ifstream& infile, string& line);

        int  readPrescribedBCs(ifstream& infile, string& line);

        int  readFluidProps(ifstream& infile, string& line);

        int  readNodalDataForOutput(ifstream& infile, string& line);

        int  readOutputDetails(ifstream& infile, string& line);

        int  readImmersedSolids(string& fname);

        int  readOutputDetailsPatch(ifstream& infile, string& line);

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  setSpecifiedDOFs_Velocity(vector<vector<bool> >&  NodeDofType);

        int  setBoundaryConditions();

        int  addBoundaryConditions();

        int  prepareInputData();

        int  prepareFluidMesh();

        int  prepareImmersedSolids();

        int  printInfo();

        int  plotGeom(int, bool, int, bool, int*);

        int  applyExternalForces();

        int  writeNodalData();

        int  writeReadResult(int, string&);

        int  findElementNumber(myPoint& pt);

        int  setLoadFactor(double fact)
        {
            loadFactor = fact;    return 0;
        }

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int   setSolver(int);

        int   prepareMatrixPattern();

        int   calcMassMatrixForExplicitDynamics();

        int   setSolverDataForSemiImplicit();

        int   setSolverDataForFullyImplicit();

        int   updateMatricesSemiImplicit();

        int   solveSchurComplementProblem();

        int   applyInterfaceTerms2D();

        bool  converged();

        bool  diverging(double);

        int  setTimeParam();

        int  timeUpdate();

        int  updateIterStep();

        int  reset();
        
        int  saveSolution();

        int  addExternalForces();

        int  computeElementErrors(int);

        int  setInitialConditions();

        int  solve();

        int  solveSemiImplicit();

        int  solveFullyImplicit();

        int  solveStep(int max_iters);
        
        int  calcStiffnessAndResidual();
        
        int  factoriseSolveAndUpdate();

        int  initialise_pardiso();

        int  factoriseAndSolve_pardiso();

        int  updateImmersedSolid();

        int  updateImmersedSolid(int id, int ndof, VectorXd&, VectorXd&);

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  postProcess();
        int  postProcessPressure();
        int  postProcessVelocity();

        int  writeReadResult(int index, string& filename, int stride);

        int  writeOutputDataPatches();
};






#endif







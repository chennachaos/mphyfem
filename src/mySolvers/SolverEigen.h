
#ifndef incl_SolverEigen_h
#define incl_SolverEigen_h


#include "headersEigen.h"
#include <vector>


enum { PARDISO_STRUCT_SYM, PARDISO_UNSYM };
enum { SOLVER_EMPTY, PATTERN_OK, INIT_OK, ASSEMBLY_OK, FACTORISE_OK };


using std::vector;


class SolverEigen
{
  public:

    VectorXd   rhsVec, rhsVec2, rhsVec3, soln, solnPrev, solnTot, var1, var1Prev, var3, var3Prev, var2, var2Prev;
    VectorXd   f1, f2, f3, reac, Fext, FintDerGrowth;
    vector<int>  assyForSoln, assyForSolnPres;

    SparseMatrixXd  mtx, matKuu, matKup, matKpu, matKpp, PreCondSchur, matKuuInv, matSchur;
    SparseMatrixXd  matA, matB, matC, matD, matAinv, matS;


    int nRow, nCol, nnz, currentStatus, algoType, nU, nP, nL, nS, update_precond, count, dispDOF, presDOF;

    bool checkIO, STABILISED, update_factorisation, MIXED_ELEMENT;

    double normPrev, normCur, normRef; // norm of solution error

    SimplicialLDLT<SparseMatrix<double> > solverSchur;
    //SuperLU<SparseMatrixXd > solverSchur;
    //ConjugateGradient<SparseMatrixXd, Lower|Upper >   solverSchur;
    //ConjugateGradient<SparseMatrixXd, Lower|Upper, IncompleteLUT<double>  >   solverSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverSchur;


    ////////////////////////////
    //
    // member functions
    //
    ///////////////////////////

    SolverEigen();

    virtual ~SolverEigen();

    virtual int initialise(int p1 = 0, int p2 = 0, int p3 = 0);

    int setSolverAndParameters();

    void setAlgorithmType(int tt)
    {  algoType = tt; return; }

    virtual void printMatrixPatternToFile();

    virtual void zeroMtx();

    virtual int free();

    virtual void printInfo();

    virtual void printMatrix(int dig=8, int dig2=4, bool gfrmt=true, int indent = 0, bool interactive = false);

    virtual int assembleMatrixAndVector(int r1, int c1, vector<int>& row, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorMixed(int var2_offset, vector<int>& forAssyVar1, vector<int>& forAssyVar2, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Ru, VectorXd& Rp);

    virtual int assembleMatrixAndVector2field(int var2_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, VectorXd& Flocal1, VectorXd& Flocal2, vector<int>& forAssyVar1, vector<int>& forAssyVar2);

    virtual int assembleMatrixAndVector3field(int var2_offset, int var3_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, MatrixXd& K13, MatrixXd& K31, MatrixXd& K33, VectorXd& Flocal1, VectorXd& Flocal2, VectorXd& Flocal3, vector<int>& forAssyVar1, vector<int>& forAssyVar2, vector<int>& forAssyVar3);

    virtual int assembleMatrixAndVectorUPexpl(int var2_offset, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Rf, VectorXd& Rp, vector<int>& forAssyVar1, vector<int>& forAssyVar2);

    virtual int assembleMatrixAndVectorElecField(vector<int>& row, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();

    void setupMatricesAndVectors();

    void  computeConditionNumber();

    int  solverSchurExplicitDirect(int frequency);

    virtual int solveArclengthSystem(int loadStep, VectorXd&  Du, double& Dl, double&  Ds, double& dl);

    virtual int solveArclengthSystemGrowthModel(int loadStep, VectorXd&  DuFull, double& Dl, double&  Ds, double& dl);

    int  solverSchurExplicitCG();



    int  solverSchurExplicitDirect();


    void solverSchurCG();

    void solverSchurGMRES();

    void solverSchurBiCGSTAB();

    void solverUzawaType1();

    void solverUzawaType2();

    void resetPrecondFlag()
    {
      update_precond = 1;
    }

    void  updatePreconditioner();

};



#endif


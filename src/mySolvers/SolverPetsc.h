
#ifndef incl_SolverPetsc_h
#define incl_SolverPetsc_h

#include "SolverEigen.h"
#include "petscksp.h"
#include "petscmat.h"
#include <Eigen/Dense>
#include <vector>

using namespace std;
using namespace Eigen;



class SolverPetsc
{
  public:

    Vec  rhsVec, solnVec, reacVec;
    Mat  mtx; // linear system matrix
    KSP  ksp; // linear solver context
    PC   pc; // preconditioner context
    vector<int>  assyForSoln, assyForSolnPres;


    PetscInt nRow, nCol, nnz, dispDOF, presDOF;

    MatInfo info;

    int currentStatus;

    bool  checkIO;

    PetscReal norm; // norm of solution error

    PetscErrorCode ierr, errpetsc;

    PetscMPIInt size;

    //PetscViewer    viewer_matx, viewer_vect;

    VectorXd   Fext;

    ////////////////////////////
    //
    // member functions
    //
    ///////////////////////////


    SolverPetsc();

    virtual ~SolverPetsc();

    virtual int initialise(int size_local, int size_global, int* diag_nnz, int* offdiag_nnz);

    int setSolverAndParameters();

    virtual bool isChildOfSolverSparse(void) { return false; }

    virtual int zeroMtx();

    virtual int free();

    virtual int printInfo();

    virtual int printMatrix(int dig=8, int dig2=4, bool gfrmt=true, int indent = 0, bool interactive = false);

    virtual double giveMatrixCoefficient(int,int);

    virtual int assembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVector(int r1, int c1, vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleVector(int r1, int c1, vector<int>& row, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorSerial(vector<int>& forAssyElem, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorParallel(vector<int>& forAssyElem, vector<int>& dof_map, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVectorMixedFormulation(int r1, int c1, vector<int>& vec1, vector<int>& vec2, MatrixXd& Klocal, VectorXd& Flocal);

    virtual int assembleMatrixAndVector2field(int var2_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, VectorXd& Flocal1, VectorXd& Flocal2, vector<int>& forAssyVar1, vector<int>& forAssyVar2);

    virtual int assembleMatrixAndVector3field(int var2_offset, int var3_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, MatrixXd& K13, MatrixXd& K31, MatrixXd& K33, VectorXd& Flocal1, VectorXd& Flocal2, VectorXd& Flocal3, vector<int>& forAssyVar1, vector<int>& forAssyVar2, vector<int>& forAssyVar3);

    virtual int factorise();

    virtual int solve();

    virtual int factoriseAndSolve();

    virtual int solveArclengthSystem(int loadStep, VectorXd&  Du, double& Dl, double&  Ds, double& dl);

};



#endif

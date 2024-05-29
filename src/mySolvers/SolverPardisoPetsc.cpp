
#include <iostream>

#include "SolverPardisoPetsc.h"
#include "FunctionsSolver.h"
#include "util.h"

#include "petscksp.h"
#include "petscmat.h"
#include "petscsys.h"


using namespace std;




SolverPardisoPetsc::SolverPardisoPetsc()
{
  //MatCreate(PETSC_COMM_WORLD, &mtx);

  //VecCreate(PETSC_COMM_WORLD, &soln);
  //VecCreate(PETSC_COMM_WORLD, &solnPrev);
  //VecCreate(PETSC_COMM_WORLD, &rhsVec);
  //VecCreate(PETSC_COMM_WORLD, &reac);
}




SolverPardisoPetsc::~SolverPardisoPetsc()
{
  phase = -1; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

}




int SolverPardisoPetsc::initialise(int numProc, int matrixType, int nrow1)
{
  nRow = nCol = nrow1;

  rhsTemp.resize(nRow, 0.0);
  solnTemp.resize(nRow, 0.0);

  if (currentStatus != PATTERN_OK)
  {
    cerr << "SolverPardisoPetsc::initialise ... prepare matrix pattern first!" << endl;
    return 2;
  }

  phase = 11; error = 0;

  int idum, mtxType = 11;

  char *tmp;

  switch (matrixType)
  {
    case PARDISO_STRUCT_SYM: mtxType =  1; break; // real and structurally symmetric

    case PARDISO_UNSYM:      mtxType = 11; break; // real and unsymmetric

    default:                 cerr << "SolverPardisoPetsc::initialise ... matrix assumed to be unsymmetric!" << endl;
  }

  SOLVER = 0;       // sparse direct solver
  MTYPE  = mtxType; // matrix type
  MAXFCT = 1;       // maximum number of factorisations of same sparsity pattern
                    //      to be kept in memory
  MNUM   = 1;       // which factorisation to use
  NRHS   = 1;       // number of right hand sides
  MSGLVL = 0;       // output message level (1 -> print statistical information)

  IPARM[0] = 0;     // PARADISO will set IPARM to default values

  IPARM[2] = max(1,numProc);  // number of processors (no default value available)

/*
  tmp = getenv("OMP_NUM_THREADS");

  if (tmp != NULL) 
  {
    sscanf(tmp,"%d", &idum);
    if (idum != IPARM[2]) prgError(1,fct,"set environment variable OMP_NUM_THREADS to numProc!");
  }
  else prgError(2,fct,"set environment variable OMP_NUM_THREADS!");
*/

  pardisoinit_(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);

  if (error != 0)
  {
    if (error == -10) cerr << "SolverPardisoPetsc::initialise ... no license file found." << endl;
    if (error == -11) cerr << "SolverPardisoPetsc::initialise ... license is expired." << endl;
    if (error == -12) cerr << "SolverPardisoPetsc::initialise ... wrong username or hostname." << endl;
  }

  cout << "\n\n PARDISO license check was successful.\n\n" << endl;

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

  PetscBool flag;

  ierr = MatGetRowIJ(mtx, 1, PETSC_FALSE, PETSC_FALSE, &nRow, &csr, &col, &flag);

  ierr = MatSeqAIJGetArray(mtx, &array);
  //ierr = MatGetArray(mtx, &array);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, &idum, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  if (error != 0)
  {
    cout << "PARDISO ERROR = " << error << endl;
    cerr << "SolverPardisoPetsc::initialise ... symbolic factorisation failed." << endl;
  }

  IPARM[5] = 1; // overwrite RHS with solution
  //IPARM[5] = 0; // do not overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  currentStatus = INIT_OK;

  return 0;
}



int SolverPardisoPetsc::factorise()
{
  if (currentStatus != ASSEMBLY_OK) 
  { 
    cerr << "SolverPARDISO::factorise ... assemble matrix first!" << endl;
    return 1;
  }

  phase = 22; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  currentStatus = FACTORISE_OK;

  return 0;
}






int  SolverPardisoPetsc::solve()
{
  if(currentStatus != FACTORISE_OK)
  {
    cerr << "SolverPARDISO::solve ... factorise matrix first!" << endl;
    return 1;
  }

  phase = 33; error = 0;

  PetscScalar *arrayRhs, *arraySoln;

  VecGetArray(rhsVec, &arrayRhs);
  VecGetArray(solnVec, &arraySoln);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, arrayRhs, arraySoln, &error, DPARM);

  //for(int ii=0; ii<nRow; ii++)
    //VecSetValue(solnVec, ii, arraySoln[ii], INSERT_VALUES);

  //VecAssemblyBegin(solnVec);
  //VecAssemblyEnd(solnVec);

  VecRestoreArray(rhsVec, &arrayRhs);
  VecRestoreArray(solnVec, &arraySoln);

  return 0;
}






int  SolverPardisoPetsc::factoriseAndSolve()
{
  char fct[] = "";

  if (currentStatus != ASSEMBLY_OK) 
  { 
    cerr << "SolverPARDISO::factoriseAndSolve ... assemble matrix first!" << endl;
    return 1;
  }

  phase = 23; error = 0;

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

  //MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);

  VecAssemblyBegin(rhsVec);
  VecAssemblyEnd(rhsVec);

  PetscScalar *arrayRhs;
  VecGetArray(rhsVec, &arrayRhs);

  int ii=0;
  for(ii=0; ii<nRow; ii++)
  {
    rhsTemp[ii] = arrayRhs[ii];
    solnTemp[ii] = 0.0;
  }

  VecRestoreArray(rhsVec, &arrayRhs);

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &rhsTemp[0], &solnTemp[0], &error, DPARM);

  VecZeroEntries(solnVec);
  for(ii=0; ii<nRow; ii++)
  {
    ierr = VecSetValue(solnVec, ii, solnTemp[ii], ADD_VALUES);
  }

  VecAssemblyBegin(solnVec);
  VecAssemblyEnd(solnVec);

  printf("Peak memory [kB] phase 1       = %d \n", IPARM[14]);
  printf("Permanent integer memory [kb]. = %d \n", IPARM[15]);
  printf("Peak real memory [kB]          = %d \n", IPARM[16]);
  printf("Number of nonzeros in LU.      = %d \n", IPARM[17]);
  printf("Gflops for LU factorization.   = %d \n", IPARM[18]);

  currentStatus = FACTORISE_OK;

  return 0;
}





int SolverPardisoPetsc::free()
{
  ierr = VecDestroy(&solnVec);CHKERRQ(ierr);
  ierr = VecDestroy(&rhsVec);CHKERRQ(ierr);
  ierr = VecDestroy(&reacVec);CHKERRQ(ierr);

  phase = -1; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, csr, col, perm, &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  PetscBool flag;

  ierr = MatRestoreRowIJ(mtx, 1, PETSC_FALSE, PETSC_FALSE, &nRow, &csr, &col, &flag);

  MatSeqAIJRestoreArray(mtx, &array);

  ierr = MatDestroy(&mtx);CHKERRQ(ierr);
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);

  currentStatus = SOLVER_EMPTY;

  return 0;
}




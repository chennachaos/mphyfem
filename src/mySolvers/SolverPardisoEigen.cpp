
#include <iostream>

#include "SolverPardisoEigen.h"
#include "FunctionsSolver.h"
#include "util.h"

using namespace std;


SolverPardisoEigen::SolverPardisoEigen()
{
  update_factorisation = false;
  return;
}




SolverPardisoEigen::~SolverPardisoEigen()
{
  phase = -1; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

}





int SolverPardisoEigen::initialise(int numProc, int matrixType, int nr)
{
  char fct[] = "SolverPARDISO::initialise";

  nRow = nCol = nr;

  soln.resize(nRow);
  soln.setZero();

  rhsVec   = soln;
  solnPrev = soln;

  perm.resize(nRow);

  //if (currentStatus != PATTERN_OK)
    //{ prgWarning(1,fct,"prepare matrix pattern first!"); return 1; }

  phase = 11; error = 0;
  int  mtxType = 11, idum;


  char *tmp;

  switch (matrixType)
  {
    case PARDISO_STRUCT_SYM: mtxType =  1; break; // real and structurally symmetric

    case PARDISO_UNSYM:      mtxType = 11; break; // real and unsymmetric

    //default:                 prgWarning(1,fct,"matrix assumed to be unsymmetric!");
  }
/*
  IPARM[3] = 31;
//  IPARM[4] = 0;  // user input permutation
//  IPARM[4] = 2;  // return the permutation
//  IPARM[5] = 1;  // overwrite RHS with solution
//  IPARM[7] = 1;  // max number of iterative refinement steps
  IPARM[9] = 4;
  IPARM[18] = 1;
  //IPARM[27] = 0; // Parallel METIS ordering
  IPARM[31] = 0; // IPARM(32) = 0 -> sparse direct solver (default) 
  //IPARM[31] = 1; // IPARM(32)=1 -> iterative method  
  //IPARM[9] = 1;
*/

  SOLVER = 0;       // sparse direct solver
  //SOLVER = 1;       // multi-recursive iterative solver
  MTYPE  = mtxType; // matrix type
  MAXFCT = 1;       // maximum number of factorisations of same sparsity pattern
                    //     to be kept in memory
  MNUM   = 1;       // which factorisation to use
  NRHS   = 1;       // number of right hand sides
  MSGLVL = 0;       // output message level (1 -> print statistical information)

  IPARM[0] = 0;     // PARADISO will set IPARM to default values
  //IPARM[0] = 1;     // user input values to IPARM
  //IPARM[1] = 2;

  tmp = getenv("OMP_NUM_THREADS");

  if (tmp != NULL) 
  {
    sscanf(tmp,"%d", &idum);
    //if (idum != IPARM[2]) prgError(1,fct,"set environment variable OMP_NUM_THREADS to numProc!");
  }
  else
    cerr << "SolverPARDISO::initialise ... set environment variable OMP_NUM_THREADS!" << endl;

  IPARM[2] = max(1,idum);  // number of processors (no default value available)
  cout <<  idum << '\t' << IPARM[2] << '\t' <<  tmp <<   endl;


  pardisoinit_(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);

  if (error != 0)
  {
    if (error == -10) cerr << "SolverPARDISO::initialise ... no license file found." << endl;
    if (error == -11) cerr << "SolverPARDISO::initialise ... license is expired." << endl;
    if (error == -12) cerr << "SolverPARDISO::initialise ... wrong username or hostname." << endl;
  }

  int  *c1, *c2, ii;

  csr.resize(nRow+1);
  col.resize(mtx.nonZeros());

  c1 = mtx.outerIndexPtr();
  c2 = mtx.innerIndexPtr();


  for(ii=0;ii<=nRow;ii++)
    csr[ii] = c1[ii] + 1;

  for(ii=0;ii<mtx.nonZeros();ii++)
    col[ii] = c2[ii] + 1;

  array = mtx.valuePtr();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  if (error != 0)
  {
    cout << "PARDISO ERROR = " << error << endl;
    cerr << "SolverPARDISO::initialise ... symbolic factorisation failed." << endl;
  }

  //IPARM[4] = 0;  // user input permutation
  //IPARM[4] = 2;  // return the permutation
  IPARM[5] = 0; // do not overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  currentStatus = INIT_OK;

  return 0;
}





int SolverPardisoEigen::factorise()
{
  if (currentStatus != ASSEMBLY_OK) 
  {
    cerr << "SolverPARDISO::factorise ... assemble matrix first!" << endl;
    return 1;
  }

  phase = 22; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  currentStatus = FACTORISE_OK;

  return 0;
}






int SolverPardisoEigen::solve()
{ 
  char fct[] = "";

  if (currentStatus != FACTORISE_OK)
  {
    cerr << "SolverPARDISO::solve ... factorise matrix first!" << endl;
    return 1;
  }

  phase = 33; error = 0;

  soln.setZero();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &soln[0], &error, DPARM);

  return 0;
}




int SolverPardisoEigen::factoriseAndSolve()
{
  if (currentStatus != ASSEMBLY_OK) 
  {
    cerr << "SolverPARDISO::factoriseAndSolve ... assemble matrix first!" << endl;
    return 1;
  }

  phase = 23; error = 0;

  soln.setZero();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &soln[0], &error, DPARM);

  if (error != 0)
  {
    cout << "PARDISO ERROR = " << error << endl;
    cerr << "SolverPARDISO::factoriseAndSolve ... factoriseAndSolve failed." << endl;
  }


//   printf("Peak memory [kB] phase 1       = %d \n", IPARM[14]);
//   printf("Permanent integer memory [kb]. = %d \n", IPARM[15]);
//   printf("Peak real memory [kB]          = %d \n", IPARM[16]);
//   printf("Number of nonzeros in LU.      = %d \n", IPARM[17]);
//   printf("Gflops for LU factorization.   = %d \n", IPARM[18]);

  currentStatus = FACTORISE_OK;

  return 0;
}



int SolverPardisoEigen::free()
{
  phase = -1; error = 0;

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &nRow, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  array = NULL;

  currentStatus = SOLVER_EMPTY;

  return 0;
}




int SolverPardisoEigen::solveArclengthSystem(int loadStep, VectorXd&  DuFull, double& Dl, double&  Ds, double& dl)
{
    double  psi = 1.0, b, A;
    double  FtF = Fext.dot(Fext);

    assert(nRow == (dispDOF+presDOF));

    VectorXd  a, du1, du2, DuShort(nRow), FextShort(nRow);
    DuShort.setZero();
    FextShort.setZero();

    for(int ii=0; ii<dispDOF; ii++)
    {
        FextShort[ii] = Fext[assyForSoln[ii]];

        DuShort[ii]   = DuFull[assyForSoln[ii]];
    }

    if(loadStep > 1)
    {
        A = DuShort.dot(DuShort) + psi*Dl*Dl*FtF - Ds*Ds;

        a = 2.0*DuShort;
        b = 2.0*psi*Dl*FtF;
    }
    else
    {
        A = 0.0;
        a = 0.0*DuShort;
        b = 1.0;
    }

    // factorise the matrix
    factorise();

    // solve the matrix system for du2
    solve();
    du2 = -1.0*soln; // minus sign is because the Residual is added to the RHS

    // solve the matrix system for du1
    rhsVec = FextShort;
    solve();
    du1 = soln;


    //cout << "aaaaaaaaa" << endl;
    //printf("A = %14.10f \t %14.10f \t %14.10f \t %14.10f \n", A, b, a.dot(du2), a.dot(du1));

    dl = (a.dot(du2) - A)/(b+a.dot(du1));

    soln = -du2 + dl*du1;

    return 0;
}




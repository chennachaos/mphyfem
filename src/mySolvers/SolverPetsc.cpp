
#include "SolverPetsc.h"
#include "petscmat.h"
#include "petscksp.h"
#include "util.h"
#include "headersEigen.h"

extern  bool debug;


SolverPetsc::SolverPetsc()
{
}


SolverPetsc::~SolverPetsc()
{
  //if (debug)  cout << " SolverPetsc destructor\n\n";
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  //cout << "SolverPetsc::~SolverPetsc() " << endl;
  //free();
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::~SolverPetsc() \n");
  //cout << "SolverPetsc::~SolverPetsc() " << endl;
}


// Initialise Petsc solver
// size_local  - number of local rows/columns in the matrix
// size_global - number of global rows/columns in the matrix
// diag_nnz - number of nonzeros in the diagonal matrix
// offdiag_nnz - number of nonzeros in the off-diagonal matrix
//
int SolverPetsc::initialise(int size_local, int size_global, int* diag_nnz, int* offdiag_nnz)
{
    nRow = nCol = size_global;

    int dummy = 50;

    // Create PETSc vector
    errpetsc = VecCreate(PETSC_COMM_WORLD, &solnVec);
    CHKERRQ(errpetsc);

    errpetsc = VecSetSizes(solnVec, size_local, size_global);
    CHKERRQ(errpetsc);

    errpetsc = VecSetFromOptions(solnVec);
    CHKERRQ(errpetsc);

    errpetsc = VecDuplicate(solnVec, &rhsVec);
    CHKERRQ(errpetsc);

    errpetsc = VecSetOption(rhsVec, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE);
    CHKERRQ(errpetsc);

    //PetscPrintf(PETSC_COMM_WORLD, " Creating PETSc matrices \n", errpetsc)

    // Create PETSc matrix
    errpetsc = MatCreate(PETSC_COMM_WORLD, &mtx);
    CHKERRQ(errpetsc);

    errpetsc = MatSetSizes(mtx, size_local, size_local, size_global, size_global);
    CHKERRQ(errpetsc);

    errpetsc = MatSetFromOptions(mtx);
    CHKERRQ(errpetsc);

    errpetsc = MatMPIAIJSetPreallocation(mtx, dummy, diag_nnz, dummy, offdiag_nnz);
    CHKERRQ(errpetsc);

    errpetsc = MatSeqAIJSetPreallocation(mtx, dummy, diag_nnz);
    CHKERRQ(errpetsc);


    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
    CHKERRQ(errpetsc);

    errpetsc = MatSetOption(mtx, MAT_NEW_NONZERO_LOCATIONS, PETSC_TRUE);
    CHKERRQ(errpetsc);

    errpetsc = MatSetOption(mtx, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);
    CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Creating KSP context ... \n\n");

    // Create the KSP context
    errpetsc = KSPCreate(PETSC_COMM_WORLD, &ksp);
    CHKERRQ(errpetsc);

    // Set the operators for the KSP context
    errpetsc = KSPSetOperators(ksp, mtx, mtx);
    CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Setting KSP context from input file ... \n\n");

    // Set KSP options from the input file
    // This is convenient as it allows to choose different options
    // from the input files instead of recompiling the code
    errpetsc = KSPSetFromOptions(ksp);    CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Creating PC context ... \n\n");

    // Get the PC context
    errpetsc = KSPGetPC(ksp, &pc);    CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n\n Setting PC context from input file ... \n\n");

    // Set PC options from the input file
    errpetsc = PCSetFromOptions(pc);    CHKERRQ(errpetsc);

    currentStatus = SOLVER_EMPTY;

    return 0;
}



int SolverPetsc::setSolverAndParameters()
{
    // set the KSP
    ///////////////////////////////////////////////

    ierr = KSPSetOperators(ksp, mtx, mtx);CHKERRQ(ierr);

    KSPSetType(ksp, KSPCG);

    ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);

    // set the PC
    ///////////////////////////////////////////////

    ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);

    ierr = KSPSetReusePreconditioner(ksp, PETSC_FALSE);

    PCSetType(pc, PCILU);

    ierr = PCSetFromOptions(pc);CHKERRQ(ierr);

    return 0;
}




int SolverPetsc::zeroMtx()
{
  //PetscPrintf(MPI_COMM_WORLD, " Matrix \n");
  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  MatZeroEntries(mtx);

  //PetscPrintf(MPI_COMM_WORLD, " rhsVec \n");
  VecAssemblyBegin(rhsVec);
  VecAssemblyEnd(rhsVec);
  VecZeroEntries(rhsVec);

  //PetscPrintf(MPI_COMM_WORLD, " reacVec \n");
  VecAssemblyBegin(reacVec);
  VecAssemblyEnd(reacVec);
  VecZeroEntries(reacVec);

  return 0;
}



int SolverPetsc::free()
{
  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");

  errpetsc = VecDestroy(&solnVec);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&rhsVec);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = VecDestroy(&reacVec);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = MatDestroy(&mtx);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = PCReset(pc);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = KSPDestroy(&ksp);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;
  errpetsc = KSPReset(ksp);CHKERRQ(errpetsc);
  //cout << " errpetsc = " << errpetsc << endl;

  //VecScatterDestroy(&ctxReac);
  //VecDestroy(&vecseqReac);

  //PetscPrintf(MPI_COMM_WORLD, "SolverPetsc::free() \n");

  return 0;
}




int SolverPetsc::printInfo()
{
  MatGetInfo(mtx, MAT_LOCAL, &info);

  PetscPrintf(MPI_COMM_WORLD, "Petsc solver:  nRow = %12d \n", nRow);
  PetscPrintf(MPI_COMM_WORLD, "               nnz  = %12d \n\n", info.nz_allocated);

  return 0;
}


int SolverPetsc::printMatrix(int dig, int dig2, bool gfrmt, int indent, bool interactive)
{
  printInfo();

  ierr = MatAssemblyBegin(mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);

  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
  //PetscViewerDrawOpen();
  //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);
  //MatView(mtx, PETSC_VIEWER_STDOUT_WORLD);
  //MatView(A,viewer);

  return 0;
}


double SolverPetsc::giveMatrixCoefficient(int row, int col)
{
  //return mtx(row,col,true);
  return 0.0;
}



int SolverPetsc::factorise()
{
  if (currentStatus != ASSEMBLY_OK) { cerr << "SolverPetsc::factorise ... assemble matrix first!" << endl; return 1; }

  if (checkIO)
  {
    // search for "nan" entries in matrix coefficients

    //if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");
  }

  currentStatus = FACTORISE_OK;

  return 0;
}


int SolverPetsc::solve()
{
  //time_t tstart, tend;
  //if(debug) PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... Matrix Assembly ...  \n\n");

  if (currentStatus != FACTORISE_OK) { cerr << "SolverPetsc::solve ... factorise matrix first!" << endl; return -1; }

  ierr = VecAssemblyBegin(rhsVec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(rhsVec); CHKERRQ(ierr);


  //VecView(rhsVec, PETSC_VIEWER_STDOUT_WORLD);

  if(debug) PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... Matrix Assembly ...  \n\n");

  ierr = MatAssemblyBegin(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mtx,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  if(debug) PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... Matrix Assembly ...  \n\n");

  //PetscViewerCreate(PETSC_COMM_WORLD, &viewer_matx);
  //PetscViewerDrawOpen();
  //PetscViewerSetFormat(viewer_matx, PETSC_VIEWER_ASCII_MATLAB);

  //MatView(mtx,PETSC_VIEWER_STDOUT_WORLD);

  if(debug) PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... rhsVec Assembly ...  \n\n");

  ierr = VecAssemblyBegin(solnVec); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(solnVec); CHKERRQ(ierr);

  if(debug) PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... soln Assembly ...  \n\n");

  //PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... vec zero ...  \n\n");

  //VecView(rhsVec, PETSC_VIEWER_STDOUT_WORLD);

  //tstart = time(0);

  //KSPSetInitialGuessNonzero(ksp, PETSC_TRUE);
  VecZeroEntries(solnVec);

  ierr = KSPSolve(ksp,rhsVec,solnVec); CHKERRQ(ierr);

  if(debug) PetscPrintf(MPI_COMM_WORLD, "  SolverPetsc::solve() ... KSP solve ...  \n\n");

  KSPConvergedReason reason;
  KSPGetConvergedReason(ksp, &reason);
  PetscInt its;
  ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);

  if(reason<0)
    PetscPrintf(MPI_COMM_WORLD, "\n Divergence in %d iterations.\n\n", its);
  else
    if(debug) PetscPrintf(MPI_COMM_WORLD, "\n Convergence in %d iterations.\n\n", its);

  //ierr = KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);//CHKERRQ(ierr);
  //VecView(solnVec, PETSC_VIEWER_STDOUT_WORLD);

  //tend = time(0);

  //cout << "SolverPetsc::solve()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

  return 0;
}






int SolverPetsc::factoriseAndSolve()
{
  if(currentStatus != ASSEMBLY_OK)
    { cerr << "SolverPetsc::factoriseAndSolve ... assemble matrix first!" << endl; return -1; }

  factorise();
  solve();

  return 0;
}




int SolverPetsc::assembleMatrixAndVector(vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj;

  int size1 = row.size();
  int size2 = col.size();

  for(ii=0;ii<size1;ii++)
  {
    VecSetValue(rhsVec, row[ii], Flocal(ii), ADD_VALUES);

    for(jj=0;jj<size2;jj++)
    {
      MatSetValue(mtx, row[ii], col[jj], Klocal(ii,jj), ADD_VALUES);
    }
  }

  return 1;
}



int SolverPetsc::assembleMatrixAndVector(int start1, int start2, vector<int>& row, vector<int>& col, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, r1, c1;

  int size1 = row.size();
  int size2 = col.size();

  for(ii=0;ii<size1;ii++)
  {
    r1 = start1 + row[ii];

    VecSetValue(rhsVec, r1, Flocal(ii), ADD_VALUES);

    for(jj=0;jj<size2;jj++)
    {
      c1 = start2 + col[jj];
      MatSetValue(mtx, r1, c1, Klocal(ii,jj), ADD_VALUES);
    }
  }

  return 0;
}




int SolverPetsc::assembleMatrixAndVectorMixedFormulation(int start, int c1, vector<int>& vec1, vector<int>& vec2, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, aa, bb, size1, size2;

  //printVector(vec1);
  //printVector(vec2);

  size1 = vec1.size();
  size2 = vec2.size();

  for(ii=0;ii<size1;ii++)
  {
    VecSetValues(rhsVec, 1, &vec1[ii], &Flocal(ii), ADD_VALUES);

    //for(jj=0;jj<size1;jj++)
      //mtx.coeffRef(vec1[ii], vec1[jj]) += Klocal(ii, jj);

    for(jj=0;jj<size2;jj++)
    {
      aa = start + vec2[jj];
      bb = size1 + jj;

      MatSetValues(mtx, 1, &vec1[ii], 1, &aa,       &(Klocal(ii,bb)), ADD_VALUES);
      MatSetValues(mtx, 1, &aa,       1, &vec1[ii], &(Klocal(bb,ii)), ADD_VALUES);
    }
  }

  for(ii=0;ii<size2;ii++)
  {
    aa = start + vec2[ii];
    VecSetValues(rhsVec, 1, &aa, &Flocal(size1+ii), ADD_VALUES);
  }
  return 0;
}



int SolverPetsc::assembleVector(int start, int c1, vector<int>& vec1, VectorXd& Flocal)
{
  int ii, jj;

  for(ii=0;ii<vec1.size();ii++)
  {
    VecSetValues(rhsVec, 1, &vec1[ii], &Flocal(ii), ADD_VALUES);
  }

  return 0;
}



int SolverPetsc::assembleMatrixAndVector2field(int var2_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, VectorXd& Flocal1, VectorXd& Flocal2, vector<int>& forAssyVar1, vector<int>& forAssyVar2)
{
  int ii, jj, row, col;

  int size1 = forAssyVar1.size();
  int size2 = forAssyVar2.size();
  //int sizeTotal = size1*size1;

  //PetscInt     row1[size1];
  //PetscScalar  array[sizeTotal];

  for(ii=0; ii<size1; ii++)
  {
    row = forAssyVar1[ii];
    if( row != -1 )
    {
      VecSetValue(rhsVec, row, Flocal1(ii), ADD_VALUES);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
          MatSetValue(mtx, row, col, K11(ii,jj), ADD_VALUES);
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          MatSetValue(mtx, row, col, K12(ii,jj), ADD_VALUES);
        }
      }
    }
  }

  for(ii=0; ii<size2; ii++)
  {
    row = forAssyVar2[ii];
    if( row != -1 )
    {
      row += var2_offset;

      VecSetValue(rhsVec, row, Flocal2(ii), ADD_VALUES);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
        {
          MatSetValue(mtx, row, col, K21(ii,jj), ADD_VALUES);
        }
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          MatSetValue(mtx, row, col, K22(ii,jj), ADD_VALUES);
        }
      }
    }
  }

  return 0;
}


int SolverPetsc::assembleMatrixAndVector3field(int var2_offset, int var3_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, MatrixXd& K13, MatrixXd& K31, MatrixXd& K33, VectorXd& Flocal1, VectorXd& Flocal2, VectorXd& Flocal3, vector<int>& forAssyVar1, vector<int>& forAssyVar2, vector<int>& forAssyVar3)
{
  int ii, jj, row, col;

  int size1 = forAssyVar1.size();
  int size2 = forAssyVar2.size();
  int size3 = forAssyVar3.size();

  //  K11, K12, K13, Flocal1
  for(ii=0; ii<size1; ii++)
  {
    row = forAssyVar1[ii];
    if( row != -1 )
    {
      VecSetValue(rhsVec, row, Flocal1(ii), ADD_VALUES);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
          MatSetValue(mtx, row, col, K11(ii,jj), ADD_VALUES);
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          MatSetValue(mtx, row, col, K12(ii,jj), ADD_VALUES);
        }
      }

      for(jj=0; jj<size3; jj++)
      {
        col = forAssyVar3[jj];
        if( col != -1 )
        {
          col += var3_offset;
          MatSetValue(mtx, row, col, K13(ii,jj), ADD_VALUES);
        }
      }
    }
  }

  // K21, K22, Flocal2
  for(ii=0; ii<size2; ii++)
  {
    row = forAssyVar2[ii];
    if( row != -1 )
    {
      row += var2_offset;

      VecSetValue(rhsVec, row, Flocal2(ii), ADD_VALUES);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
        {
          MatSetValue(mtx, row, col, K21(ii,jj), ADD_VALUES);
        }
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          MatSetValue(mtx, row, col, K22(ii,jj), ADD_VALUES);
        }
      }
    }
  }

  // K31, K33, Flocal3
  for(ii=0; ii<size3; ii++)
  {
    row = forAssyVar3[ii];
    if( row != -1 )
    {
      row += var3_offset;

      VecSetValue(rhsVec, row, Flocal3(ii), ADD_VALUES);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
        {
          MatSetValue(mtx, row, col, K31(ii,jj), ADD_VALUES);
        }
      }

      for(jj=0; jj<size3; jj++)
      {
        col = forAssyVar3[jj];
        if( col != -1 )
        {
          col += var3_offset;
          MatSetValue(mtx, row, col, K33(ii,jj), ADD_VALUES);
        }
      }
    }
  }

  return 0;
}





int SolverPetsc::assembleMatrixAndVectorSerial(vector<int>& forAssyElem, MatrixXd& Klocal, VectorXd& Flocal)
{
  int  size1 = forAssyElem.size();

  MatrixXdRM Klocal2 = Klocal;

  VecSetValues(rhsVec, size1, &forAssyElem[0], &Flocal[0], ADD_VALUES);
  MatSetValues(mtx,    size1, &forAssyElem[0], size1, &forAssyElem[0], &Klocal2(0,0), ADD_VALUES);

  return 0;
}






int SolverPetsc::assembleMatrixAndVectorParallel(vector<int>& forAssyElem, vector<int>& dof_map, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, r, kk=0;

  int  size1 = forAssyElem.size();
  int  size2 = size1*size1;

  PetscInt  row1[size1];
  PetscScalar  array[size2];

  for(ii=0;ii<size1;ii++)
  {
    r = dof_map[forAssyElem[ii]];

    VecSetValue(rhsVec, r, Flocal(ii), ADD_VALUES);

    row1[ii] = r;

    for(jj=0;jj<size1;jj++)
    {
      array[kk++] = Klocal(ii,jj);
    }
  }

  MatSetValues(mtx, size1, row1, size1, row1, array, ADD_VALUES);

  return 0;
}







int SolverPetsc::solveArclengthSystem(int loadStep, VectorXd&  DuFull, double& Dl, double&  Ds, double& dl)
{
    double  psi = 1.0, b, A;
    double  FtF = Fext.dot(Fext);

    assert(nRow == (dispDOF+presDOF));

    VectorXd  a, du1(nRow), du2(nRow), DuShort(nRow), FextShort(nRow);
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
    //du2 = -1.0*soln; // minus sign is because the Residual is added to the RHS

    PetscScalar *arrayTemp;

    VecGetArray(solnVec, &arrayTemp);

    for(int ii=0; ii<nRow; ii++)
        du2[ii] = -arrayTemp[ii];

    VecRestoreArray(solnVec, &arrayTemp);


    // solve the matrix system for du1
    //rhsVec = FextShort;

    for(int ii=0; ii<nRow; ii++)
    {
        VecSetValue(rhsVec, ii, FextShort[ii], INSERT_VALUES);
    }

    solve();
    //du1 = soln;

    VecGetArray(solnVec, &arrayTemp);

    for(int ii=0; ii<nRow; ii++)
        du1[ii] = arrayTemp[ii];

    VecRestoreArray(solnVec, &arrayTemp);


    //cout << "aaaaaaaaa" << endl;
    printf("A = %14.10f \t %14.10f \t %14.10f \t %14.10f \n", A, b, a.dot(du2), a.dot(du1));

    dl = (a.dot(du2) - A)/(b+a.dot(du1));

    //soln = -du2 + dl*du1;

    for(int ii=0; ii<nRow; ii++)
    {
        A = -du2[ii] + dl*du1[ii];

        VecSetValue(solnVec, ii, A, INSERT_VALUES);
    }

    ierr = VecAssemblyBegin(solnVec); CHKERRQ(ierr);
    ierr = VecAssemblyEnd(solnVec); CHKERRQ(ierr);

    return 0;
}










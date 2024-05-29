#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#define EIGEN_SUPERLU_SUPPORT


#include "SolverEigen.h"
#include "SolverPetsc.h"
#include "util.h"
#include <chrono>

//#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>


using namespace std;
using namespace Eigen;



SolverEigen::SolverEigen()
{
  STABILISED = false;

  update_precond = 1;
  update_factorisation = true;
  count = 0;
}


SolverEigen::~SolverEigen()
{
  free();
}


int SolverEigen::initialise(int p1, int p2, int p3)
{
   nRow = nCol = p3;

   soln.resize(nRow);
   soln.setZero();

   rhsVec   = soln;
   solnPrev = soln;
   solnTot  = soln;

  return 0;
}


int SolverEigen::setSolverAndParameters()
{
    ///////////////////////
    // Create the linear solver and set various options
    ///////////////////////


    ///////////////////////
    //Set operators. Here the matrix that defines the linear system
    //also serves as the preconditioning matrix.
    ///////////////////////

    return 0;
}


void SolverEigen::zeroMtx()
{
  mtx *= 0.0;
  rhsVec.setZero();

  return;
}



int SolverEigen::free()
{
  return 0;
}


void SolverEigen::printInfo()
{
  //cout << "Eigen solver:  nRow = " << nRow << "\n";
  //cout << "               nnz  = " <<  << "\n\n"; 
  //printVector(rhsVec);
  //rhsVec.setZero();
  //cout << " nRow = " << nRow << endl;
  //cout << mtx << endl;

  return;
}


void SolverEigen::printMatrixPatternToFile()
{
    /*
    ofstream fout("matrix-pattern2.dat");

    if(fout.fail())
    {
      cout << " Could not open the Output file" << endl;
      exit(1);
    }

    fout << totalDOF << setw(10) << totalDOF << endl;

    for(int k=0; k<mtx.outerSize(); ++k)
    for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      printf("%9d \t %9d \n", it.row(), it.col());

    //for(int ii=0;ii<nRow;ii++)
      //for(int jj=0;jj<nRow;jj++)
        //printf("%5d \t %5d \t %12.6f \n",mtx.coeffRef(ii,jj));

    fout.close();
    */

   FILE * pFile;
   int n;
   char name [100];

   cout << mtx.nonZeros() << endl;
   //pFile = fopen ("Stokes.dat","w");
   pFile = fopen ("Poisson.dat","w");

    fprintf(pFile, "%9d \t %9d \t %9d \n", nRow, nCol, 0 );

    for(int k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        //if(it.row() == it.col())
          fprintf(pFile, "%9d \t %9d \t %20.16f \n", it.row(), it.col(), it.value());
      }
    }

   fclose (pFile);

   pFile = fopen ("rhsVec.dat","w");

    for(int k=0; k<nRow; ++k)
      fprintf(pFile, "%9d \t %14.8f \n", k, rhsVec(k));

   fclose (pFile);

   pFile = fopen ("refsoln.dat","w");

    for(int k=0; k<nRow; ++k)
      fprintf(pFile, "%9d \t %14.8f \n", k, soln(k));

   fclose (pFile);
  return;
}


void SolverEigen::printMatrix(int dig, int dig2, bool gfrmt, int indent, bool interactive)
{
  printInfo();
  
  cout << mtx << endl;
  printf("\n\n");

  return;
}



int SolverEigen::factorise()
{
  if (currentStatus != ASSEMBLY_OK)
  {
    cerr << " SolverEigen::factorise ... assemble matrix first!" << endl;
    return 1;
  }

  if (checkIO)
  {
    // search for "nan" entries in matrix coefficients

    //if (prgNAN(mtx.x.x,NE)) prgError(1,fct,"nan matrix coefficient!");
  }

  currentStatus = FACTORISE_OK;

  return 0;
}


int  SolverEigen::solve()
{
  char fct[] = "SolverEigen::solve";

  time_t tstart, tend;

  if (currentStatus != FACTORISE_OK)
  {
    cerr << " SolverEigen::factorise ... factorise matrix first!" << endl;
    return 1;
  }

  //algoType = 2;

  //cout << mtx << endl;


//#ifdef EIGEN_PARALLELIZE
//#ifdef EIGEN_HAS_OPENMP
//printf("Eigen parallellize is on \n");
//#else
//printf("Eigen parallellize is off \n");
//#endif

  if(algoType == 1)
  {
    //cout << " Solving with Eigen::SimplicialLDLT " << endl;

    SimplicialLDLT<SparseMatrix<double> > solver;

    //printf("nnz =  %d \n ", mtx.nonZeros() );

    tstart = time(0);
    //SuperLU<SparseMatrixXd > solver;

    solver.compute(mtx);

    soln = solver.solve(rhsVec);

    //printVector(soln);
    //
    //computeConditionNumber();

    //VectorXd x0 = VectorXd::LinSpaced(nRow, 0.0, 1.0);

    //double  cond = myCondNum(mtx, x0, 50, solver);

    //printf("\n Matrix condition number = %12.6E \n\n\n", cond);

    //soln = solver.solveWithGuess(rhsVec, soln);
    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
    tend = time(0);
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
  }
 
  if( algoType == 2 )
  {
    cout << " Solving with Eigen::BiCGSTAB " << endl;

    //ConjugateGradient<SparseMatrixXd, Lower|Upper >   solver;
    //ConjugateGradient<SparseMatrixXd, Lower|Upper, IncompleteLUT<double> >   solver;

    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;
    //BiCGSTAB<SparseMatrixXd > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;
    //GMRES<SparseMatrixXd> solver;

    //solver.preconditioner().setDroptol(1.0e-2);
    solver.preconditioner().setFillfactor(2);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solver.setMaxIterations(2000);
    solver.setTolerance(1.0e-8);

    //tstart = time(0);
    auto  time1 = chrono::steady_clock::now();

    /*
    mtx.uncompress();
    cout << " Getting the diagonal part " << endl;
    VectorXd  diagVec1(nRow),  diagVec2(nRow),  diagVec3(nRow),  solnGuess(nRow);
    diagVec1.setZero();
    diagVec2.setZero();
    for(int ii=0; ii<nRow; ++ii)
    {
      diagVec1(ii) = mtx.coeffRef(ii,ii);
    }
    cout << " Getting the diagonal part " << endl;
    for(int ii=0; ii<(nRow-1); ++ii)
    {
      diagVec2(ii)   = mtx.coeffRef(ii,ii+1);
      diagVec3(ii+1) = mtx.coeffRef(ii+1,ii);
    }
    cout << " Solving bidiagonal part " << endl;
    solnGuess(nRow-1) = rhsVec(nRow-1)/diagVec1(nRow-1);
    for(int ii=nRow-2; ii>0; --ii)
    {
      solnGuess(ii) = (rhsVec(ii)-diagVec2(ii)*solnGuess(ii+1))/diagVec1(ii);
    }
    cout << " Solving main system " << endl;

    mtx.makeCompressed();
    */

    //cout << " iiiiiiiiiiiiiii " << endl;
    solver.compute(mtx);
    //cout << " iiiiiiiiiiiiiii " << endl;

    //printf("\n\n");
    //printVector(rhsVec);

    //soln = solver.solveWithGuess(rhsVec, soln);
    soln = solver.solve(rhsVec);

    solnTot += soln;
    //printf("\n\n");
    //printVector(soln);

    //tend = time(0); 
    //printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    auto  time2 = chrono::steady_clock::now();
    auto  duration = chrono::duration_cast<chrono::milliseconds>(time2-time1).count();

    cout << " Error       = " << solver.error() << endl;
    cout << " Iterations  = " << solver.iterations() << endl;
    cout << " Time (ms)   = " << duration << endl;
    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
  }

  if( algoType == 5 )
  {
    solverUzawaType1();
    //solverSchurGMRES();
  }

  solnPrev = soln;

  //printf("\n\n");
  //printVector(soln);

  //for(int ii=0;ii<nRow;ii++)
    //RHS[ii] = soln(ii);

  //cout << "SolverEigen::solve()  took "<< difftime(tend, tstart) <<" second(s)."<< endl;

  if(checkIO)
  {
    // search for "nan" entries in solution vector

    //if (prgNAN(RHS,N)) prgError(1,fct,"nan entry in solution vector!");
  }

  return 0;
}






int SolverEigen::factoriseAndSolve()
{
  if(currentStatus != ASSEMBLY_OK)
  {
    cerr << " assemble matrix first! " << endl;
    return 1;
  }

  factorise();

  return solve();
}




int SolverEigen::assembleMatrixAndVector(int start, int c1, vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, aa, bb, size1, r, c;

  //printVector(forAssy);

  //cout << "update_factorisation = " << update_factorisation << endl;

  size1 = forAssy.size();

  for(ii=0; ii<size1; ii++)
  {
    aa = forAssy[ii];
    if( aa != -1 )
    {
      r = start + aa;
      rhsVec[r] += Flocal(ii);

      for(jj=0; jj<size1; jj++)
      {
        bb = forAssy[jj];
        if( bb != -1 )
          mtx.coeffRef(r, start+bb) += Klocal(ii,jj);
      }
    }
  }

  return 0;
}



int SolverEigen::assembleMatrixAndVectorMixed(int var2_offset, vector<int>& forAssyVar1, vector<int>& forAssyVar2, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Fu, VectorXd& Fp)
{
  int ii, jj, row, col;

  int size1 = forAssyVar1.size();
  int size2 = forAssyVar2.size();

  for(ii=0; ii<size1; ii++)
  {
    row = forAssyVar1[ii];
    if( row != -1 )
    {
      rhsVec[row] += Fu(ii);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
          mtx.coeffRef(row, col) += Kuu(ii,jj);
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          mtx.coeffRef(row, col) += Kup(ii,jj);
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

      rhsVec[row] += Fp(ii);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
        {
          mtx.coeffRef(row, col) += Kpu(ii,jj);
        }
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          mtx.coeffRef(row, col) += Kpp(ii,jj);
        }
      }
    }
  }

  return 0;
}


int SolverEigen::assembleMatrixAndVector2field(int var2_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, VectorXd& Flocal1, VectorXd& Flocal2, vector<int>& forAssyVar1, vector<int>& forAssyVar2)
{
  int ii, jj, row, col;

  int size1 = forAssyVar1.size();
  int size2 = forAssyVar2.size();

  for(ii=0; ii<size1; ii++)
  {
    row = forAssyVar1[ii];
    if( row != -1 )
    {
      rhsVec[row] += Flocal1(ii);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
          mtx.coeffRef(row, col) += K11(ii,jj);
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          mtx.coeffRef(row, col) += K12(ii,jj);
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

      rhsVec[row] += Flocal2(ii);

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
        {
          mtx.coeffRef(row, col) += K21(ii,jj);
        }
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          mtx.coeffRef(row, col) += K22(ii,jj);
        }
      }
    }
  }

  return 0;
}



int SolverEigen::assembleMatrixAndVector3field(int var2_offset, int var3_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, MatrixXd& K13, MatrixXd& K31, MatrixXd& K33, VectorXd& Flocal1, VectorXd& Flocal2, VectorXd& Flocal3, vector<int>& forAssyVar1, vector<int>& forAssyVar2, vector<int>& forAssyVar3)
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
      rhsVec[row] += Flocal1[ii];

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
          mtx.coeffRef(row, col) += K11(ii,jj);
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          mtx.coeffRef(row, col) += K12(ii,jj);
        }
      }

      for(jj=0; jj<size3; jj++)
      {
        col = forAssyVar3[jj];
        if( col != -1 )
        {
          col += var3_offset;
          mtx.coeffRef(row, col) += K13(ii,jj);
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

      rhsVec[row] += Flocal2[ii];

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
        {
          mtx.coeffRef(row, col) += K21(ii,jj);
        }
      }

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          col += var2_offset;
          mtx.coeffRef(row, col) += K22(ii,jj);
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

      rhsVec[row] += Flocal3[ii];

      for(jj=0; jj<size1; jj++)
      {
        col = forAssyVar1[jj];
        if( col != -1 )
        {
          mtx.coeffRef(row, col) += K31(ii,jj);
        }
      }

      for(jj=0; jj<size3; jj++)
      {
        col = forAssyVar3[jj];
        if( col != -1 )
        {
          col += var3_offset;
          mtx.coeffRef(row, col) += K33(ii,jj);
        }
      }
    }
  }

  return 0;
}



int SolverEigen::assembleMatrixAndVectorUPexpl(int var2_offset, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Fu, VectorXd& Fp, vector<int>& forAssyVar1, vector<int>& forAssyVar2)
{
  int ii, jj, row, col;

  int size1 = forAssyVar1.size();
  int size2 = forAssyVar2.size();

  for(ii=0; ii<size1; ii++)
  {
    row = forAssyVar1[ii];
    if( row != -1 )
    {
      rhsVec[row] += Fu(ii);

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          matKup.coeffRef(row, col) += Kup(ii,jj);
          matKpu.coeffRef(col, row) += Kpu(jj,ii);
        }
      }
    }
  }

  for(ii=0; ii<size2; ii++)
  {
    row = forAssyVar2[ii];
    if( row != -1 )
    {
      rhsVec2[row] += Fp(ii);

      for(jj=0; jj<size2; jj++)
      {
        col = forAssyVar2[jj];
        if( col != -1 )
        {
          matKpp.coeffRef(row, col) += Kpp(ii,jj);
        }
      }
    }
  }

  return 0;
}


int SolverEigen::assembleMatrixAndVectorElecField(vector<int>& forAssy, MatrixXd& Klocal, VectorXd& Flocal)
{
  int ii, jj, r, c, size1 = forAssy.size();

  for(ii=0; ii<size1; ii++)
  {
    r = forAssy[ii];
    if( r != -1 )
    {
      rhsVec[r] += Flocal(ii);

      for(jj=0; jj<size1; jj++)
      {
        c = forAssy[jj];
        //cout <<  ii << '\t' << r <<  '\t' << jj << '\t' <<  c <<  endl;
        if( c != -1 )
          mtx.coeffRef(r, c) += Klocal(ii,jj);
      }
    }
  }

  return 0;
}




int SolverEigen::solveArclengthSystem(int loadStep, VectorXd&  DuFull, double& Dl, double&  Ds, double& dl)
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


    cout << "aaaaaaaaa" << endl;
    printf("A = %14.10f \t %14.10f \t %14.10f \t %14.10f \n", A, b, a.dot(du2), a.dot(du1));

    dl = (a.dot(du2) - A)/(b+a.dot(du1));

    soln = -du2 + dl*du1;

    return 0;
}



int SolverEigen::solveArclengthSystemGrowthModel(int loadStep, VectorXd&  DuFull, double& Dl, double&  Ds, double& dl)
{
    double  psi = 1.0, b, A;
    double  FtF = Fext.dot(Fext);
    VectorXd  a, du1, du2, DuShort(nRow), FintDerGrowthShort(nRow);

    for(int ii=0; ii<dispDOF; ii++)
    {
        FintDerGrowthShort[ii] = FintDerGrowth[assyForSoln[ii]];

        DuShort[ii]            = DuFull[assyForSoln[ii]];
    }

    for(int ii=0; ii<presDOF; ii++)
    {
        FintDerGrowthShort[ii] = FintDerGrowth[assyForSoln[ii]];
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
    rhsVec = FintDerGrowthShort;
    solve();
    du1 = soln;

    dl = (a.dot(du2) - A)/(b+a.dot(du1));

    soln = -du2 + dl*du1;

    return 0;
}




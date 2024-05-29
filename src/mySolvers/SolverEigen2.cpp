
#include "SolverEigen.h"
#include "util.h"

//#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>


using namespace std;
using namespace Eigen;


void SolverEigen::computeConditionNumber()
{
/*
    MatrixXd  globalK;
    
    globalK.resize(nRow, nCol);
    globalK.setZero();

    int k, ii, jj;

    for(k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();
        
        //cout << ii << '\t' << jj << '\t' << it.value() << endl;

        globalK.coeffRef(ii, jj) = it.value();
      }
    }


    VectorXd sing_vals = globalK.jacobiSvd().singularValues();

    //printf("\n Matrix condition number = %12.6E \n", sing_vals(0)/sing_vals(sing_vals.size()-1) );
    
    //printf("\n Minimum eigenvalue = %12.6f \n", sing_vals.minCoeff() );
    //printf("\n Minimum eigenvalue = %12.6f \n", sing_vals.maxCoeff() );
    printf("\n Matrix condition number = %12.6E \n", sing_vals.maxCoeff() / sing_vals.minCoeff() );
    printf("\n\n\n\n");
*/

    //myMatrix::OneNormEst estimator(mtx.rows(), 4);
    //double norm_a;
    
    //estimator.ANorm(mtx, mtx, norm_a);

    //cout << " norm_a = " << norm_a << endl;
//

  //VectorXd x0 = VectorXd::LinSpaced(nRow, 0.0, 1.0);

  //SuperLU<SparseMatrixXd>  solver;

  //double  cond = myCondNum(mtx, x0, 50, solver);

  //printf("\n Matrix condition number = %12.6E \n", cond);

  //myCondNumMatlab(mtx);

  return;
}





int  SolverEigen::solverSchurExplicitDirect(int frequency)
{
    //cout << " SolverEigen::solverSchurExplicitDirect() " << endl;
    //  Kuu*x1 + Kup*x2 = f1
    //  Kpu*x1 + Kpp*x2 = f2
    //
    //  Kuu*x1 = f1 - Kup*x2
    //  S*x2 = f2 - Kpu*KuuInv*f1
    //
    //  S = Kpp - Kpu*KuuInv*Kup,  Schur complement

    //cout << matKup << endl;
    //cout << matKpu << endl;

    // compute the Schur complement
    //if( (count % frequency) == 0 )
    //{
      //cout <<  "Updating the Schur complement solver " <<  endl;
      matSchur = matKpp - matKpu*(matKuuInv*matKup);
      matSchur.makeCompressed();

      solverSchur.compute(matSchur);
    //}

    //cout << matS << endl;

    //cout << " matKuuInv ... " << matKuuInv.rows() << '\t' << matKuuInv.cols() << endl;
    //cout << " matKup    ... " << matKup.rows() << '\t' << matKup.cols() << endl;
    //cout << " matKpu    ... " << matKpu.rows() << '\t' << matKpu.cols() << endl;
    //cout << " matKpp    ... " << matKpp.rows() << '\t' << matKpp.cols() << endl;
    //cout << " rhsVec    ... " << rhsVec.rows() << endl;
    //cout << " rhsVec2   ... " << rhsVec2.rows() << endl;

    //solver.preconditioner().setDroptol(1.0e-2);
    //solver.preconditioner().setFillfactor(1);

    //auto  time1 = chrono::steady_clock::now();

    //printVector(rhsVec);
    //printVector(rhsVec2);

    VectorXd  r2 = rhsVec2 - matKpu*(matKuuInv*rhsVec);
    //printVector(r2);
    var2 = solverSchur.solve(r2);
    //var2 = solverSchur.solveWithGuess(r2, var2);
    //printVector(var2);

    var1 = matKuuInv*(rhsVec - matKup*var2);
    //printVector(var1);

    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

    count++;

    return 0;
}






void SolverEigen::setupMatricesAndVectors()
{
    int nU, nP, ii, jj, niter, k;
    double alpha, TOL, rho, beta, fact, f3norm;
    bool preCond;

    nU = 1808;
    nP = 255;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    k = mtx.nonZeros();
    matKuu.resize(nU,nU);    matKuu.reserve(k*2/3);
    matKup.resize(nU,nP);    matKup.reserve(k/3);
    matKpu.resize(nP,nU);    matKpu.reserve(k/3);
    matKpp.resize(nP,nP);    matKpu.reserve(k/3);

    PreCondSchur.resize(nP, nP);
    PreCondSchur.reserve(nP*nP);

    for(k=0; k<mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(mtx,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();

        if( ii < nU )
        {
          if( jj < nU )
            matKuu.coeffRef(ii, jj) = it.value();
          else
            matKup.coeffRef(ii, jj-nU) = it.value();
        }
        else
        {
          if( jj < nU )
            matKpu.coeffRef(ii-nU, jj) = it.value();
          else
            matKpp.coeffRef(ii-nU, jj-nU) = it.value();
        }
      }
    }

    //for(ii=0;ii<nU;ii++)
      //tempMat.coeffRef(ii,ii) = 1.0/matKuu.coeffRef(ii,ii);

    //PreCondSchur = matKup*tempMat*matKup.transpose();
    //PreCondSchur = matKup*matKup.transpose();

    /*
    for(ii=0;ii<nP;ii++)
      PreCondSchur.coeffRef(ii,ii) = 0.0;

    for(k=0; k<matKup.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(matKup, k); it; ++it)
      {
        ii   = it.row();
        jj   = it.col();
        fact = it.value();

        PreCondSchur.coeffRef(ii,ii) += (fact*fact/matKuu.coeffRef(jj,jj));
      }
    }

    for(k=0; k<PreCondSchur.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(PreCondSchur, k); it; ++it)
        cout << it.row() << '\t' << it.col() << '\t' << it.value() << endl;
    }
    */

    f1.resize(nU);
    f2.resize(nP);
    f1.setZero();
    f2.setZero();

    for(ii=0;ii<nU;ii++)
      f1(ii) = rhsVec(ii);

    for(ii=0;ii<nP;ii++)
      f2(ii) = rhsVec(nU+ii);

    var1.resize(nU);
    var2.resize(nP);
    var1.setZero();
    var2.setZero();

    var1Prev = var1;
    var2Prev = var2;

    return;
}




void SolverEigen::solverUzawaType1()
{

    return;
}



void SolverEigen::solverUzawaType2()
{
    setupMatricesAndVectors();

    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k;
    double alpha, TOL, rho, beta, normPrev, error;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    
    //nU = 10918;
    //nP = 5459;
    //nL = 404;
    //nL = 0;

    VectorXd  tempVec, f4, r2, r2prev, z, q;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matKuu);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-10);

    char        tmp[200];
    //MyString    tmpStr;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    cout << f1.norm() << '\t' << f2.norm() << '\t' << f3.norm() << endl;

    preCond = false;
    //preCond = true;
    TOL = 1.0e-5;
    niter = 200;
    alpha = 30.0;
    beta  = 0.001;

    tstart = time(0);
  
/*
    r2 = f1 - matKup.transpose()*var2 - matKpu.transpose()*var3;
    normPrev =  f2.norm() + 100.0;

    for(ii=0;ii<niter;ii++)
    {
      printf("%5d \t %12.6f \t %12.6f \t %12.6f \n", ii, var1.norm(), var2.norm(), var3.norm());

      var1 = solver.solveWithGuess(r2, var1Prev);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
      //var2 += alpha * solverPCSchur.solve(matKup*var1 - f2);
      var2 += alpha * (matKup*var1 - f2);
      var3 += beta  * (matKpu*var1 - f3);

      //r2 = f1 - matKup.transpose()*var2 - matKpu.transpose()*var3;
      
      error = r2.norm();

      tempVec = var1 - var1Prev;
      //printf("%5d \t %12.6f \t %12.6f \n", k, var1.norm(), var2.norm());

      tmpStr.free();
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, tempVec.norm());
      sprintf(tmp," \t %6d \t %12.6E \t %12.6E \n", ii, error, tempVec.norm());
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      if( error < TOL )
        break;

      var1Prev = var1;
    }
*/
//
    for(ii=0;ii<niter;ii++)
    {
     // //var1 = solver.solve(f1 - matKup.transpose()*var2);
      //var1 = solver.solveWithGuess(f1 - matKup.transpose()*var2, var1Prev);
      cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
      ////var2 += alpha * (matKup*var1 - matKpp*var2);
      //var2 += alpha * (matKup*var1 - f2);
      ////var2 += alpha * solverPCSchur.solve(matKup*var1 - f2);

      tempVec = var1 - var1Prev;
      //printf("%5d \t %12.6f \t %12.6f \n", k, var1.norm(), var2.norm());

      //tmpStr.free();
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, tempVec.norm());
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, var1.norm());
      //tmpStr.append(tmp);

      //prgWriteToTFile(tmpStr);

      //if( tempVec.norm() < TOL )
        //break;

      var1Prev = var1;
      //printf("%5d \t %12.6f \t %12.6f \n", ii, u.norm(), p.norm());
    }
//

    VectorXd  var1b(nU), var1c(nU), var2b(nP), var2c(nP);
    var1b.setZero();
    var1c.setZero();

    var2b.setZero();
    var2c.setZero();

/*
    // with Aitken acceleration

    var1 = solver.solve(f1 - matKup.transpose()*var2);
    var2 += alpha * (matKup*var1 - f2);

    for(ii=0;ii<niter;ii++)
    {
      var1b = solver.solve(f1 - matKup.transpose()*var2);
      //var1b = solver.solveWithGuess(f1 - matKup.transpose()*var2, var1Prev);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
      //var2 += alpha * (matKup*var1 - matKpp*var2);
      var2b += alpha * (matKup*var1b - f2);
      //var2 += alpha * solverPCSchur.solve(matKup*var1 - f2);

      //var1c = solver.solveWithGuess(f1 - matKup.transpose()*var2b, var1Prev);
      var1c = solver.solve(f1 - matKup.transpose()*var2b);
      var2c += alpha * (matKup*var1c - f2);

      for(k=0;k<nU;k++)
        var1(k) = var1(k) - (var1b(k)-var1(k))*(var1b(k)-var1(k))/(var1c(k)-2.0*var1b(k)+var1(k));
      //var1 = var1c;

      for(k=0;k<nP;k++)
        var2(k) = var2(k) - (var2b(k)-var2(k))*(var2b(k)-var2(k))/(var2c(k)-2.0*var2b(k)+var2(k));

      //tempVec = var1 - var1Prev;
      //printf("%5d \t %12.6f \t %12.6f \n", k, var1.norm(), var2.norm());

      tmpStr.free();
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, tempVec.norm());
      sprintf(tmp," \t %6d \t %12.6E \n", ii, var1c.norm());
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      //if( tempVec.norm() < TOL )
        //break;

      var1Prev = var1;
      //printf("%5d \t %12.6f \t %12.6f \n", ii, u.norm(), p.norm());
    }
*/
    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
/*
    TOL = 1.0e-8;
    niter = 1000;
    alpha = 20.0;

    f4 = matKup*solver.solve(f1) - f2;
    
    r2 = f4 - matKup*solver.solve(matKup.transpose()*var2);
    r2prev = r2;
    
    normPrev = r2.norm() + 10.0;

    for(ii=0;ii<niter;ii++)
    {
      //z = solverPCSchur.solve(r2);
      z = r2;

      var2 += alpha * z;

      r2 = f4 - matKup*solver.solve(matKup.transpose()*var2);
//
      if(r2.norm() < normPrev)
      {
        //r2prev = r2;
        alpha *= 1.1;
        normPrev = r2.norm();
      }
      else
      {
        r2 = r2prev;
        alpha *= 0.5;
      }
      //printf(" \t %6d \t %12.6E \n", ii, alpha);
//
//
      if(r2.norm() > 1.0e-4)
        var2 += alpha * r2;
      else
      {
        var2b = var2  + alpha * r2;
        var2c = var2b + alpha * (f4 - matKup*solver.solve(matKup.transpose()*var2b));

        for(k=0;k<nP;k++)
          var2(k) = var2(k) - (var2b(k)-var2(k))*(var2b(k)-var2(k))/(var2c(k)-2.0*var2b(k)+var2(k));
      }
//

      //printf("%5d \t %12.6f \t %12.6f \n", ii, var2b.norm(), var2c.norm());

      tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \t %12.6E \n", ii, r2.norm(), alpha);
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, r2.norm());
      //sprintf(tmp," \t %6d \t %12.6E \n", ii, var2.norm());
      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

      //if( tempVec.norm() < TOL )
        //break;
    }

    var1 = solver.solve(f1 - matKup.transpose()*var2);
    cout << " Number of iterations = " << (ii+1) << endl;

*/

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);

    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

    /*
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver2;
    GMRES<SparseMatrixXd, IncompleteLUT<double> > solver2;

    solver2.compute(mtx);
    solver2.setMaxIterations(100);
    solver2.setTolerance(1.0e-16);

    //solver2.preconditioner().setDroptol(1.0e-15);
    solver2.preconditioner().setFillfactor(4);

    soln = solver2.solveWithGuess(rhsVec, soln);
    */

  return;
}







void SolverEigen::solverSchurCG()
{
    setupMatricesAndVectors();
    
    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k;
    double alpha, TOL, rho, rhoPrev, beta, fact;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    nL = 0;
  
    VectorXd  tempVec, f4, r, Ap, z, p;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matKuu);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-16);

    char        tmp[200];
    //MyString    tmpStr;

    preCond = false;
    preCond = true;
    TOL = 1.0e-8;
    niter = 300;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    tstart = time(0);

    //f4 = matKup*solver.solve(f1) - f2;

    normRef = f4.norm();

    cout << f1.norm() << '\t' << f2.norm() << '\t' << normRef << endl;

    //r = f4 - matKup*solver.solve(matKup.transpose()*var2);

    if(preCond)
      z = solverPCSchur.solve(r);
    else
      z = r;
    
    p = z;

    for(ii=0;ii<niter;ii++)
    {
      rho = r.dot(z);

      //Ap = matKup*solver.solve(matKup.transpose()*p);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

      alpha = rho / (p.dot(Ap) );
      
      var2  += alpha * p; // update the solution variable (var2 here)
      r  -= alpha * Ap;   // update the residual

      fact = r.norm()/normRef;
      
      //tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \t %12.6E \t %12.6E \t %12.6E \t %12.6E \n", ii, fact, p.norm(), rho, alpha, beta);
      //tmpStr.append(tmp);

      //prgWriteToTFile(tmpStr);

      //printf("%5d \t %12.6f \n", ii, var2.norm());

      if( fact < TOL )
        break;
      
      if(preCond)
        z = solverPCSchur.solve(r);
      else
        z = r;

      beta = z.dot(r)/rho;
      p = z + beta*p;
    }

    cout << " Number of iterations = " << (ii+1) << endl;

    var1 = solver.solve(f1 - matKup.transpose()*var2);
    cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
    
    var1Prev = var1;
    var2Prev = var2;
    var3Prev = var3;

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);
  
    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

  return;
}


void SolverEigen::solverSchurBiCGSTAB()
{
    setupMatricesAndVectors();

    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k;
    double alpha, TOL, rho, beta, fact, omega, rhoPrev, error;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    nL = 0;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    VectorXd  tempVec, f4, r, rtilde, v, t, s, shat, z, p, phat;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matKuu);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-16);

    char        tmp[200];
    //MyString    tmpStr;

    preCond = false;
    //preCond = true;
    TOL = 1.0e-6;
    niter = 35;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    p.resize(nP);
    p.setZero();
    v = p;
    t = p;

    tstart = time(0);

    ////f4 = matKup*solver.solve(f1) - f2;

   // r = f4 - matKup*solver.solve(matKup.transpose()*var2);

    rho = rhoPrev = alpha = beta = omega = 1.0;

    normRef = f4.norm();

    rtilde = r;

    cout << f1.norm() << '\t' << f2.norm() << '\t' << f4.norm() << endl;

    for(ii=0;ii<niter;ii++)
    {
      rho = rtilde.dot(r);
      
      if(CompareDoubles(rho, 0.0))
        break;

      beta = (rho/rhoPrev)/(alpha/omega);
      p = r + beta*(p-omega*v);

      if(preCond)
        phat = solverPCSchur.solve(p);
      else
        phat = p;

      //v = matKup*solver.solve(matKup.transpose()*phat);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

      alpha = rho / (rtilde.dot(v) );
      
      s = r - alpha*v;

      if(s.norm() < TOL)
      {
        var2 = var2 + alpha*phat;
        break;
      }

      if(preCond)
        shat = solverPCSchur.solve(s);
      else
        shat = s;

      //t = matKup*solver.solve(matKup.transpose()*shat);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

      //printf("%5d \t %12.6f \t %12.6f \t %12.6f \t %12.6f \n", ii, s.norm(), t.norm(), p.norm(), var2.norm());
      printf("%5d \t %12.6E \t %12.6f \t %12.6f \t %12.6f \n", ii, rho, alpha, beta, omega);

      omega = t.dot(s)/t.dot(t);

      if( omega == 0.0 )
        break;

      //cout << rho << '\t' << alpha << '\t' << omega << endl;

      var2 = var2 + alpha*phat + omega*shat;
      r = s - omega*t;
      error = r.norm()/normRef;

      //tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \n", ii, error);
      //tmpStr.append(tmp);

      //prgWriteToTFile(tmpStr);
      
      if( error < TOL )
        break;

      rhoPrev = rho;
    }

    int flag;

    if( error <= TOL || s.norm() <= TOL )  // converged
    {
      flag =  0;
      if ( s.norm() <= TOL )
        error = s.norm() / normRef;
    }
    else if( omega == 0.0 ) // breakdown
      flag = -2;
    else if( rho == 0.0 )
      flag = -1;
    else        // no convergence
      flag = 1;

    cout << " Number of iterations = " << (ii+1) << endl;
    cout << " flag " << flag << endl;

    var1 = solver.solve(f1 - matKup.transpose()*var2);
    //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
    
    var1Prev = var1;
    var2Prev = var2;
    var3Prev = var3;

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);
  
    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

  return;
}


void SolverEigen::solverSchurGMRES()
{
    setupMatricesAndVectors();

    time_t tstart, tend;

    int nU, nP, nL, ii, jj, niter, k, m;
    double alpha, TOL, rho, beta, fact;
    bool preCond;

    nP = nRow/3;
    nU = nP*2;
    nL = 0;

    VectorXd  tempVec, f4, r, r1, r2, Ap, Ap1, Ap2, z, z1, z2, p, p1, p2;

    cout << " nU and nP " << nU << '\t' << nP << endl;

    //BiCGSTAB<SparseMatrixXd> solver;
  
    BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //GMRES<SparseMatrixXd, IncompleteLUT<double> > solver;

    //ConjugateGradient<SparseMatrixXd>  solver;

    solver.compute(matKuu);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-16);

    GMRES<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solverPCSchur;
    //ConjugateGradient<SparseMatrixXd> solverPCSchur;
  
    //solver.preconditioner().setDroptol(1.0e-15);
    //solver.preconditioner().setFillfactor(3);

    //solver.set_restart(50);
    //solver.setEigenv(10);

    solverPCSchur.compute(PreCondSchur);

    solverPCSchur.setMaxIterations(1000);
    solverPCSchur.setTolerance(1.0e-16);

    char        tmp[200];
    //MyString    tmpStr;

    preCond = false;
    //preCond = true;
    TOL = 1.0e-7;
    niter = 400;

    //for(ii=0;ii<nU;ii++)
      //var1(ii) = solnPrev(ii) ;

    //for(ii=0;ii<nP;ii++)
      //var2(ii) = solnPrev(nU+ii);

    tstart = time(0);

    ////f4 = matKup*solver.solve(f1) - f2;

    //r = f4 - matKup*solver.solve(matKup.transpose()*var2);
    
    rho = alpha = beta = 1.0;

    p.resize(nP);
    p.setZero();

    normRef = f4.norm();

    cout << f1.norm() << '\t' << f2.norm() << '\t' << f4.norm() << endl;
    
    m = 30;
    VectorXd  s(m+1), cs(m+1), sn(m+1);
    double  w;
    MatrixXd  v(nP, m+1);

    for(ii=0;ii<niter;ii++)
    {
      Ap1 = matKuu*p1 + matKup.transpose()*p2;
      Ap2 = matKup*p1;

      rho = r1.dot(z1) + r2.dot(z2);

      alpha = rho / (p1.dot(Ap1)+p2.dot(Ap2));

      var1  += alpha * p1;
      var2  += alpha * p2; // update the solution variable (var2 here)

      r1 -= alpha * Ap1;   // update the residual
      r2 -= alpha * Ap2;

      fact = var1.norm();

      //tmpStr.free();
      sprintf(tmp," \t %6d \t %12.6E \n", ii, fact);
      //tmpStr.append(tmp);

      //prgWriteToTFile(tmpStr);

      //printf("%5d \t %12.6f \n", k, fact);

      if( fact < TOL )
        break;

      z1 = r1;
      //z2 = -matKup*solver.solve(r1) + r2;

      beta = (z1.dot(r1)+z2.dot(r2))/rho;

      p1 = z1 + beta*p1;
      p2 = z2 + beta*p2;
    }

    cout << " Number of iterations = " << (ii+1) << endl;

    tend = time(0); 
    printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );
    
    var1Prev = var1;
    var2Prev = var2;
    var3Prev = var3;

    for(ii=0;ii<nU;ii++)
      soln(ii) = var1(ii);
  
    for(ii=0;ii<nP;ii++)
      soln(nU+ii) = var2(ii);

    for(ii=0;ii<nL;ii++)
      soln(nU+nP+ii) = var3(ii);

  return;
}



void  SolverEigen::updatePreconditioner()
{
  ////////////////////// 
  //
  //  preconditioner
  //
  ////////////////////// 
 
  //precond1.setDroptol(1.0e-3);
  //precond1.setFillfactor(2);
  //precond1.compute(mtx);

    time_t tstart, tend;
    tstart = time(0);

/*
  precond.free();

  precond.setDroptol(1.0e-3);
  precond.setFillfactor(3);

  precond.compute(mtx);
*/
  tend = time(0); 
  printf("ILUT preconditioner took %8.4f second(s) \n ", difftime(tend, tstart) );

  //mtx.makeCompressed();

  //DiagonalPreconditioner<double>  precond(mtx);

  update_precond++;

  return;

}



int  SolverEigen::solverSchurExplicitCG()
{
    int nU = matA.rows();
    int nP = matB.rows();

    int  ii, jj, niter=500, k;
    double alpha, TOL= 1.0e-8, rho, rhoPrev, beta, fact;
    bool preCond=true;
    //bool preCond=false;

    //MatrixXd  BAinvBt = matB*(matAinv*matB.transpose());
    //MatrixXd  BAinvBt_diag_inv(nP, nP);

    //BAinvBt_diag_inv.setZero();
    //for(ii=0; ii<nP; ii++)
      //BAinvBt_diag_inv(ii,ii) = 1.0/BAinvBt(ii,ii);

    VectorXd  r1, r2, Ap, z, p, f4, r;

    f4 = matB*(matAinv*rhsVec) - rhsVec2;

    double  normRef = f4.norm();

    //cout << f1.norm() << '\t' << f2.norm() << '\t' << normRef << endl;

    r = f4 - matB*(matAinv*(matB.transpose()*var2));

    //if(preCond)
      //z = BAinvBt_diag_inv*r;
    //else
      z = r;
    
    p = z;

    for(ii=0;ii<niter;ii++)
    {
      rho = r.dot(z);

      Ap = matB*(matAinv*(matB.transpose()*p));
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;

      alpha = rho / (p.dot(Ap) );

      var2  += alpha * p; // update the solution variable (var2 here)
      r   -= alpha * Ap;   // update the residual

      fact = r.norm()/normRef;
      
      //printf("%6d \t %12.6E \t %12.6E \t %12.6E \t %12.6E \t %12.6E \n", ii, fact, p.norm(), rho, alpha, beta);

      if( fact < TOL )
        break;

      //if(preCond)
        //z = BAinvBt_diag_inv*r;
      //else
        z = r;

      beta = z.dot(r)/rho;
      p = z + beta*p;
    }
    //cout << " Number of iterations = " << (ii+1) << endl;

    var1 = matAinv*(rhsVec - matB.transpose()*var2);

    return 0;
}






int  SolverEigen::solverSchurExplicitDirect()
{
    //cout << " SolverEigen::solverSchurExplicitDirect() " << endl;
    //  A*x1 + B*x2 = f1
    //  C*x1 + D*x2 = f2
    //
    //  A*x1 = f1 - B*x2
    //  S*x2 = f2 - C*Ainv*f1
    //
    //  S = D - C*Ainv*B, Schur complement

    //cout << matB << endl;
    //cout << matC << endl;

  //if(update_precond)
  //{
    // compute the Schur complement
    matS = matD - matC*(matAinv*matB);
    matS.makeCompressed();

    /*
    ofstream fout("matrix-pattern2-Kpp.dat");

    if(fout.fail())
    {
      cout << " Could not open the Output file" << endl;
      exit(1);
    }

    fout << matD.rows() << '\t' << matD.rows() << '\t' <<  0 <<  endl;

    for(int k=0; k<matD.outerSize(); ++k)
    for(SparseMatrixXd::InnerIterator it(matD,k); it; ++it)
      fout << it.row() <<  '\t' <<  it.col() << '\t' <<  it.value() <<   endl;

    fout.close();

    ofstream fout2("matrix-pattern2-Schur.dat");

    if(fout2.fail())
    {
      cout << " Could not open the Output file" << endl;
      exit(1);
    }

    fout2 << matS.rows() << '\t' << matS.rows() << '\t' <<  0 <<  endl;

    for(int k=0; k<matS.outerSize(); ++k)
    for(SparseMatrixXd::InnerIterator it(matS,k); it; ++it)
      fout2 << it.row() <<  '\t' <<  it.col() << '\t' <<  it.value() <<   endl;

    fout2.close();
    */

    //cout << matS << endl;

    //cout << " matAinv ... " << matAinv.rows() << '\t' << matAinv.cols() << endl;
    //cout << " matB    ... " << matB.rows() << '\t' << matB.cols() << endl;
    //cout << " matC    ... " << matC.rows() << '\t' << matC.cols() << endl;
    //cout << " matD    ... " << matD.rows() << '\t' << matD.cols() << endl;
    //cout << " rhsVec  ... " << rhsVec.rows() << endl;
    //cout << " rhsVec2 ... " << rhsVec2.rows() << endl;

    //SimplicialLDLT<SparseMatrix<double> > solver;
    //SuperLU<SparseMatrixXd > solver;

    //ConjugateGradient<SparseMatrixXd, Lower|Upper >   solver;
    //ConjugateGradient<SparseMatrixXd, Lower|Upper, IncompleteLUT<double> >   solver;

    BiCGSTAB<SparseMatrixXd> solver;
    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //solver.preconditioner().setDroptol(1.0e-3);
    //solver.preconditioner().setFillfactor(2);

    solver.setMaxIterations(1000);
    solver.setTolerance(1.0e-8);

    //auto  time1 = chrono::steady_clock::now();

    //solverSLU.compute(matS);
    solver.compute(matS);

    update_precond = 0;
  //}

    //printVector(rhsVec);
    //printVector(rhsVec2);

    //printVector(rhsVec);
    //cout << " iiiiiiiiiii " << endl;
    VectorXd  r2 = rhsVec2 - matC*(matAinv*rhsVec);
    //printVector(r2);
    //cout << " iiiiiiiiiii " << endl;
    //var2 = solver.solveWithGuess(r2, var2);
    //var2 = solverSLU.solve(r2);
    var2 = solver.solve(r2);

    //printVector(var2);

    //cout << " Error       = " << solver.error() << endl;
    //cout << " Iterations  = " << solver.iterations() << endl;
    //cout << " Time (ms)   = " << duration << endl;

    var1 = matAinv*(rhsVec - matB*var2);

    //printVector(var1);

    return 0;
}










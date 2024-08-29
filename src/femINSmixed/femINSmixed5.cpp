
#include "femINSmixed.h"
#include "util.h"
#include "ElementBase.h"
#include <chrono>
#include "KimMoinFlow.h"
#include "BasisFunctionsBernstein.h"
#include "MyTime.h"
#include "TimeFunction.h"

extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime                 myTime;
extern  bool  debug;




int femINSmixed::calcMassMatrixForExplicitDynamics()
{
    cout << " femINSmixed::calcMassMatrixForExplicitDynamics() STARTED " << endl;

    int  dd, ee, ii, jj;
    vector<int> nodeNums, nodeNumsPres, globalDOFnums;

    VectorXd  Mlocal1(200), Mlocal2(30);

    // Compute global mass matrices
    //The Mass is assumed to be lumped so that the mass matrix is diagonal

    globalMassVelo.resize(nNode_Velo*ndim);
    globalMassVelo.setZero();

    globalMassPres.resize(nNode_Velo);
    globalMassPres.setZero();

    Mlocal1.setZero();
    Mlocal2.setZero();

    for(ee=0; ee<nElem_global; ++ee)
    {
        // compute mass matrix and assemble it
        elems[ee]->MassMatrices(nodeCoords, fluidProperties, Mlocal1, Mlocal2);

        nodeNums = elems[ee]->nodeNums;
        globalDOFnums = elems[ee]->globalDOFnums;

        //Assemble the element mass to the global mass
        for(ii=0; ii<globalDOFnums.size(); ++ii)
        {
            globalMassVelo(globalDOFnums[ii])   +=  Mlocal1(ii);
        }

        nodeNumsPres = elems[ee]->nodeNumsPres;

        for(ii=0; ii<npElemPres; ++ii)
        {
            globalMassPres(nodeNumsPres[ii]) += Mlocal2(ii);
        }
    }

    //printVector(globalMassVelo);
    //printVector(globalMassPres);

    cout << " femINSmixed::calcMassMatrixForExplicitDynamics() FINISHED " << endl;

    return 0;
}




int femINSmixed::setSolverDataForSemiImplicit()
{
    cout << " setSolverDataForSemiImplicit() ... STARTED " << endl;


    int  ee, ii, jj, size1, size2, row, col;
    vector<int>  vecIntTemp(10);

    // inverse of Kuu matrix
    matMuuInv.setZero();
    matMuuInv.resize(totalDOF_Velo, totalDOF_Velo);
    matMuuInv.reserve(totalDOF_Velo);

    for(ii=0; ii<totalDOF_Velo; ii++)
    {
      matMuuInv.coeffRef(ii,ii) = 0.0;
    }

    VectorXi  nnzVec(totalDOF_Velo);

    // Kup matrix
    matKup.setZero();
    matKup.resize(totalDOF_Velo, totalDOF_Pres+totalDOF_Lambda);
    matKup.reserve(ceil(totalDOF_Pres*totalDOF_Velo*0.1));
    jj = ceil(totalDOF_Pres*0.12);
    cout << " jj = " << jj << endl;

    for(ii=0; ii<totalDOF_Velo; ii++)
      nnzVec(ii) = jj;
    matKup.reserve(nnzVec);

    // Kpu matrix
    matKpu.setZero();
    matKpu.resize(totalDOF_Pres+totalDOF_Lambda, totalDOF_Velo);
    matKpu.reserve(ceil(totalDOF_Pres*totalDOF_Velo*0.1));
    matKpu.reserve(nnzVec);

    // Schur complement matrix
    matSchur.setZero();
    matSchur.resize(totalDOF_Pres+totalDOF_Lambda, totalDOF_Pres+totalDOF_Lambda);
    matSchur.reserve(ceil(totalDOF_Pres*totalDOF_Pres*0.2));
    matSchur.reserve(nnzVec);


    for(ee=0; ee<nElem_global; ee++)
    {
        size1 = elems[ee]->forAssyVecVelo.size();
        size2 = elems[ee]->forAssyVecPres.size();

        for(ii=0; ii<size2; ii++)
        {
          row = elems[ee]->forAssyVecPres[ii];

          if(row != -1)
          {
            for(jj=0; jj<size1; jj++)
            {
              col = elems[ee]->forAssyVecVelo[jj];

              if(col != -1)
              {
                matKup.coeffRef(col, row) = 0.0;
                matKpu.coeffRef(row, col) = 0.0;
              }
            }
            for(jj=0; jj<size2; jj++)
            {
              col = elems[ee]->forAssyVecPres[jj];

              if(col != -1)
              {
                matSchur.coeffRef(row, col) = 0.0;
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
    } //for(ee=0;)


    cout << " setSolverDataForSemiImplicit () ... adding matrix entries for Lambdas " << endl;

    ImmersedIntegrationElement *lme;

    for(int bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      for(int iie=0; iie<ImmersedBodyObjects[bb]->getNumberOfElements(); iie++)
      {
        lme = ImmersedBodyObjects[bb]->ImmIntgElems[iie];

        size1 = lme->forAssyVec.size();

        ee = ImmersedElements[iie]->elemNums[0];

        size2 = elems[ee]->forAssyVecVelo.size();

        cout << " Lambdas " << iie << '\t' << size1 << '\t' << size2 << endl;

        printVector(ImmersedElements[iie]->forAssyVec);
        printVector(elems[ee]->forAssyVecVelo);

        for(ii=0; ii<size1; ii++)
        {
          row = ImmersedElements[iie]->forAssyVec[ii];

          if(row != -1)
          {
            row += totalDOF_Pres;

            for(jj=0; jj<size2; jj++)
            {
              col = elems[ee]->forAssyVecVelo[jj];

              if(col != -1)
              {
                matKpu.coeffRef(row, col) = 0.0;
                matKup.coeffRef(col, row) = 0.0;
              }
            }

            for(jj=0; jj<size1; jj++)
            {
              col = ImmersedElements[iie]->forAssyVec[jj];
              if(col != -1)
              {
                col += totalDOF_Pres;

                matSchur.coeffRef(row, col) = 0.0;
              }
            }

          }//if(row != -1)
        } //for(ii=0;)
      }
    } //for(ee=0;)


    matMuuInv.makeCompressed();
    matKup.makeCompressed();
    matKpu.makeCompressed();
    matSchur.makeCompressed();

    cout << " setSolverDataForSemiImplicit() ... ENDED " << endl;

    return 0;
}




int femINSmixed::updateMatricesSemiImplicit()
{
    int ii, jj, row, col, size1, size2;

    MatrixXd  Kup(npElemVelo*ndof, 4);

    matKup *= 0.0;
    matKpu *= 0.0;

    for(int ee=0; ee<nElem_global; ee++)
    {
        elems[ee]->StiffnessForSemiImpl(elemData, timeData, Kup);

        //printMatrix(Kup);

        //Assemble the matrix

        size1 = elems[ee]->forAssyVecVelo.size();
        size2 = elems[ee]->forAssyVecPres.size();

        for(ii=0; ii<size1; ii++)
        {
            row = elems[ee]->forAssyVecVelo[ii];
            if( row != -1 )
            {
                for(jj=0; jj<size2; jj++)
                {
                    col = elems[ee]->forAssyVecPres[jj];
                    if( col != -1 )
                    {
                        matKup.coeffRef(row, col) -= Kup(ii,jj);
                        matKpu.coeffRef(col, row) += Kup(ii,jj);
                    }
                }
            }
        }
    } //LoopElem

    cout << " adding matrix entries for Lambdas " << endl;

    myPoint coords_global, coords_local;
    vector<double>  Nf(9), dNf_du1(9), dNf_du2(9);
    MatrixXd  Khorz(2, 18);
    Khorz.setZero();

    for(int ime=0; ime<nImmersedElems; ime++)
    {
        size1 = ImmersedElements[ime]->forAssyVec.size();

        int elnum = ImmersedElements[ime]->elemNums[0];

        size2 = elems[elnum]->forAssyVecVelo.size();

        //cout << " Lambdas " << ime << '\t' << size1 << '\t' << size2 << endl;

        //elems[elnum]->findLocalCoordinates(nodeCoords, ImmersedElements[ime]->coords, coords_local);

        BernsteinBasisFunsQuad(2, coords_local(0), coords_local(1), &Nf[0], &dNf_du1[0], &dNf_du2[0]);

        for(ii=0; ii<9; ii++)
        {
            Khorz(0, ii*2)   = Nf[ii];
            Khorz(1, ii*2+1) = Nf[ii];
        }

        //printMatrix(Khorz);

        for(ii=0; ii<size1; ii++)
        {
          row = ImmersedElements[ime]->forAssyVec[ii];

          if(row != -1)
          {
            row += totalDOF_Pres;

            //cout << "row = " << ii << '\t' << row << endl;

            for(jj=0; jj<size2; jj++)
            {
              col = elems[elnum]->forAssyVecVelo[jj];

              if(col != -1)
              {
                matKpu.coeffRef(row, col) += Khorz(ii, jj);
                matKup.coeffRef(col, row) += Khorz(ii, jj);
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
    } //for(ee=0;)



    return 0;
}


//
int  femINSmixed::solveSchurComplementProblem()
{
    //cout << " SolverEigen::solverSchurExplicitDirect() ... STARTED " << endl;
    //  Kuu*x1 + Kup*x2 = f1
    //  Kpu*x1 + Kpp*x2 = f2
    //
    //  Kuu*x1 = f1 - Kup*x2
    //  S*x2 = f2 - Kpu*KuuInv*f1
    //
    //  S = Kpp - Kpu*KuuInv*Kup,  Schur complement

    // compute the Schur complement
    //if( (count % frequency) == 0 )
    //{
      //cout <<  "Updating the Schur complement solver " <<  endl;
//       matSchur = matKpp - matKpu*(matMuuInv*matKup);
//       matSchur.makeCompressed();
//
//       solverSchur.compute(matSchur);
    //}

    //cout << matS << endl;

    //auto  time1 = chrono::steady_clock::now();

    //printVector(rhsVecVelo);
    //printVector(rhsVecPres);

    VectorXd  r2 = matKpu*(matMuuInv*rhsVecVelo) - amDgammaDt*rhsVecPres;
    //printVector(r2);
    //solverSchur.compute(matSchur);
    presIncr = solverSchur.solve(r2);
    //presIncr = solverSchur.solveWithGuess(r2, presIncr);
    //printVector(presIncr);

    for(int ii=0; ii<totalDOF_Lambda; ii++)
    {
        lambdasIncr[ii] = presIncr[totalDOF_Pres+ii];
    }

    veloIncr = (matMuuInv/amDgammaDt)*(rhsVecVelo - matKup*presIncr);
    //printVector(veloIncr);

    //cout << solverSchur.info() << '\t' << solverSchur.error() << '\t' << solverSchur.iterations() << endl;

    //cout << " SolverEigen::solverSchurExplicitDirect() ... FINISHED " << endl;

    return 0;
}
//



/*
int  femINSmixed::solveSchurComplementProblem()
{
    cout << " SolverEigen::solverSchurExplicitDirect() ... STARTED " << endl;

    //printf("\n \t   Number of Velocity DOF    =  %5d\n\n", totalDOF_Velo);
    //printf("\n \t   Number of Pressure DOF    =  %5d\n\n", totalDOF_Pres);
    //printf("\n \t   Number of Immersed DOF    =  %5d\n\n", totalDOF_Lambda);
    //printf("\n \t   Number of Solid DOF       =  %5d\n\n", totalDOF_Solid);
    //printf("\n \t   Total number of DOF       =  %5d\n\n", totalDOF);

    int  ii, jj, k;

    //SparseMatrixXd  matK;

    matK.resize(totalDOF, totalDOF);                        //    matK.reserve(totalDOF*totalDOF*0.2);

    VectorXi  nnzVec(totalDOF);

    jj = ceil(totalDOF_Velo*0.1);
    //cout << " jj = " << jj << endl;

    for(ii=0; ii<totalDOF; ii++)
      nnzVec(ii) = jj;

    matK.reserve(nnzVec);

    for(k=0; k<matMuuInv.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(matMuuInv,k); it; ++it)
      {
        ii = it.row();
        jj = it.col();

        matK.coeffRef(ii, jj) = amDgammaDt/it.value();
      }
    }


    for(k=0; k<matKpu.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(matKpu,k); it; ++it)
      {
        ii = it.row()+totalDOF_Velo;
        jj = it.col();

        matK.coeffRef(ii, jj) = it.value();
        matK.coeffRef(jj, ii) = it.value();
      }
    }
    matK.makeCompressed();

    //cout << "bbbbbbbbbbbbbbbb" << endl;

    //VectorXd  rhsTemp(totalDOF),  solnTemp(totalDOF);
    rhsVec.resize(totalDOF),  soln.resize(totalDOF);

    for(ii=0; ii<totalDOF_Velo; ii++)
      rhsVec(ii) = rhsVecVelo(ii);

    for(ii=0; ii<totalDOF_Pres+totalDOF_Lambda; ii++)
      rhsVec(totalDOF_Velo+ii) = rhsVecPres(ii);

    //cout << "ccccccccccccccccc" << endl;

    //SuperLU<SparseMatrixXd > solverSchur;
    //SparseLU<SparseMatrix<double> > solverSchur;
    //solverSchur.compute(matK);

    //solnTemp = solverSchur.solve(rhsTemp);

    initialise_pardiso();
    factoriseAndSolve_pardiso();

    //cout << "ccccccccccccccccc" << endl;

    if (veloIncr.rows() != totalDOF_Velo)
        veloIncr.resize(totalDOF_Velo);

    if (presIncr.rows() != totalDOF_Pres)
        presIncr.resize(totalDOF_Pres);

    if (lambdasIncr.rows() != totalDOF_Lambda)
        lambdasIncr.resize(totalDOF_Lambda);

    for(ii=0; ii<totalDOF_Velo; ii++)
    {
      veloIncr[ii] = soln[ii];
    }

    //cout << "bbbbbbbbbbbbbbbb" << endl;

    for(ii=0; ii<totalDOF_Pres; ii++)
    {
      presIncr[ii] = soln[totalDOF_Velo+ii];
    }

    for(ii=0; ii<totalDOF_Lambda; ii++)
    {
        lambdasIncr[ii] = soln[totalDOF_Velo+totalDOF_Pres+ii];
    }

    //printVector(veloIncr);
    //printVector(presIncr);
    //printVector(lambdasIncr);


    //cout << " SolverEigen::solverSchurExplicitDirect() ... FINISHED " << endl;

    return 0;
}
*/






int  femINSmixed::solveSemiImplicit()
{
    cout << " Solving with the Semi-Implicit Scheme " << endl;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    velo.setZero();      veloPrev.setZero();      veloCur.setZero();
    veloDot.setZero();   veloDotPrev.setZero();   veloDotCur.setZero();
    pres.setZero();      presPrev.setZero();      presCur.setZero();
    presDot.setZero();   presDotPrev.setZero();   presDotCur.setZero();

    int  stepsCompleted=1;
    int  aa, bb, dd, ee, ii, jj, kk, count, row, col, ind, n1, n2;

    double  norm_velo=1.0, norm_pres=1.0, norm_velo_rel=1.0, norm_pres_rel=1.0;
    double  dtCrit=1.0e10, fact, fact1, fact2;
    double  dtgamma11, dtgamma12, norm_velo_diff, norm_pres_diff;
    double  am = 1.0;
    //double  gamma = 0.5+am;
    double  gamma = 1.0;

    vector<int>  nodeNums, nodeNumsPres, globalDOFnums;

    VectorXd  FlocalVelo(npElemVelo*ndof), FlocalPres(npElemPres);
    VectorXd  TotalForce(3);


    double tstart = 0.0;//omp_get_wtime();

    matMuuInv *= 0.0;
    for(ii=0; ii<totalDOF_Velo; ii++)
    {
        matMuuInv.coeffRef(ii,ii) = 1.0/globalMassVelo[assyForSolnVelo[ii]];
        //cout << ii << '\t' << matMuuInv.coeffRef(ii,ii) << endl;
    }
    matMuuInv.makeCompressed();

    updateMatricesSemiImplicit();

    matSchur = matKpu*(matMuuInv*matKup);
    matSchur.makeCompressed();

    /*
    cout << " matKpu " << endl;
    cout << matKpu << endl;
    cout << endl;    cout << endl;    cout << endl;
    cout << "matKup " << endl;
    cout << matKup << endl;
    cout << endl;    cout << endl;    cout << endl;
    cout << "matSchur "<< endl;
    cout << matSchur << endl;
    cout << endl;    cout << endl;    cout << endl;

    cout << "matSchur.rows()" << '\t' << "matSchur.cols()" << endl;
    cout << matSchur.rows() << '\t' << matSchur.cols() << endl;
    */

    //solverSchur.preconditioner().setDroptol(1.0e-2);
    //solverSchur.preconditioner().setFillfactor(2);
    //solverSchur.setTolerance(1.0e-8);

    solverSchur.compute(matSchur);

    //velo=veloApplied;
    //setInitialConditions();
    postProcess();

    timeNow=0.0, loadFactor=0.0;


    // resize and initialise rhs vectors
    // global full size vector
    rhsVecVeloTemp.resize(nNode_Velo*ndim);
    rhsVecVeloTemp.setZero();

    if( (npElemVelo == 7) || (npElemVelo == 11) || (npElemVelo == 20) )
      rhsVecPresTemp.resize(npElemPres*nElem_global);
    else
      rhsVecPresTemp.resize(nNode_Velo);
    rhsVecPresTemp.setZero();

    // used for the solution
    rhsVecVelo.resize(totalDOF_Velo);
    rhsVecVelo.setZero();

    rhsVecPres.resize(totalDOF_Pres+totalDOF_Lambda);
    rhsVecPres.setZero();


    //Time loop
    while( (myTime.cur <= (timeFinal-EPSILON)) && (stepsCompleted <= stepsMax) )
    {
        rhsVecVeloTemp.setZero();
        rhsVecPresTemp.setZero();

        //#pragma omp parallel private(ee, ii, jj, kk, dd) default(shared)
        //{
          //Loop over elements and compute the RHS and time step
          dtCrit=1.0e10;
          //#pragma omp for reduction(min : dtCrit) schedule(dynamic,100) //shared(ndim, nElem, rhsVecVelo, rhsVecPres, elemConn)
          for(ee=0; ee<nElem_global; ee++)
          {
            fact = timeNow - dt;

            dtCrit = min(dtCrit, elems[ee]->ResidualIncNavStokesSemiImpl(nodeCoords, fluidProperties, timeData, velo, veloPrev, veloDot, veloDotPrev, pres, presPrev, FlocalVelo, FlocalPres) );

            //printVector(FlocalVelo);
            //printVector(FlocalPres);

            nodeNums      = elems[ee]->nodeNums;
            globalDOFnums = elems[ee]->globalDOFnums;
            nodeNumsPres  = elems[ee]->nodeNumsPres;

            //int ii, jj, kk, dd;
            //Assemble the element vector
            //#pragma omp critical
            //{
              for(ii=0; ii<globalDOFnums.size(); ii++)
              {
                //cout << ii << '\t' << globalDOFnums[ii] << '\t' << FlocalVelo(ii) << endl;
                rhsVecVeloTemp(globalDOFnums[ii]) += FlocalVelo(ii);
              }

              for(ii=0; ii<nodeNumsPres.size(); ii++)
              {
                //cout << ii << '\t' << nodeNumsPres[ii] << '\t' << FlocalPres(ii) << endl;
                rhsVecPresTemp(nodeNumsPres[ii])  += FlocalPres(ii);
              }
            //}
          } //LoopElem

          // Actual time, critical time step scaled by used-specified CFL number
          myTime.dt = dtCrit*CFL;

          myTime.update();

          // update time functions
          for(auto& tmf : timeFunctions)
            tmf->update();

          //double  time1 = omp_get_wtime();
          //printf(" ==================================================================== \n");
          //printf(" Time step number     =  %d  \n", stepsCompleted);
          //printf(" Time step size       =  %f  \n", myTime.dt);
          //printf(" Current time         =  %f  \n", myTime.cur);
          //printf(" ==================================================================== \n");

          // set Dirichlet boundary conditions
          setBoundaryConditions();

          amDgammaDt = am/(gamma*myTime.dt);

          if(totalDOF_Lambda > 0)
            applyInterfaceTerms2D();

          // Add specified nodal force
          //addExternalForces(loadFactor);

          //rhsVecVelo.setZero();
          //rhsVecPres.setZero();

          fact = am/gamma - 1.0;
          for(ii=0; ii<totalDOF_Velo; ii++)
          {
            jj = assyForSolnVelo[ii];

            rhsVecVelo[ii]   =  rhsVecVeloTemp[jj];
            rhsVecVelo[ii]  +=  (fact*globalMassVelo[jj])*veloDotPrev[jj];
          }

          for(ii=0; ii<totalDOF_Pres; ii++)
          {
            rhsVecPres[ii]  =  rhsVecPresTemp[assyForSolnPres[ii]];
          }

          for(ii=0; ii<totalDOF_Lambda; ii++)
          {
            rhsVecPres[totalDOF_Pres+ii]  =  rhsVecPresTemp[totalDOF_Pres+ii];
          }


          solveSchurComplementProblem();

          // compute the solution at t_{n+1}

          for(ii=0; ii<totalDOF_Velo; ii++)
          {
            velo[assyForSolnVelo[ii]]   +=  veloIncr[ii];
          }

          for(ii=0; ii<totalDOF_Pres; ii++)
          {
            pres[assyForSolnPres[ii]]   +=  presIncr[ii];
          }

          for(ii=0; ii<totalDOF_Lambda; ii++)
          {
            lambdas[ii]   +=  lambdasIncr[ii];
          }

          // apply boundary conditions
          int ii, dof;
          for(ii=0; ii<dofs_specified_velo.size(); ii++)
          {
            dof = dofs_specified_velo[ii];
            velo[dof] += veloApplied[dof] ;
          }

          for(ii=0; ii<dofs_specified_pres.size(); ii++)
          {
            dof = dofs_specified_pres[ii];
            pres[dof] += presApplied[dof] ;
          }

          // calculate acceleration
          veloDot = (1.0/gamma/myTime.dt)*(velo-veloPrev) + ((gamma-1.0)/gamma)*veloDotPrev;

          // compute the norms and store the variables
          norm_velo = 0.0;
          norm_velo_diff = 0.0;
          for(ii=0; ii<velo.size(); ii++)
          {
            norm_velo += velo[ii]*velo[ii];

            fact1 = velo[ii]-veloPrev[ii];
            norm_velo_diff += fact1*fact1;
          }

          norm_pres = 0.0;
          norm_pres_diff = 0.0;
          for(ii=0; ii<pres.size(); ii++)
          {
            norm_pres += pres[ii]*pres[ii];

            fact1 = pres[ii]-presPrev[ii];
            norm_pres_diff += fact1*fact1;
          }
        //}

        if( std::isnan(norm_velo) || std::isnan(norm_pres) )
        {
          cerr << " NAN encountered in the solution ... " << endl;
          cerr << " Program has been terminated ... " << endl;
          exit(-1);
        }

        if(stepsCompleted > 1)
        {
          norm_velo_rel = sqrt(norm_velo_diff/norm_velo);
          norm_pres_rel = sqrt(norm_pres_diff/norm_pres);
        }

        // store the variables
        //veloPrev3  = veloPrev2;
        veloPrev2   = veloPrev;
        veloPrev    = velo;
        veloDotPrev = veloDot;

        //presPrev3  = presPrev2;
        presPrev2  = presPrev;
        presPrev   = pres;

        if( (stepsCompleted%outputFreq == 0) || (norm_velo_rel < conv_tol) )
        {
          cout << " Steps completed = " << stepsCompleted << '\t'
             << "Time step = " << myTime.dt << '\t'
             << "Current time = " << myTime.cur << '\t'
             << " velocity norm = " << norm_velo_rel << '\t'
             << " pressure norm = " << '\t' << norm_pres_rel << endl;


            postProcess();
        }

        reacVec = rhsVecVeloTemp;
        writeOutputDataPatches();

        if(norm_velo_rel < conv_tol)
        {
          cout << " Solution convdataged below the specified tolerance " << endl;
          break;
        }

        stepsCompleted++;
    } //Time loop


    //double tend = omp_get_wtime(), duration=tend-tstart;
    double tend = 0.0, duration=tend-tstart;

    cout << " \n \n Total time taken = " << duration << " seconds \n\n" << endl;

    return 0;
}





#include "femINSmixed.h"
#include "util.h"
#include "ElementBase.h"
#include <chrono>
#include "KimMoinFlow.h"
#include "BasisFunctionsBernstein.h"
#include "MyTime.h"
#include "TimeFunction.h"

extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime                 myTime;


int femINSmixed::setSolver(int slv)
{
    //Eigen::initParallel();
    Eigen::setNbThreads(0);

    prepareMatrixPattern();

    cout << " slv = " << slv << endl;

    if(slv == IMPLICIT)
    {
      setSolverDataForFullyImplicit();
    }
    else if(slv ==  SEMI_IMPLICIT)
    {
      calcMassMatrixForExplicitDynamics();
      setSolverDataForSemiImplicit();
    }
    else if(slv ==  EXPLICIT)
    {
      calcMassMatrixForExplicitDynamics();
    }
    else
    {
      cerr << "femINSmixed::setSolver ... slv not supported " <<  endl;
      exit(-1);
    }

    return 0;
}




int femINSmixed::prepareMatrixPattern()
{
    cout <<  "\n     femINSmixed::prepareMatrixPattern()  .... STARTED ...\n" <<  endl;

    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    cout << " nNode_Velo     = " << '\t' << nNode_Velo  << endl;
    cout << " nNode_Pres     = " << '\t' << nNode_Pres  << endl;

    cout << " nDBC_Velo      = " << '\t' << nDBC_Velo  << endl;
    cout << " nDBC_Pres      = " << '\t' << nDBC_Pres  << endl;

    vector<vector<int> >  IDvelo;
    vector<vector<bool> >  NodeTypeVelo;


    NodeTypeVelo.resize(nNode_Velo);
    IDvelo.resize(nNode_Velo);

    int  ee, ii, jj, kk, nn, dof, ind;

    for(ii=0;ii<nNode_Velo;++ii)
    {
      NodeTypeVelo[ii].resize(ndim);
      IDvelo[ii].resize(ndim);

      for(jj=0;jj<ndim;++jj)
      {
        NodeTypeVelo[ii][jj] = false;
        IDvelo[ii][jj] = -1;
      }
    }

    // fix the specified Dirichlet BCs
    for(ii=0; ii<nDBC_Velo; ++ii)
    {
      NodeTypeVelo[DirichletBCsVelo[ii][0]][DirichletBCsVelo[ii][1]] = true;
    }

    totalDOF_Velo = 0;
    for(ii=0;ii<nNode_global;++ii)
    {
      for(jj=0;jj<ndim;++jj)
      {
        if(!NodeTypeVelo[ii][jj])
        {
          IDvelo[ii][jj] = totalDOF_Velo++;
          assyForSolnVelo.push_back(ii*ndim+jj);
        }
      }
    }


    if(npElemVelo == 7)
      ind = nNode_Pres;
    else
      ind = nNode_Velo;

    vector<bool>  NodeTypePres(ind,  false);
    vector<int>   IDpres(ind,  -1);

    // fix the pressure at all the mid nodes
    if(npElemVelo != 7)
    {
      for(ii=0; ii<nNode_Velo; ++ii)
      {
        if(midNodeData[ii][0])
          NodeTypePres[ii] = true;
      }

      if( (ndim == 2) && (npElemVelo == 9) )
      {
        for(ee=0; ee<nElem_global; ++ee)
        {
          NodeTypePres[elemConn[ee][8]] = true;
        }
      }
    }

    // fix the specified pressure DOFs
    for(ii=0; ii<nDBC_Pres; ++ii)
    {
      NodeTypePres[DirichletBCsPres[ii][0]] = true;
    }

    totalDOF_Pres = 0;
    if(npElemVelo == 7)
    {
      for(ii=0;ii<nNode_Pres;++ii)
      {
        if(!NodeTypePres[ii])
        {
          IDpres[ii] = totalDOF_Pres++;
          assyForSolnPres.push_back(ii);
        }
      }
    }
    else
    {
      for(ii=0;ii<nNode_Velo;++ii)
      {
        if(!NodeTypePres[ii])
        {
          IDpres[ii] = totalDOF_Pres++;
          assyForSolnPres.push_back(ii);
        }
      }
    }

    cout << " Mesh statimeIntegrationSchemetics .....\n" << endl;
    cout << " nElem          = " << '\t' << nElem_global << endl;
    cout << " nNode_Velo     = " << '\t' << nNode_Velo  << endl;
    cout << " nNode_Pres     = " << '\t' << nNode_Pres  << endl;
    cout << " npElemVelo     = " << '\t' << npElemVelo << endl;
    cout << " npElemPres     = " << '\t' << npElemPres << endl;
    cout << " ndof           = " << '\t' << ndof << endl;
    cout << " Velocity DOF   = " << '\t' << totalDOF_Velo << endl;
    cout << " Pressure DOF   = " << '\t' << totalDOF_Pres << endl;


    for(ee=0;ee<nElem_global;++ee)
    {
      npElemVelo = elems[ee]->nodeNums.size();

      ind = ndim*npElemVelo;

      elems[ee]->forAssyVecVelo.resize(ind);

      for(ii=0; ii<npElemVelo; ++ii)
      {
        ind = ndim*ii;

        kk = elems[ee]->nodeNums[ii];

        for(jj=0;jj<ndim;++jj)
        {
          elems[ee]->forAssyVecVelo[ind+jj] = IDvelo[kk][jj];
        }
      }

      elems[ee]->forAssyVecPres.resize(npElemPres);

      if(npElemVelo == 7)
      {
        for(ii=0; ii<npElemPres; ++ii)
        {
          kk = elems[ee]->nodeNumsPres[ii];

          elems[ee]->forAssyVecPres[ii] = IDpres[kk];
        }
      }
      else
      {
        for(ii=0; ii<npElemPres; ++ii)
        {
          kk = elems[ee]->nodeNums[ii];

          elems[ee]->forAssyVecPres[ii] = IDpres[kk];
        }
      }

      //printVector(elems[ee]->forAssyVecVelo);
      //printVector(elems[ee]->forAssyVecPres);
    }


    bool pp=false;
    //pp=true;
    if(pp)
    {
       printf("   ID array \n\n");
       for(ii=0;ii<nNode_global;++ii)
       {
          for(jj=0;jj<ndim;++jj)
            cout << '\t' << IDvelo[ii][jj];
          cout << endl;
       }
       printf("\n\n\n");

       printf("  assyForSolnVelo array \n\n");
       for(ii=0;ii<totalDOF_Velo;++ii)
       {
          cout << assyForSolnVelo[ii] << endl;
       }
       printf("\n\n\n");
    }

    printf("\n element DOF values initialised \n\n");

    // matrix pattern needs to be prepared only for the implicit solver

    // remove data objects

    //ID.clear();
    //forAssyMat.clear();

    prepareImmersedSolids();

    totalDOF = totalDOF_Velo + totalDOF_Pres + totalDOF_Lambda + totalDOF_Solid;

    printf("\n \t   Number of Velocity DOF    =  %5d\n\n", totalDOF_Velo);
    printf("\n \t   Number of Pressure DOF    =  %5d\n\n", totalDOF_Pres);
    printf("\n \t   Number of Immersed DOF    =  %5d\n\n", totalDOF_Lambda);
    printf("\n \t   Number of Solid DOF       =  %5d\n\n", totalDOF_Solid);
    printf("\n \t   Total number of DOF       =  %5d\n\n", totalDOF);

    printf("\n     femINSmixed::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return 1;
}







int femINSmixed::addExternalForces()
{
    int  nn, dof, ii, ind;
    double specVal=0.0;

    VectorXd  vecTemp, Flocal;
    vecTemp.setZero();

    // specified nodal forces
    for(ii=0;ii<nodeForcesData.size();++ii)
    {
      nn  = (int) (nodeForcesData[ii][0] - 1);
      dof = (int) (nodeForcesData[ii][1] - 1);
      specVal = nodeForcesData[ii][2];

      ind = nn*ndof+dof;

      vecTemp[ind] += specVal*loadFactor;
    }
    //printVector(vecTemp);

    return 0;
}




// set the off-diagonal terms for the solver
int femINSmixed::setSolverDataForFullyImplicit()
{
    cout <<  " femINSmixed::setSolverDataForFullyImplicit() ... STARTED " << endl;

    int  ee, ii, jj, size1, size2, row, col;
    vector<int>  vecIntTemp(10);

    cout << " Total DOF   = " << '\t' << totalDOF << endl;

    rhsVec.resize(totalDOF);


    matK.setZero();
    matK.resize(totalDOF, totalDOF);

    VectorXi  nnzVec(totalDOF);

    jj = ceil(totalDOF*0.1);
    jj = 1000;

    for(ii=0; ii<totalDOF; ii++)
      nnzVec(ii) = jj;

    matK.reserve(nnzVec);

    for(ee=0; ee<nElem_global; ee++)
    {
        size1 = elems[ee]->forAssyVecVelo.size();
        size2 = elems[ee]->forAssyVecPres.size();

        //printVector(elems[ee]->forAssyVecVelo);
        //printVector(elems[ee]->forAssyVecPres);

        for(ii=0; ii<size1; ii++)
        {
          row = elems[ee]->forAssyVecVelo[ii];

          if(row != -1)
          {
            for(jj=0; jj<size1; jj++)
            {
              col = elems[ee]->forAssyVecVelo[jj];

              //cout << ii << '\t' << jj << '\t' << row << '\t' << col << endl;

              if(col != -1)
              {
                matK.coeffRef(row, col) = 0.0;
              }
            }

            for(jj=0; jj<size2; jj++)
            {
              col = elems[ee]->forAssyVecPres[jj];

              if(col != -1)
              {
                col += totalDOF_Velo;

                matK.coeffRef(row, col) = 0.0;
                matK.coeffRef(col, row) = 0.0;
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
    } //for(ee=0;)

    cout << " setSolverDataForFullyImplicit () ... adding matrix entries for Lambdas " << endl;
    
    int totalDOF_VeloPres = totalDOF_Velo+totalDOF_Pres;

    
    ImmersedIntegrationElement *lme;
    int elnum, bb, ime;

    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      for(ime=0; ime<ImmersedBodyObjects[bb]->getNumberOfElements(); ime++)
      {
        lme = ImmersedBodyObjects[bb]->ImmIntgElems[ime];

        //cout << bb << '\t' << aa << '\t' << lme->isActive() << '\t' << lme->gausspoints.size() << endl;

        if( lme->isActive() )
        {
          size1 = lme->forAssyVec.size();

          for(int gp=0; gp<lme->gausspoints.size(); gp++)
          {
            elnum = lme->elemNums[gp];
            size2 = elems[elnum]->forAssyVecVelo.size();

            //cout << " Lambdas " << ime << '\t' << ee << '\t' << size1 << '\t' << size2 << endl;

            //printVector(lme->forAssyVec);
            //printVector(elems[elnum]->forAssyVecVelo);

            for(ii=0; ii<size1; ii++)
            {
              row = lme->forAssyVec[ii];

              if(row != -1)
              {
                row += (totalDOF_VeloPres);

                for(jj=0; jj<size2; jj++)
                {
                  col = elems[elnum]->forAssyVecVelo[jj];

                  if(col != -1)
                  {
                    matK.coeffRef(row, col) = 0.0;
                    matK.coeffRef(col, row) = 0.0;
                  }
                }
              }//if(row != -1)
            } //for(ii=0;)
          }
        }
      }
    }

    matK.makeCompressed();

    //initialise_pardiso();
    solver.analyzePattern(matK);
    solver.factorize(matK);
    //solver.compute(matK);

    cout <<  " femINSmixed::setSolverDataForFullyImplicit() ... ENDED " << endl;

    return 0;
}


int  femINSmixed::solveFullyImplicit()
{
    cout << " Solving with the Fully-Implicit Scheme " << endl;

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////

    int  stepsCompleted=0;
    int  nsize_velo = rhsVecVelo.rows();
    int  nsize_pres = rhsVecPres.rows();
    int  aa, bb, ee, ii, jj, kk, count, row, col, ind, n1, n2, size1, size2;

    double  fact, fact1, fact2;

    VectorXd  td(100), reacVec(nNode_Velo*ndim+nNode_Pres);
    VectorXd  TotalForce(3);
    VectorXd  FlocalVelo(30), FlocalPres(10);
    ind = npElemVelo*ndim;
    MatrixXd  Kuu(ind,ind), Kup(ind,npElemPres);

    vector<int>  vecTempInt;

    cout << " timeIntegrationScheme = " << timeIntegrationScheme << '\t' << spectralRadius << endl;
 
    SetTimeParametersFluid(timeIntegrationScheme, spectralRadius, dt, td);

    KimMoinFlowUnsteadyNavierStokes  analy(elemData[0], elemData[1], 1.0);

    double tstart = 0.0;// omp_get_wtime();

    timeNow=0.0;

    //setInitialConditions();
    //Time loop
    for(int tstep=0; tstep<stepsMax; tstep++)
    {
        timeNow = timeNow + dt;
 
        cout << " Time = " << timeNow << endl;
 
        if(timeNow > timeFinal)
          break;
 
        //
        if(timeNow <= 1.0)
        {
          loadFactor = timeNow;
          //loadFactor = stepsCompleted/5000.0;
        }
        else
        {
          loadFactor = 1.0;
        }
        //
        //loadFactor = 1.0;

        // set boundary conditions
        assignBoundaryConditions();

        // add boundary conditions
        //addBoundaryConditions();

        // update Immersed Solids
        //updateImmersedSolid();

        // recompute matrix pattern and update the solver
        if(IB_MOVED)
        {
          prepareImmersedSolids();
          setSolverDataForFullyImplicit();
        }

        cout << " Iteration loop " << endl;

        for(int iter=0; iter<iterationsMax; iter++)
        {
            cout << " Iteration ... " << iter << endl;

            veloDot    = td[9]*velo + td[10]*veloPrev + td[15]*veloDotPrev ;

            veloCur    = td[2]*velo    + (1.0-td[2])*veloPrev; // velocity
            presCur    = td[2]*pres    + (1.0-td[2])*presPrev; // pressure
            lambdasCur = td[2]*lambdas + (1.0-td[2])*lambdasPrev; // Lagrange multipliers
            veloDotCur = td[1]*veloDot + (1.0-td[1])*veloDotPrev;


            matK *= 0.0;
            rhsVec.setZero();
            reacVec.setZero();

            cout << " Element loop " << endl;

            // loop over elements and compute matrix and residual
            for(ee=0; ee<nElem_global; ee++)
            {
                //elems[ee]->StiffnessAndResidualFullyImplicit(nodeCoords, elemData, &td[0], veloPrev, veloPrev2, veloCur, veloDotCur, presCur, Kuu,  Kup, FlocalVelo, FlocalPres, dt, timeNow);

                //printMatrix(Kuu);
                //printMatrix(Kup);
                //printVector(FlocalVelo);
                //printVector(FlocalPres);

                size1 = elems[ee]->forAssyVecVelo.size();
                size2 = elems[ee]->forAssyVecPres.size();

                //cout << "Applying boundary conditions " << ee << endl;
                if(iter == 0)
                {
                    // apply boundary conditions

                    // applied velocity dof
                    for(ii=0; ii<size1; ii++)
                    {
                        //printVector(elems[ee]->forAssyVecVelo);
                        //printVector(elems[ee]->forAssyVecPres);
                        //printVector(elems[ee]->globalDOFnums);

                        aa = elems[ee]->forAssyVecVelo[ii];

                        if(aa == -1) // this DOF has a prescibed value
                        {
                            fact = veloApplied[elems[ee]->globalDOFnums[ii]];

                            for(jj=0; jj<size1; jj++)
                            {
                                if( elems[ee]->forAssyVecVelo[jj] != -1 )
                                {
                                    FlocalVelo(jj) -= Kuu(jj, ii) * fact;
                                }
                            }

                            for(jj=0; jj<size2; jj++)
                            {
                                if( elems[ee]->forAssyVecPres[jj] != -1 )
                                {
                                    FlocalPres(jj) -= Kup(ii, jj) * fact;
                                }
                            }
                        }
                    }

                    // applied pressure
                    for(ii=0; ii<size2; ii++)
                    {
                        aa = elems[ee]->forAssyVecPres[ii];
                        if(aa == -1)                       // this DOF has a prescribed value
                        {
                            fact = presApplied[elems[ee]->nodeNums[ii]];

                            for(jj=0; jj<size1; jj++)
                            {
                                if( elems[ee]->forAssyVecVelo[jj] != -1 )
                                {
                                    FlocalVelo(jj) -= Kup(jj, ii) * fact;
                                }
                            }
                        }
                    }
                } // if(iter == 0)
                //cout << "Applying boundary conditions " << ee << endl;

                //cout << "Assembling matrices and vectors " << ee << endl;

                // assemble matrices and vectors

                for(ii=0; ii<size1; ii++)
                {
                  row = elems[ee]->forAssyVecVelo[ii];

                  reacVec[elems[ee]->globalDOFnums[ii]] += FlocalVelo[ii];


                  if(row != -1)
                  {
                    rhsVec(row)  += FlocalVelo(ii);

                    for(jj=0; jj<size1; jj++)
                    {
                      col = elems[ee]->forAssyVecVelo[jj];

                      if(col != -1)
                      {
                        matK.coeffRef(row, col) += Kuu(ii, jj);
                      }
                    }

                    for(jj=0; jj<size2; jj++)
                    {
                      col = elems[ee]->forAssyVecPres[jj];

                      if(col != -1)
                      {
                        col += totalDOF_Velo;

                        matK.coeffRef(row, col) += Kup(ii, jj);
                        matK.coeffRef(col, row) += Kup(ii, jj);
                      }
                    }
                  }//if(row != -1)
                } //for(ii=0;)

                //cout << "bbbbbbbbbbbbbbbb.... " << ee << endl;

                for(ii=0; ii<size2; ii++)
                {
                  row = elems[ee]->forAssyVecPres[ii];

                  if(row != -1)
                  {
                    row +=  totalDOF_Velo;

                    rhsVec(row)  += FlocalPres(ii);
                  }//if(row != -1)
                } //for(ii=0;)

                //cout << "ccccccccccccccc.... " << ee << endl;
            } //Element Loop

            
            cout << "Adding Lagrange multipliers " << endl;

            if(totalDOF_Lambda > 0)
              applyInterfaceTerms2D();
            
            cout << "Adding boundary conditions " << endl;

            // add boundary conditions
            if(iter == 0)
              addBoundaryConditions();

            rhsNorm = rhsVec.norm();

            printf(" Iteration : %2d ... RHS norm = %E \n", iter, rhsNorm);
            
            //printVector(rhsVec);

            if(rhsNorm < 1.0e-8)
            {
                //cout << " Solution converged below the specified tolerance " << endl;
                break;
            }
            else
            {
                //cout << "Solving the matrix system " << endl;

                solver.compute(matK);
                soln = solver.solve(rhsVec);

                //factoriseAndSolve_pardiso();

                //printVector(rhsVec);
                //printVector(soln);
                printf(" Solution completed \n");

                for(ii=0; ii<totalDOF_Velo; ii++)
                {
                  velo[assyForSolnVelo[ii]]   +=  soln[ii];
                }

                for(ii=0; ii<totalDOF_Pres; ii++)
                {
                  pres[assyForSolnPres[ii]]   +=  soln[totalDOF_Velo+ii];
                }
                
                int totalDOF_VeloPres = totalDOF_Velo + totalDOF_Pres;

                for(ii=0; ii<totalDOF_Lambda; ii++)
                {
                  lambdas[ii]   +=  soln[totalDOF_VeloPres+ii];
                }
                //printVector(velo);
                //printf("\n\n");
                //printVector(pres);
                //printf("\n\n");
                //printVector(lambdas);
                //for(ii=0; ii<totalDOF_Lambda/2; ii++)
                    //cout << ii << '\t' << lambdas(2*ii) << '\t' << lambdas(2*ii+1) << endl;
                //printf("\n\n");
            }
        } //Iteration Loop

        cout << "Postprocessing " << endl;
        postProcess();
        cout << "aaaaaaaaaaaaaaaaa " << endl;

        double TotalForce[2] = {0.0, 0.0};
        for(ii=0; ii<nOutputFaceLoads; ++ii)
        {
          TotalForce[0] +=  reacVec[outputEdges[ii][0]*2];
          TotalForce[1] +=  reacVec[outputEdges[ii][0]*2+1];
        }

        //fout_convdata <<  timeNow << '\t' << TotalForce[0] << '\t' << TotalForce[1] << endl;
        cout << endl; cout << endl;

        cout << "bbbbbbbbbbbbbbbb " << endl;

        veloPrev2    =  veloPrev;
        veloPrev     =  velo;
        veloDotPrev  =  veloDot;
        presPrev     =  pres;
        lambdasPrev  =  lambdas;

        cout << "aaaaaaaaaaaaaaaaa " << endl;
    } //Time loop

    //double tend = omp_get_wtime(), duration=tend-tstart;
    double tend = 0.0, duration=tend-tstart;

    cout << " \n \n Total time taken = " << duration << " seconds \n\n" << endl;


    /*
    cout << " Computing errors \n " << endl;
    double totalError = 0.0;
    //cout << " index = " << index << endl;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(ee=0; ee<nElem; ee++)
      {
        //Compute the element force vector, including residual force
        totalError += elems[ee]->CalculateError(nodeCoords, elemData, timeData, velo, veloDot, pres, 5.0, index);
      }

      totalError = sqrt(totalError);

      if(index == 0)
        printf(" \n\n \t L2 Error in X-velocity = %12.6E \n\n " , totalError);
      else if(index == 1)
        printf(" \n\n \t L2 Error in Y-velocity = %12.6E \n\n " , totalError);
      else if(index == 2)
        printf(" \n\n \t L2 Error in pressure   = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t H1 Error in velocity   = %12.6E \n\n " , totalError);
    }
    */


    return 0;
}




int femINSmixed::updateImmersedSolid()
{
    if( ImmersedBodyObjects.size() > 0)
    {
      for(int bb=0; bb<ImmersedBodyObjects.size(); bb++)
      {
        ImmersedBodyObjects[0]->updatePointPositions(timeNow);
      }
      IB_MOVED = true;
    }

    return 0;
}





int femINSmixed::updateImmersedSolid(int id, int ndof, VectorXd& disp, VectorXd& velo)
{
    ImmersedBodyObjects[id]->updatePointPositions(ndof, disp);

    IB_MOVED = true;
    
    return 0;
}



int femINSmixed::calcStiffnessAndResidual()
{
    int  stepsCompleted=0;
    int  nsize_velo = rhsVecVelo.rows();
    int  nsize_pres = rhsVecPres.rows();
    int  aa, bb, ee, ii, jj, kk, count, row, col, ind, n1, n2, size1, size2;

    double  fact, fact1, fact2;

    VectorXd  reacVec(nNode_Velo*ndim+nNode_Pres);
    VectorXd  FlocalVelo(30), FlocalPres(10);
    ind = npElemVelo*ndim;
    MatrixXd  Kuu(ind,ind), Kup(ind,npElemPres);

    vector<int>  vecTempInt;

    matK *= 0.0;
    rhsVec.setZero();
    reacVec.setZero();

    cout << " Element loop " << endl;

    // loop over elements and compute matrix and residual
    for(ee=0; ee<nElem_global; ee++)
    {
        //elems[ee]->StiffnessAndResidualFullyImplicit(nodeCoords, elemData, &td[0], veloPrev, veloPrev2, veloCur, veloDotCur, presCur, Kuu,  Kup, FlocalVelo, FlocalPres, dt, loadFactor);

        //printMatrix(Kuu);
        //printMatrix(Kup);
        //printVector(FlocalVelo);
        //printVector(FlocalPres);

        size1 = elems[ee]->forAssyVecVelo.size();
        size2 = elems[ee]->forAssyVecPres.size();

        //cout << "Applying boundary conditions " << ee << endl;
        if(firstIteration)
        {
            // apply boundary conditions

            // applied velocity dof
            for(ii=0; ii<size1; ii++)
            {
                //printVector(elems[ee]->forAssyVecVelo);
                //printVector(elems[ee]->forAssyVecPres);
                //printVector(elems[ee]->globalDOFnums);

                aa = elems[ee]->forAssyVecVelo[ii];

                if(aa == -1) // this DOF has a prescibed value
                {
                    fact = veloApplied[elems[ee]->globalDOFnums[ii]];

                    for(jj=0; jj<size1; jj++)
                    {
                        if( elems[ee]->forAssyVecVelo[jj] != -1 )
                        {
                            FlocalVelo(jj) -= Kuu(jj, ii) * fact;
                        }
                    }

                    for(jj=0; jj<size2; jj++)
                    {
                        if( elems[ee]->forAssyVecPres[jj] != -1 )
                        {
                            FlocalPres(jj) -= Kup(ii, jj) * fact;
                        }
                    }
                }
            }

            // applied pressure
            for(ii=0; ii<size2; ii++)
            {
                aa = elems[ee]->forAssyVecPres[ii];
                if(aa == -1)                       // this DOF has a prescribed value
                {
                    fact = presApplied[elems[ee]->nodeNums[ii]];

                    for(jj=0; jj<size1; jj++)
                    {
                        if( elems[ee]->forAssyVecVelo[jj] != -1 )
                        {
                            FlocalVelo(jj) -= Kup(jj, ii) * fact;
                        }
                    }
                }
            }
        } // if(iter == 0)
        //cout << "Applying boundary conditions " << ee << endl;

        //cout << "Assembling matrices and vectors " << ee << endl;

        // assemble matrices and vectors

        for(ii=0; ii<size1; ii++)
        {
            row = elems[ee]->forAssyVecVelo[ii];

            reacVec[elems[ee]->globalDOFnums[ii]] += FlocalVelo[ii];


            if(row != -1)
            {
                rhsVec(row)  += FlocalVelo(ii);

                for(jj=0; jj<size1; jj++)
                {
                    col = elems[ee]->forAssyVecVelo[jj];

                    if(col != -1)
                    {
                        matK.coeffRef(row, col) += Kuu(ii, jj);
                    }
                }

                for(jj=0; jj<size2; jj++)
                {
                    col = elems[ee]->forAssyVecPres[jj];

                    if(col != -1)
                    {
                        col += totalDOF_Velo;

                        matK.coeffRef(row, col) += Kup(ii, jj);
                        matK.coeffRef(col, row) += Kup(ii, jj);
                    }
                }
            }//if(row != -1)
        } //for(ii=0;)

        //cout << "bbbbbbbbbbbbbbbb.... " << ee << endl;

        for(ii=0; ii<size2; ii++)
        {
            row = elems[ee]->forAssyVecPres[ii];

            if(row != -1)
            {
                row +=  totalDOF_Velo;

                rhsVec(row)  += FlocalPres(ii);
            }//if(row != -1)
        } //for(ii=0;)

        //cout << "ccccccccccccccc.... " << ee << endl;
    } //Element Loop

            
    cout << "Adding Lagrange multipliers " << endl;

    if(totalDOF_Lambda > 0)
        applyInterfaceTerms2D();
            
    cout << "Adding boundary conditions " << endl;

    // add boundary conditions
    if(firstIteration)
        addBoundaryConditions();

    rhsNorm = rhsVec.norm();

    printf(" RHS norm = %E \n", rhsNorm);


    return 0;
}



int femINSmixed::factoriseSolveAndUpdate()
{
    solver.compute(matK);
    soln = solver.solve(rhsVec);

    //factoriseAndSolve_pardiso();

    //printVector(rhsVec);
    //printVector(soln);
    printf(" Solution completed \n");

    int ii;
    for(ii=0; ii<totalDOF_Velo; ii++)
    {
        velo[assyForSolnVelo[ii]]   +=  soln[ii];
    }

    for(ii=0; ii<totalDOF_Pres; ii++)
    {
        pres[assyForSolnPres[ii]]   +=  soln[totalDOF_Velo+ii];
    }
                
    int totalDOF_VeloPres = totalDOF_Velo + totalDOF_Pres;

    for(ii=0; ii<totalDOF_Lambda; ii++)
    {
        lambdas[ii]   +=  soln[totalDOF_VeloPres+ii];
    }
    //printVector(velo);
    //printf("\n\n");
    //printVector(pres);
    //printf("\n\n");
    //printVector(lambdas);
    //for(ii=0; ii<totalDOF_Lambda/2; ii++)
        //cout << ii << '\t' << lambdas(2*ii) << '\t' << lambdas(2*ii+1) << endl;

    return 0;
}



int femINSmixed::solveStep(int max_iters)
{
    for(int iter=0; iter<max_iters; iter++)
    {
      firstIteration = (iter == 0);

      cout << " aaaaaaaaaaa " << endl;
      calcStiffnessAndResidual();

      cout << " bbbbbbbbbbb " << endl;
      if( converged() )
        break;

      cout << " ccccccccccc " << endl;
      factoriseSolveAndUpdate();

      cout << " ddddddddddd " << endl;
      updateIterStep();
    }
    
    if( !converged() )
      return 0;
  
  return 1;
}




int femINSmixed::calcMassMatrixForExplicitDynamics()
{
    cout << " femINSmixed::calcMassMatrixForExplicitDynamics() STARTED " << endl;

    int  dd, ee, ii, jj;

    VectorXd  Flocal1(npElemVelo*ndim), Flocal2(npElemPres);

    // Compute global mass matrices
    //The Mass is assumed to be lumped so that the mass matrix is diagonal

    globalMassVelo.setZero();
    globalMassPres.setZero();

    for(ee=0; ee<nElem_global; ++ee)
    {
        // compute mass matrix and assemble it
        //elems[ee]->MassMatrices(nodeCoords, elemData, Flocal1, Flocal2);

        //Assemble the element mass to the global mass
        for(ii=0; ii<npElemVelo; ++ii)
        {
            jj = elemConn[ee][ii]*ndim ;

            for(dd=0; dd<ndim; dd++)
              globalMassVelo(jj+dd)   +=  Flocal1(ii);
        }

        for(ii=0; ii<npElemPres; ++ii)
        {
            globalMassPres(elemConn[ee][ii]) += Flocal2(ii);
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
        //elems[ee]->StiffnessForSemiImpl(elemData, timeData, Kup);

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
                        matKup.coeffRef(row, col) += Kup(ii,jj);
                        matKpu.coeffRef(col, row) -= Kup(ii,jj);
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

//     printVector(rhsVecVelo);
//     printVector(rhsVecPres);

    VectorXd  r2 = -amDgammaDt*rhsVecPres + matKpu*(matMuuInv*rhsVecVelo);
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
//     printVector(veloIncr);

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
    int  nsize_velo = rhsVecVelo.rows();
    int  nsize_pres = rhsVecPres.rows();
    int  aa, bb, dd, ee, ii, jj, kk, count, row, col, ind, n1, n2;

    double  norm_velo=1.0, norm_pres=1.0;
    double  dtCrit=1.0e10, fact, fact1, fact2;
    double  dtgamma11, dtgamma12, norm_velo_diff, norm_pres_diff;
    double  am = 1.0;
    double  gamma = 0.5+am;


    VectorXd  FlocalVelo(npElemVelo*ndof), FlocalPres(npElemPres);
    VectorXd  TotalForce(3);

    rhsVecVeloTemp.resize(nsize_velo), rhsVecPresTemp.resize(nsize_pres);

    rhsVecVelo.resize(totalDOF_Velo);
    rhsVecPres.resize(totalDOF_Pres+totalDOF_Lambda);


    double tstart = 0.0;//omp_get_wtime();

    matMuuInv *= 0.0;
    for(ii=0; ii<totalDOF_Velo; ii++)
    {
        matMuuInv.coeffRef(ii,ii) = 1.0/globalMassVelo[assyForSolnVelo[ii]];
    }
    matMuuInv.makeCompressed();

    updateMatricesSemiImplicit();

    matSchur = matKpu*(matMuuInv*matKup);
    matSchur.makeCompressed();

    //cout << matKpu << endl;
    //cout << endl;    cout << endl;    cout << endl;
    //cout << matKup << endl;
    //cout << endl;    cout << endl;    cout << endl;
    //cout << matSchur << endl;
    //cout << endl;    cout << endl;    cout << endl;

    cout << "matSchur.rows()" << '\t' << "matSchur.cols()" << endl;
    cout << matSchur.rows() << '\t' << matSchur.cols() << endl;

    //solverSchur.preconditioner().setDroptol(1.0e-2);
    //solverSchur.preconditioner().setFillfactor(2);
    //solverSchur.setTolerance(1.0e-8);


    solverSchur.compute(matSchur);


    //velo=veloApplied;
    postProcess();

    timeNow=0.0, loadFactor=0.0;

    //setInitialConditions();
    //loadFactor = 1.0;
    //Time loop
    while( (stepsCompleted < stepsMax ) && (timeNow < timeFinal) )
    {
        //
        if(stepsCompleted < 1000)
        {
          loadFactor = 0.5*(1-cos(PI*stepsCompleted/1000));
          //loadFactor = stepsCompleted/5000.0;
        }
        else
        {
          loadFactor = 1.0;
        }
        //
        loadFactor = 1.0;

        //double  time1 = omp_get_wtime();

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

            //dtCrit = min(dtCrit, elems[ee]->ResidualIncNavStokesSemiImpl(nodeCoords, elemData, timeData, velo, veloPrev, veloDot, veloDotPrev, pres, presPrev, FlocalVelo, FlocalPres, fact) );

            //printVector(FlocalVelo);
            //printVector(FlocalPres);
            //int ii, jj, kk, dd;
            //Assemble the element vector
            //#pragma omp critical
            {
              for(ii=0; ii<npElemVelo; ii++)
              {
                jj = ndim*ii;
                kk = ndim*elemConn[ee][ii];

                for(dd=0; dd<ndim; dd++)
                  rhsVecVeloTemp(kk+dd) += FlocalVelo(jj+dd);
              }
              for(ii=0; ii<npElemPres; ii++)
              {
                rhsVecPresTemp(elemConn[ee][ii])   += FlocalPres(ii);
              }
            }
          } //LoopElem
          dt = dtCrit*CFL;

          amDgammaDt = am/(gamma*dt);

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
          for(ii=0; ii<nDBC_Velo; ++ii)
          {
            n1 = DirichletBCsVelo[ii][0];
            n2 = DirichletBCsVelo[ii][1];

            jj = n1*ndim+n2;

            velo[jj]    = veloApplied[jj] * loadFactor;
          }

          for(ii=0; ii<nDBC_Pres; ++ii)
          {
            jj = DirichletBCsPres[ii][0];

            pres[jj]    = DirichletBCsPres[ii][2] * loadFactor;
          }

          double  fact1 = 1.0/(gamma*dt), fact2 = (gamma-1.0)/gamma;

          veloDot = fact1*(velo-veloPrev) + fact2*veloDotPrev;

          // compute the norms and store the variables
          norm_velo = 0.0;
          norm_velo_diff = 0.0;
          for(ii=0; ii<nsize_velo; ii++)
          {
            norm_velo += velo[ii]*velo[ii];

            fact1 = velo[ii]-veloPrev[ii];
            norm_velo_diff += fact1*fact1;

            // store the variables
            //veloPrev3  = veloPrev2;
            veloPrev2[ii]  = veloPrev[ii];
            veloPrev[ii]  = velo[ii];
            veloDotPrev[ii]  = veloDot[ii];
          }

          norm_pres = 0.0;
          norm_pres_diff = 0.0;
          for(ii=0; ii<nsize_pres; ii++)
          {
            norm_pres += pres[ii]*pres[ii];

            fact1 = pres[ii]-presPrev[ii];
            norm_pres_diff += fact1*fact1;

            //presPrev3  = presPrev2;
            presPrev2[ii]  = presPrev[ii];
            presPrev[ii]     = pres[ii];
          }
        //}

        if( std::isnan(norm_velo) || std::isnan(norm_pres) )
        {
          cerr << " NAN encountered in the solution ... " << endl;
          cerr << " Program has been terminated ... " << endl;
          exit(-1);
        }

        norm_velo = sqrt(norm_velo_diff/norm_velo);
        norm_pres = sqrt(norm_pres_diff/norm_pres);
        if(stepsCompleted == 1)
        {
          norm_velo = 1.0;
          norm_pres = 1.0;
        }

        if( (stepsCompleted%outputFreq == 0) || (norm_velo < conv_tol) )
        {
            cout << " stepsCompleted = " << stepsCompleted << '\t' << " timeNow = " << timeNow << endl;
            cout << " velocity difference norm = " << '\t' << norm_velo << endl;

            postProcess();

            TotalForce.setZero();
            for(ii=0; ii<outputEdges.size(); ii++)
            {
              //cout << ii << '\t' << outputEdges[ii][0] << '\t' << outputEdges[ii][1] << endl;
              //elems[outputEdges[ii][0]]->CalculateForces(outputEdges[ii][1], nodeCoords, elemData, timeData, velo, pres, TotalForce);
              int nn = outputEdges[ii][0];
              TotalForce[0] += rhsVecVeloTemp[nn*2];
              TotalForce[1] += rhsVecVeloTemp[nn*2+1];
            }

            //fout_convdata << timeNow << '\t' << stepsCompleted << '\t' << norm_velo << '\t' << norm_pres ;
            //fout_convdata << '\t' << TotalForce(0) << '\t' << TotalForce(1) << '\t' << TotalForce(2) << endl;
        }

        if(norm_velo < conv_tol)
        {
          cout << " Solution convdataged below the specified tolerance " << endl;
          break;
        }

        stepsCompleted = stepsCompleted + 1;
        timeNow = timeNow + dt;
    } //Time loop


    postProcess();

    //double tend = omp_get_wtime(), duration=tend-tstart;
    double tend = 0.0, duration=tend-tstart;

    cout << " \n \n Total time taken = " << duration << " seconds \n\n" << endl;

    return 0;
}






int femINSmixed::computeElementErrors(int ind)
{
    cout << " Computing errors \n " << endl;

    VectorXd  solnVTK(nNode_Velo*(ndim+1));
    int  n1, n2, dd;
    for(int ii=0; ii<nNode_global; ++ii)
    {
      n1 = ii*ndof;
      n2 = ii*ndim;

      for(dd=0; dd<ndim; dd++)
        solnVTK(n1+dd) = velo(n2+dd);

      solnVTK(n1+ndim) = pres(ii);
    }

    double totalError = 0.0, timeNow;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(int ee=0; ee<nElem_global; ++ee)
      {
        //totalError += elems[ee]->CalculateError(nodeCoords, elemData, timeData, solnVTK, veloDot, pres, timeNow, index);
      }

      totalError = sqrt(totalError);

      if(index == 0)
        printf(" \n\n \t L2 Error in X-velocity = %12.6E \n\n " , totalError);
      else if(index == 1)
        printf(" \n\n \t L2 Error in Y-velocity = %12.6E \n\n " , totalError);
      else if(index == 2)
        printf(" \n\n \t L2 Error in pressure   = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t H1 Error in velocity   = %12.6E \n\n " , totalError);
    }

    return 0;
}




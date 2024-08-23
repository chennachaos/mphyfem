
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


    vector<vector<int> >  NodeDofArray(nNode_Velo, vector<int>(ndim, -1)), LM;
    vector<vector<bool> >  NodeDofType(nNode_Velo, vector<bool>(ndim, false));

    setSpecifiedDOFs_Velocity(NodeDofType);


    int  ee, ii, jj, kk, nn, dof, ind;

    totalDOF_Velo = 0;
    for(ii=0;ii<nNode_global;++ii)
    {
      for(jj=0;jj<ndim;++jj)
      {
        if(!NodeDofType[ii][jj])
        {
          NodeDofArray[ii][jj] = totalDOF_Velo++;
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
          elems[ee]->forAssyVecVelo[ind+jj] = NodeDofArray[kk][jj];
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
            cout << '\t' << NodeDofArray[ii][jj];
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

    int  stepsCompleted=1;
    int  nsize_velo = rhsVecVelo.rows();
    int  nsize_pres = rhsVecPres.rows();
    int  aa, bb, ee, ii, jj, kk, count, row, col, ind, n1, n2, size1, size2;
    bool convergedFlagPrev = false, convergedFlag=false;

    double  fact, fact1, fact2;

    td.resize(100), reacVec.resize(nNode_Velo*ndim+nNode_Pres);

    VectorXd  TotalForce(3);
    VectorXd  FlocalVelo(30), FlocalPres(10);
    ind = npElemVelo*ndim;
    MatrixXd  Kuu(ind,ind), Kup(ind,npElemPres);

    vector<int>  vecTempInt;

    cout << " timeIntegrationScheme = " << timeIntegrationScheme << '\t' << spectralRadius << endl;
 
    SetTimeParametersFluid(timeIntegrationScheme, spectralRadius, myTime.dt, td);

    KimMoinFlowUnsteadyNavierStokes  analy(fluidProperties[0], fluidProperties[1], 1.0);

    double tstart = 0.0;// omp_get_wtime();

    timeNow=0.0;

    //setInitialConditions();

    //Time loop
    while( (myTime.cur <= (timeFinal-EPSILON)) && (stepsCompleted <= stepsMax) )
    {
        // do a time update: reset variables and flags, and prediction step of fields
        timeUpdate();

        printf(" ==================================================================== \n");
        printf(" Time step number     =  %d  \n", stepsCompleted);
        printf(" Time step size       =  %f  \n", myTime.dt);
        printf(" Current time         =  %f  \n", myTime.cur);
        printf(" ==================================================================== \n");

        convergedFlagPrev = convergedFlag;
        convergedFlag = false;

        rhsNormPrev = rhsNorm = -1.0;

        // update Immersed Solids
        //updateImmersedSolid();

        // recompute matrix pattern and update the solver
        if(IB_MOVED)
        {
          prepareImmersedSolids();
          setSolverDataForFullyImplicit();
        }

        cout << " Iteration loop " << endl;

        for(int iter=1; iter<=iterationsMax; iter++)
        {
            cout << " Iteration ... " << iter << endl;

            firstIteration = (iter == 1);

            // Compute the velocity and acceleration and the respective values at n+af, n+am
            updateIterStep();
            //errpetsc = MPI_Barrier(MPI_COMM_WORLD);

            matK *= 0.0;
            rhsVec.setZero();
            reacVec.setZero();

            cout << " Element loop " << endl;

            // loop over elements and compute matrix and residual
            for(ee=0; ee<nElem_global; ee++)
            {
                elems[ee]->StiffnessAndResidualFullyImplicit(nodeCoords, fluidProperties, &td[0], veloPrev, veloPrev2, veloCur, veloDotCur, presCur, Kuu,  Kup, FlocalVelo, FlocalPres);

                //printMatrix(Kuu);
                //printMatrix(Kup);
                //printVector(FlocalVelo);
                //printVector(FlocalPres);

                size1 = elems[ee]->forAssyVecVelo.size();
                size2 = elems[ee]->forAssyVecPres.size();

                //cout << "Applying boundary conditions " << ee << endl;
                if(firstIteration)
                {
                    //printVector(elems[ee]->forAssyVecVelo);
                    //printVector(elems[ee]->forAssyVecPres);
                    //printVector(elems[ee]->globalDOFnums);

                    // applied velocity dof
                    for(ii=0; ii<size1; ii++)
                    {
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
                } // if(firstIteration)
                //cout << "Applying boundary conditions " << ee << endl;

                //printVector(FlocalVelo);
                //printVector(FlocalPres);
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

            //printVector(rhsVec);

            rhsNorm = rhsVec.norm();

            printf(" Iteration : %2d ... RHS norm = %E \n", iter, rhsNorm);
            
            //printVector(rhsVec);

            if(rhsNorm < 1.0e-8)
            {
                //cout << " Solution converged below the specified tolerance " << endl;
                convergedFlag = true;
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

        // if the residual is converged, then save the DOFs vectors
        if( convergedFlag )
        {
            if( stepsCompleted%outputFreq == 0 )
              postProcess();

            writeNodalData();

            writeOutputDataPatches();

            saveSolution();

            myTime.stck();

            stepsCompleted++;

            fileCount++;
        }
        else
        {
            myTime.cut();

            reset();
        }

        cout << endl; cout << endl;

    } //Time loop

    //double tend = omp_get_wtime(), duration=tend-tstart;
    double tend = 0.0, duration=tend-tstart;

    cout << " \n \n Total time taken = " << duration << " seconds \n\n" << endl;

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




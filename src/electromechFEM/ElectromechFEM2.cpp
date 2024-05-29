
#include "SolverPardisoEigen.h"
#include "ElectromechFEM.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "ElementBase.h"
#include "SolutionData.h"
#include "Files.h"
#include "FunctionsProgram.h"
#include <chrono>

extern MyTime myTime;
extern List<TimeFunction> timeFunction;
extern Files files;



void ElectromechFEM::setSolver(int slv, int *parm, bool cIO)
{
    //Eigen::initParallel();
    Eigen::setNbThreads(0);

    int numProc=1;

  //if(IMPLICIT_SOLVER)
  //{
    switch(slv)
    {
        case  4: // SolverEigen ..........................

            //SOLVER_TYPE = SOLVER_TYPE_EIGEN;

            solverEigen = (SolverEigen*) new SolverEigen;

            solverEigen->setAlgorithmType(parm[0]);

            prepareMatrixPattern();

            if(totalDOF > 0)
            {
              if(solverEigen->initialise(0,0,totalDOF) != 0)
                return;
            }
            //solver->setSolverAndParameters();

            solverEigen->printInfo();

        break;

        case  5: // PARDISO(sym) with Eigen
        case  6: // PARDISO(unsym) with Eigen

            //SOLVER_TYPE = SOLVER_TYPE_EIGEN;

            solverEigen = (SolverEigen*) new SolverPardisoEigen;

            if (parm != NULL) numProc = parm[0]; else numProc = 1;

            cout <<  numProc <<  endl;

            prepareMatrixPattern();

            if(slv == 5)
            {
              if (solverEigen->initialise(numProc, PARDISO_STRUCT_SYM, totalDOF) != 0)
                return;
            }
            if(slv == 6)
            {
              if (solverEigen->initialise(numProc, PARDISO_UNSYM, totalDOF) != 0)
                return;
            }

        break;

        default: // invalid slv ...................

             cout << " this solver has not been implemented yet!\n\n";

        break;
    }
  //}
  //else
    //prepareMatrixPattern();


    solverOK = true;

    if(solverEigen != NULL)
      solverEigen->checkIO = cIO;

    if( (tis > 0) )
      setInitialConditions();

    return;
}




int ElectromechFEM::prepareMatrixPattern()
{
    printf("\n     ElectromechFEM::prepareMatrixPattern()  .... STARTED ...\n");

    int  r, c, r1, c1, count=0, count1=0, count2=0, iii, e, ind, nsize;
    int *tt, *tt1, *tt2,  val1, val2, n1, n2, a, b, ll, pp, nnz;
    int  ind1, ind2, ee, ii, jj, kk, e1, e2;

    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    IEN.resize(nElem);

    for(ee=0;ee<nElem;ee++)
    {
      IEN[ee] = elems[ee]->nodeNums;
    }

    NodeType.resize(nNode);
    ID.resize(nNode);
    for(ii=0;ii<nNode;ii++)
    {
      NodeType[ii].resize(ndof);
      ID[ii].resize(ndof);
      for(jj=0;jj<ndof;jj++)
      {
        NodeType[ii][jj] = false;
        ID[ii][jj] = -1;
      }
    }

    for(ii=0;ii<DirichletBCs.size();ii++)
    {
      NodeType[DirichletBCs[ii][0]][DirichletBCs[ii][1]] = true;
    }

    dispDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        if(!NodeType[ii][jj])
          ID[ii][jj] = dispDOF++;
      }
    }

    LM.resize(nElem);
    for(ee=0;ee<nElem;ee++)
    {
      npElem = IEN[ee].size();

      ind = ndof*npElem;
      LM[ee].resize(ind);

      for(ii=0;ii<npElem;ii++)
      {
        ind = ndof*ii;

        kk = IEN[ee][ii];

        for(jj=0;jj<ndof;jj++)
        {
          LM[ee][ind+jj] = ID[kk][jj];
        }
      }
    }

    assyForSoln.resize(dispDOF);
    count = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        if( ID[ii][jj] != -1)
          assyForSoln[count++] = ii*ndof + jj;
      }
    }

    for(ee=0; ee<nElem; ee++)
    {
      elems[ee]->forAssyVec = LM[ee];
      elems[ee]->prepareElemData();
    }

    for(ee=0; ee<nElemFaces; ee++)
    {
      elemsFaces[ee]->forAssyVec.clear();

      for(ii=0; ii<elemsFaces[ee]->nodeNums.size(); ii++)
      {
        ind1 = ii*ndof;
        kk = elemsFaces[ee]->nodeNums[ii];

        for(jj=0;jj<ndof;jj++)
        {
          elemsFaces[ee]->forAssyVec.push_back(ID[kk][jj]);
        }
      }
    }

    prepareDataForPressure();
    prepareDataForElectricPotential();
    prepareDataForTemperature();

    totalDOF = dispDOF + presDOF + epotDOF + tempDOF;

    if(ARCLENGTH)
      totalDOF += 1;

    cout << " nElem    = " << nElem    << endl;
    cout << " nNode    = " << nNode    << endl;
    cout << " npElem   = " << npElem   << endl;
    cout << " ndof     = " << ndof     << endl;
    cout << " dispDOF  = " << dispDOF  << endl;
    cout << " presDOF  = " << presDOF  << endl;
    cout << " epotDOF  = " << epotDOF  << endl;
    cout << " tempDOF  = " << tempDOF  << endl;
    cout << " totalDOF = " << totalDOF << endl;

    printf("\n element DOF values initialised \n\n");
    printf("\n Preparing matrix pattern \n\n");

    if(IMPLICIT_SOLVER)
    {
        forAssyMat.clear();
        forAssyMat.resize(totalDOF);
        for(ii=0; ii<totalDOF; ii++)
        {
          forAssyMat[ii].reserve(500);
        }

        for(ee=0;ee<nElem;ee++)
        {
          tt = &(LM[ee][0]);
          nsize = LM[ee].size();

          for(ii=0;ii<nsize;ii++)
          {
            r = tt[ii];

            if(r != -1)
            {
              for(jj=0;jj<nsize;jj++)
              {
                if(tt[jj] != -1)
                {
                  forAssyMat[r].push_back(tt[jj]);
                }
              }
            }
          }
        }

        if(presDOF > 0)
        {
          cout << "Adding matrix entries for the displacement-pressure coupling " << endl;

          int  *ttP,  *ttU;
          int  sizeP, sizeU, row, col, offset=dispDOF;

          for(ee=0;ee<nElem;ee++)
          {
            sizeP = elems[ee]->forAssyVecPres.size();
            ttP   = &(elems[ee]->forAssyVecPres[0]);

            sizeU = elems[ee]->forAssyVec.size();
            ttU   = &(elems[ee]->forAssyVec[0]);

            for(ii=0; ii<sizeP; ii++)
            {
              row = ttP[ii];

              if(row != -1)
              {
                row += offset;

                // Kpu and Kup
                for(jj=0; jj<sizeU; jj++)
                {
                  col = ttU[jj];
                  if(col != -1)
                  {
                    forAssyMat[row].push_back(col);
                    forAssyMat[col].push_back(row);
                  }
                }

                // Kpp
                for(jj=0; jj<sizeP; jj++)
                {
                  col = ttP[jj];
                  if(col != -1)
                  {
                    forAssyMat[row].push_back(col+offset);
                  }
                }
              }
            }
          }
        }

        if(epotDOF > 0)
        {
          cout << "Adding matrix entries for the displacement-electricpotential coupling " << endl;

          int  *ttF,  *ttU;
          int  sizeF, sizeU, row, col, offset=dispDOF+presDOF;

          for(ee=0;ee<nElem;ee++)
          {
            sizeF = elems[ee]->forAssyVecEpot.size();
            ttF   = &(elems[ee]->forAssyVecEpot[0]);

            sizeU = elems[ee]->forAssyVec.size();
            ttU   = &(elems[ee]->forAssyVec[0]);

            for(ii=0; ii<sizeF; ii++)
            {
              row = ttF[ii];

              if(row != -1)
              {
                row += offset;

                // Kfu and Kuf
                for(jj=0; jj<sizeU; jj++)
                {
                  col = ttU[jj];
                  if(col != -1)
                  {
                    forAssyMat[row].push_back(col);
                    forAssyMat[col].push_back(row);
                  }
                }

                // Kff
                for(jj=0; jj<sizeF; jj++)
                {
                  col = ttF[jj];
                  if(col != -1)
                  {
                    forAssyMat[row].push_back(col+offset);
                  }
                }
              }
            }
          }
        }

        if(tempDOF > 0)
        {
          cout << "Adding matrix entries for the displacement-temperature coupling " << endl;

          int  *ttT,  *ttU;
          int  sizeT, sizeU, row, col, offset=dispDOF+presDOF+epotDOF;

          for(ee=0;ee<nElem;ee++)
          {
            sizeT = elems[ee]->forAssyVecTemp.size();
            ttT   = &(elems[ee]->forAssyVecTemp[0]);

            sizeU = elems[ee]->forAssyVec.size();
            ttU   = &(elems[ee]->forAssyVec[0]);

            for(ii=0; ii<sizeT; ii++)
            {
              row = ttT[ii];

              if(row != -1)
              {
                row += offset;

                // Ktu and Kut
                for(jj=0; jj<sizeU; jj++)
                {
                  col = ttU[jj];
                  if(col != -1)
                  {
                    forAssyMat[row].push_back(col);
                    forAssyMat[col].push_back(row);
                  }
                }

                // Ktt
                for(jj=0; jj<sizeT; jj++)
                {
                  col = ttT[jj];
                  if(col != -1)
                  {
                    forAssyMat[row].push_back(col+offset);
                  }
                }
              }
            }
          }
        }


        if(ARCLENGTH)
        {
          jj = totalDOF-1;

          for(ii=0; ii<jj; ii++)
          {
            forAssyMat[ii].push_back(jj);
            forAssyMat[jj].push_back(ii);
          }
          forAssyMat[jj].push_back(jj);

          solnIncrStep.resize(totalDOF);
          solnIncrIter.resize(totalDOF);
          rhsVec.resize(totalDOF);

          solnIncrStep.setZero();
          solnIncrIter.setZero();
          rhsVec.setZero();

          lambIncrStep = 0.0;
        }

        printf("\n Preparing matrix pattern DONE \n\n");

        VectorXi  nnzVec(totalDOF);

        nnz = 0;
        for(ii=0;ii<totalDOF;ii++)
        {
          findUnique(forAssyMat[ii]);
          nnzVec[ii] = forAssyMat[ii].size();
          nnz += nnzVec[ii];
        }

        cout << " nnz " << nnz << endl;


        cout << " Eigen based solver " << totalDOF << endl;
        if(totalDOF > 0)
        {
          solverEigen->mtx.setZero();

          solverEigen->mtx.resize(totalDOF, totalDOF);
          solverEigen->mtx.reserve(nnz);
          solverEigen->mtx.reserve(nnzVec);

          for(ii=0; ii<totalDOF; ii++)
          {
            for(jj=0;jj<forAssyMat[ii].size();jj++)
            {
              solverEigen->mtx.coeffRef(ii, forAssyMat[ii][jj]) = 0.0;
            }
          }

          solverEigen->mtx.makeCompressed();
        }
    }
    else
      setSolverDataForSemiImplicit();


    solverEigen->currentStatus = PATTERN_OK;

    // remove data objects
    nodePosData.clear();
    elemConn.clear();
    IEN.clear();
    LM.clear();
    ID.clear();
    forAssyMat.clear();

    printf("\n     ElectromechFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return 1;
}



// set the off-diagonal terms for the solver
void ElectromechFEM::setSolverDataForSemiImplicit()
{
    int  idd = SolnData.ElemProp[0]->id;
    int  ee, ii, jj, size1, size2, row, col;
    vector<int>  vecIntTemp(10);

    // inverse of Kuu matrix
    solverEigen->matKuuInv.setZero();
    solverEigen->matKuuInv.resize(dispDOF, dispDOF);
    solverEigen->matKuuInv.reserve(dispDOF);

    for(ii=0; ii<dispDOF; ii++)
    {
      solverEigen->matKuuInv.coeffRef(ii,ii) = 0.0;
    }

    VectorXi  nnzVec(dispDOF);

    // Kup matrix
    solverEigen->matKup.setZero();
    solverEigen->matKup.resize(dispDOF, presDOF);
    solverEigen->matKup.reserve(ceil(presDOF*dispDOF*0.1));
    jj = ceil(presDOF*0.1);
    cout << " jj = " << jj << endl;

    for(ii=0; ii<dispDOF; ii++)
      nnzVec(ii) = jj;
    solverEigen->matKup.reserve(nnzVec);


    // Kpu matrix
    solverEigen->matKpu.setZero();
    solverEigen->matKpu.resize(presDOF, dispDOF);
    solverEigen->matKpu.reserve(ceil(presDOF*dispDOF*0.1));
    solverEigen->matKpu.reserve(nnzVec);

    // Kpp matrix
    solverEigen->matKpp.setZero();
    solverEigen->matKpp.resize(presDOF, presDOF);
    solverEigen->matKpp.reserve(ceil(presDOF*presDOF*0.1));
    solverEigen->matKpp.reserve(nnzVec);


    for(ee=0; ee<nElem; ee++)
    {
        size1 = elems[ee]->forAssyVec.size();
        size2 = elems[ee]->forAssyVecPres.size();

        for(ii=0; ii<size2; ii++)
        {
          row = elems[ee]->forAssyVecPres[ii];

          if(row != -1)
          {
            for(jj=0; jj<size1; jj++)
            {
              col = elems[ee]->forAssyVec[jj];

              if(col != -1)
              {
                solverEigen->matKup.coeffRef(col, row) = 0.0;
                solverEigen->matKpu.coeffRef(row, col) = 0.0;
              }
            }

            for(jj=0; jj<size2; jj++)
            {
              col = elems[ee]->forAssyVecPres[jj];

              if(col != -1)
              {
                solverEigen->matKpp.coeffRef(row, col) = 0.0;
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
    } //for(ee=0;)

    solverEigen->matKuuInv.makeCompressed();
    solverEigen->matKup.makeCompressed();
    solverEigen->matKpu.makeCompressed();
    solverEigen->matKpp.makeCompressed();

    // Sparse matrix for the electric field

    cout << "Preparing sparse matrix for the electric field problem" << endl;

    rhsElecField.resize(epotDOF);
    spmtxElecField.setZero();
    spmtxElecField.resize(epotDOF, epotDOF);
    spmtxElecField.reserve(ceil(epotDOF*epotDOF*0.1));
    jj = ceil(epotDOF*0.2);
    cout << " jj = " << jj << endl;

    for(ii=0; ii<epotDOF; ii++)
      nnzVec(ii) = jj;
    spmtxElecField.reserve(nnzVec);


    int  *ttF, *ttU, sizeF;

    for(ee=0;ee<nElem;ee++)
    {
        sizeF = elems[ee]->forAssyVecEpot.size();
        ttF   = &(elems[ee]->forAssyVecEpot[0]);

        for(ii=0; ii<sizeF; ii++)
        {
            row = ttF[ii];

            if(row != -1)
            {
              for(jj=0; jj<sizeF; jj++)
              {
                col = ttF[jj];
                if(col != -1)
                {
                  spmtxElecField.coeffRef(row, col) = 0.0;
                }
              }
            }
        }
    }

    return;
}




int ElectromechFEM::solveStep(int niter)
{
    if(totalDOF == 0)
      return 0;

    for(int iter=0; iter<niter; iter++)
    {
      calcStiffnessAndResidual();

      if( converged() )
        break;

      factoriseSolveAndUpdate();

      updateIterStep();
    }

    if( !converged() )
      return -1;

    return 0;
}


/*
int ElectromechFEM::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
    if(totalDOF == 0)
      return 0;


    if(ARCLENGTH)
      calcStiffnessAndResidualAL();
    else
      calcStiffnessAndResidualNR();

    return 0;
}
*/

/*
int ElectromechFEM::calcStiffnessAndResidual(int arclengthflag, bool zeroMtx, bool zeroRes)
{
    //cout << "     ElectromechFEM: generating coefficient Matrices ...\n\n";

    int  ee, ii, jj, nsize, size1, size2, row;

    if(solverEigen == NULL)
    {
      COUT << "You need to select a solver first!\n\n";
      return 1;
    }

    arclengthflag = false;
    //cout << "arclengthflag = " << arclengthflag << endl;

    if(firstIter)
    {
      rNorm = -1.0;
    }

    MatrixXd  Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp;
    VectorXd  FlocalU, FlocalF, FlocalP;
    //MatrixXd  Kuu(30, 30), Kuf(30, 10), Kfu(10,30), Kff(10,10), Kup, Kpu, Kpp;
    //VectorXd  FlocalU(30), FlocalF(10), FlocalP(10);


    SolnData.reac.setZero();

    solverEigen->zeroMtx();

    //printVector(solver->rhsVec);

    if(MIXED_ELEMENT)
    {
      int offset2 = dispDOF+epotDOF;
      // loop over all the elements
      for(ee=0;ee<nElem;ee++)
      {
        //cout << "       elem... : " << (ee+1) << endl;

        elems[ee]->calcStiffnessAndResidualMixed(Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp, FlocalU, FlocalF, FlocalP, firstIter);
        elems[ee]->calcLoadVector(FlocalU);

        if(firstIter)
          elems[ee]->applyDirichletBCs(1, Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp, FlocalU, FlocalF, FlocalP);

        solverEigen->assembleMatrixAndVectorUFP(dispDOF, offset2, elems[ee]->forAssyVec, elems[ee]->forAssyVecEpot, elems[ee]->forAssyVecPres, Kuu, Kuf, Kfu,  Kff, Kup, Kpu, Kpp, FlocalU, FlocalF, FlocalP);
      }
    }
    else
    {
      // loop over all the elements
      for(ee=0;ee<nElem;ee++)
      {
        //cout << "       elem... : " << (ee+1) << endl;
        elems[ee]->calcStiffnessAndResidual(Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, firstIter);
        elems[ee]->calcLoadVector(FlocalU);

        if(!arclengthflag)
        {
          if(firstIter)
            elems[ee]->applyDirichletBCs(1, Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF);
        }

        solverEigen->assembleMatrixAndVectorUF(dispDOF, elems[ee]->forAssyVec, elems[ee]->forAssyVecEpot, Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF);

        if(firstIter && arclengthflag)
        {
          FlocalU.setZero();
          FlocalF.setZero();

          elems[ee]->applyDirichletBCs(1, Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF);

          int size1 = elems[ee]->forAssyVec.size();
          int size2 = elems[ee]->forAssyVecEpot.size();

          for(ii=0; ii<size1; ii++)
          {
            row = elems[ee]->forAssyVec[ii];
            if( row != -1 )
            {
              rhsVec[row] += FlocalU(ii);
            }
          }

          for(ii=0; ii<size2; ii++)
          {
            row = elems[ee]->forAssyVecEpot[ii];
            if( row != -1 )
            {
              rhsVec[dispDOF+row] += FlocalF(ii);
            }
          }
        }
      }
    }

    //cout << " BBBBBBBBBBBBB " << endl;
    //cout << solver->mtx << endl;
    //printf("\n solverEigen->rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec);

    addExternalForces();

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec);

    //printf("\n solverEigen->rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    if(firstIter)
    {
      int ii, nn, dof, ind;

      for(ii=0; ii<DirichletBCs.size(); ii++)
      {
        nn  = (int) (DirichletBCs[ii][0]);
        dof = (int) (DirichletBCs[ii][1]);

        ind = nn*ndof+dof;

        //cout << nn << '\t' << ind << '\t' << SolnData.var1applied[ind] << endl;

        SolnData.var1[ind] += SolnData.var1applied[ind] ;
      }
      for(ii=0; ii<DirichletBCs_Epot.size(); ii++)
      {
        nn  = (int) (DirichletBCs_Epot[ii][0]);

        //cout << nn << '\t' << ind << '\t' << SolnData.var3applied[nn] << endl;

        SolnData.var3[nn] += SolnData.var3applied[nn] ;
      }
    }

    //cout << solverEigen->mtx << endl;
    if(ARCLENGTH)
    {
    if(arclengthflag)
    {
      double  psi = 1.0,  psi2 = psi*psi;
      double rhsVecMagn = 0.0, solnIncrStepMagn = 0.0, lfact = lamb;
      dl = myTime.dt;
      //dl = 1.0;

      int ind = totalDOF-1;

      for(ii = 0; ii<ind; ii++)
      {
        solverEigen->mtx.coeffRef(ind, ii)  +=  2.0*solnIncrStep(ii);
        solverEigen->mtx.coeffRef(ii,  ind) -=  rhsVec(ii);

        solverEigen->rhsVec(ii) -=  lamb*rhsVec(ii);

        rhsVecMagn +=  rhsVec(ii)*rhsVec(ii);

        solnIncrStepMagn +=  solnIncrStep(ii)*solnIncrStep(ii);
      }
      solverEigen->mtx.coeffRef(ind, ind) += (2.0*psi2)*lambIncrStep*rhsVecMagn;

      solverEigen->rhsVec(ind) -=  (solnIncrStepMagn + psi2*lambIncrStep*lambIncrStep*rhsVecMagn - dl*dl);
    }
    else
    {
      int ind = totalDOF-1;

      solverEigen->mtx.coeffRef(ind, ind) += 1.0;
      solverEigen->rhsVec(ind) =  0.0;
    }
    }



    firstIter = false;
    rNormPrev = rNorm;

    rNorm = solverEigen->rhsVec.norm();
    solverEigen->currentStatus = ASSEMBLY_OK;

    COUT << "ElectromechFEM"; printf("  %11.4e\n",rNorm);

    return 0;
}
*/



int ElectromechFEM::calcStiffnessAndResidual(int arclengthflag, bool zeroMtx, bool zeroRes)
{
    //cout << "     ElectromechFEM: generating coefficient Matrices ...\n\n";

    int  ee, ii, jj, nsize, size1, size2, row;

    if(solverEigen == NULL)
    {
      cout << "You need to select a solver first!\n\n" << endl;
      return 1;
    }

    if(firstIter)
    {
      rNorm = -1.0;
    }

    MatrixXd  Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp;
    VectorXd  FlocalU, FlocalF, FlocalP;
    //MatrixXd  Kuu(30, 30), Kuf(30, 10), Kfu(10,30), Kff(10,10), Kup, Kpu, Kpp;
    //VectorXd  FlocalU(30), FlocalF(10), FlocalP(10);


    SolnData.reac.setZero();

    solverEigen->zeroMtx();

    //printVector(solver->rhsVec);

    if(MIXED_ELEMENT)
    {
      int offset2 = dispDOF+presDOF;
      // loop over all the elements
      for(ee=0;ee<nElem;ee++)
      {
        //cout << "       elem... : " << (ee+1) << endl;

        elems[ee]->calcStiffnessAndResidualMixed(Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp, FlocalU, FlocalF, FlocalP, firstIter);
        elems[ee]->calcLoadVector(FlocalU);

        if(firstIter)
          elems[ee]->applyDirichletBCs3field(1, Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp, FlocalU, FlocalF, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecEpot, elems[ee]->forAssyVecPres, SolnData.var1applied, SolnData.var3applied, SolnData.var2applied);

        solverEigen->assembleMatrixAndVector3field(dispDOF, offset2, Kuu, Kup, Kpu, Kpp, Kuf, Kfu, Kff, FlocalU, FlocalP, FlocalF, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres, elems[ee]->forAssyVecEpot);
      }
    }
    else
    {
      // loop over all the elements
      if( tempDOF == 0 )
      {
        for(ee=0;ee<nElem;ee++)
        {
          //cout << "       elem... : " << (ee+1) << endl;
          elems[ee]->calcStiffnessAndResidual(Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, firstIter);
          elems[ee]->calcLoadVector(FlocalU);

          if(firstIter)
            elems[ee]->applyDirichletBCs2field(1, Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, elems[ee]->forAssyVec, elems[ee]->forAssyVecEpot, SolnData.var1applied, SolnData.var3applied);

          solverEigen->assembleMatrixAndVector2field(dispDOF, Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, elems[ee]->forAssyVec, elems[ee]->forAssyVecEpot);
        }
      }
      if( epotDOF == 0 )
      {
        for(ee=0;ee<nElem;ee++)
        {
          //cout << "       elem... : " << (ee+1) << endl;
          elems[ee]->calcStiffnessAndResidual(Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, firstIter);
          elems[ee]->calcLoadVector(FlocalU);

          if(firstIter)
            elems[ee]->applyDirichletBCs2field(1, Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, elems[ee]->forAssyVec, elems[ee]->forAssyVecTemp, SolnData.var1applied, SolnData.var4applied);

          solverEigen->assembleMatrixAndVector2field(dispDOF, Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, elems[ee]->forAssyVec, elems[ee]->forAssyVecTemp);
        }
      }
    }

    //cout << " BBBBBBBBBBBBB " << endl;
    //cout << solver->mtx << endl;
    //printf("\n solverEigen->rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec);

    addExternalForces();

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec);

    //printf("\n solverEigen->rhsVec norm = %12.6E \n", solverEigen->rhsVec.norm());

    if(SolnData.var1Incr.rows() != (nNode*ndim)) SolnData.var1Incr.resize(nNode*ndim);
    SolnData.var1Incr.setZero();

    if(firstIter)
    {
      int ii, nn, dof, ind;

      for(ii=0; ii<DirichletBCs.size(); ii++)
      {
        nn  = (int) (DirichletBCs[ii][0]);
        dof = (int) (DirichletBCs[ii][1]);

        ind = nn*ndof+dof;

        //SolnData.var1[ind] += SolnData.var1applied[ind] ;

        SolnData.var1Incr[ind] += SolnData.var1applied[ind] ;
      }

      for(ii=0; ii<DirichletBCs_Pres.size(); ii++)
      {
        nn  = (int) (DirichletBCs_Pres[ii][0]);

        SolnData.var2[nn] += SolnData.var2applied[nn] ;
      }

      for(ii=0; ii<DirichletBCs_Epot.size(); ii++)
      {
        nn  = (int) (DirichletBCs_Epot[ii][0]);

        SolnData.var3[nn] += SolnData.var3applied[nn] ;
      }
    }

    firstIter = false;
    rNormPrev = rNorm;

    rNorm = solverEigen->rhsVec.norm();
    solverEigen->currentStatus = ASSEMBLY_OK;

    printf(" ElectromechFEM %11.4e \n",rNorm);

    return 0;
}



int ElectromechFEM::factoriseSolveAndUpdate()
{
    //cout << " ElectromechFEM::factoriseSolveAndUpdate " << endl;

    time_t tstart, tend;

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec); printf("\n\n\n");
    //tstart = time(0);

    solverEigen->factoriseAndSolve();

    //tend = time(0);
    //printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    //printVector(solverEigen->soln);

    // update solution vector
    for(int ii=0; ii<dispDOF; ii++)
    {
      SolnData.var1Incr[assyForSoln[ii]] = solverEigen->soln[ii];
    }
    SolnData.var1 += SolnData.var1Incr;


    if(MIXED_ELEMENT)
    {
      for(int ii=0; ii<presDOF; ii++)
      {
        SolnData.var2[assyForSolnPres[ii]] += solverEigen->soln[dispDOF+ii];
      }
    }

    int epot_offset = dispDOF + presDOF;
    for(int ii=0; ii<epotDOF; ii++)
    {
      SolnData.var3[assyForSolnEpot[ii]] += solverEigen->soln[epot_offset+ii];
    }


    if(ARCLENGTH)
    {
        solnIncrIter = solverEigen->soln;
        solnIncrStep += solnIncrIter;

        //lambIncrIter = solverEigen->soln[totalDOF-1];
        lambIncrStep += lambIncrIter;
        lamb         += lambIncrIter;
        //cout << "lambIncrIter = " << lambIncrIter << '\t' <<  lamb <<  endl;
    }


    // solve pressure variable in the mixed formulation
    // only for the constant pressure elements, TRIA6/1, TET10/1

    int  idd = SolnData.ElemProp[0]->id;
    if( (idd == 1005) || (idd == 1010) || (idd == 1011) || (idd == 1013) || (idd == 1045) || (idd == 1050) || (idd == 1051) || (idd == 1053) )
    {
      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->solveForPressure();
      }
    }

    //cout << " ElectromechFEM::factoriseSolveAndUpdate ... DONE " << endl;

    return 0;
}





void ElectromechFEM::addExternalForces()
{
    //cout <<  " ElectromechFEM::addExternalForces() " <<  endl;

    int  ee, nn, dof, ii, jj, ind, ind1, ind2;
    double specVal=0.0, loadfactor=0.0, fact1=0.0;

    MatrixXd  Klocal(30, 30);
    VectorXd  vecTemp(nNode*ndof), vecTempEpot(nNode), Flocal;
    vecTemp.setZero();
    vecTempEpot.setZero();

    // forces due to face/edge loads
    for(ee=0; ee<nElemFaces; ee++)
    {
      //cout << " jjjjjjjjjj " << '\t' << ee << endl;

      //elemsFaces[ee]->calcLoadVector(Flocal);
      elemsFaces[ee]->calcStiffnessAndResidual(Klocal, Flocal);

      //printVector(elemsFaces[ee]->forAssyVec);
      //printMatrix(Klocal);
      //printVector(Flocal);

      solverEigen->assembleMatrixAndVectorElecField(elemsFaces[ee]->forAssyVec, Klocal, Flocal);

      /*
      for(ii=0; ii<elemsFaces[ee]->nodeNums.size(); ii++)
      {
        ind1 = elemsFaces[ee]->nodeNums[ii]*ndof;
        ind2 = ii*ndof;
        for(jj=0; jj<ndof; jj++)
        {
          vecTemp(ind1+jj) += Flocal(ind2+jj);
        }
      }
      */
    }

    // specified nodal forces
    for(ii=0;ii<nodeForcesData.size();ii++)
    {
      nn  = (int) (nodeForcesData[ii][0] - 1);
      dof = (int) (nodeForcesData[ii][1] - 1);
      specVal = nodeForcesData[ii][2];

      //cout << nn << '\t' << dof << '\t' << specVal << endl;

      if(dof < ndim)
      {
        ind = nn*ndof+dof;

        vecTemp[ind] += specVal;
      }
      else
      {
        vecTempEpot[nn] += specVal;
      }
    }
    //printVector(vecTemp);

    //loadfactor = 1.0;
    loadfactor = timeFunction[0].prop;
    //cout << "fact=" << myTime.cur << '\t' << fact << endl;

    //if ( timeFunction[0].prop < 1.0)
      //loadfactor = timeFunction[0].prop;
    //else
      //loadfactor = 0.0;


    if(IMPLICIT_SOLVER)
    {
        for(ii=0; ii<dispDOF; ii++)
          solverEigen->rhsVec[ii] += (loadfactor*vecTemp[assyForSoln[ii]]);

        int epot_offset = dispDOF + presDOF;
        for(ii=0; ii<epotDOF; ii++)
          solverEigen->rhsVec[epot_offset+ii] += (loadfactor*vecTempEpot[assyForSolnEpot[ii]]);
    }
    else
    {
      //ForceVectorExternal.setZero();
      ForceVectorExternal = vecTemp;
      //printVector(ForceVectorExternal);
      //for(ii=0;ii<totalDOF;ii++)
        //rhsVec[ii] += (fact*vecTemp[assyForSoln[ii]]);
    }

    vecTemp.resize(0);

    //cout << "ElectromechFEM::addExternalForces() ... ENDED" << endl;

    return;
}




void  ElectromechFEM::computeElementErrors(int index)
{
    int  ii, ee, count=0, dd, domTemp;

    cout << " index = " << index << endl;

    totalError = 0.0;

    if(index < 4) // L2 or H1 norm based errors
    {
      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->calcError(index);

        totalError += elems[ee]->getError();
      }

      totalError = sqrt(totalError);

      if(index == 0)
        printf(" \n\n \t Displacement Error = %12.6E \n\n " , totalError);
      else if(index == 1)
        printf(" \n\n \t Stress Error = %12.6E \n\n " , totalError);
      else if(index == 2)
        printf(" \n\n \t Pressure Error = %12.6E \n\n " , totalError);
      else
        printf(" \n\n \t H1 Error = %12.6E \n\n " , totalError);
    }
    else if(index == 10) // total energy
    {
      VectorXd  energyElem(3), energyGlobal(3);
      energyGlobal.setZero();

      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->computeEnergy(0, 0, energyElem);

        energyGlobal += energyElem;
      }
      energyGlobal[2] = energyGlobal[0] + energyGlobal[1];

      char        tmp[200];
      MyString    tmpStr;

      sprintf(tmp," \t %12.6E \t %12.6E \t %12.6E", energyGlobal[0], energyGlobal[1], energyGlobal[2]);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
    }
    else if(index == 11) // total momentum
    {
      VectorXd  momElem(6), momGlobal(6);
      momGlobal.setZero();

      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->computeMomentum(0, 0, momElem);

        momGlobal += momElem;
      }

      char        tmp[200];
      MyString    tmpStr;

      sprintf(tmp," \t %12.6E \t %12.6E \t %12.6E", momGlobal[0], momGlobal[1], momGlobal[5]);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);

    }
    else
    {
      cout << " Compute element errors not defined for thie index value " << endl;
    }
    //printf(" \n\n \t totalError = %12.6E \n\n " , totalError);

    return;
}






int ElectromechFEM::solveExplicitStep(int* inputdata)
{
    // just in case this function is called by mistake while using the implicit solver
    if(IMPLICIT_SOLVER)
      return 0;

    auto  time1 = chrono::steady_clock::now();


    int  bb, ee, ii, jj, nn, n1, n2, n3, row, cc, size1, size2;
    int  dof, ind, ind1, ind2, iter, niter=10;
    int  resln[3]={1,1,1};
    double timeFact=0.0, fact=0.0, fact1=0.0, vec[10];

    int stepsMax = inputdata[0];
    int outputfreq = inputdata[1];
    int stepsCompleted = 1;
    double  dt   = myTime.dt, dtCrit;
    double  timeNow=0.0;
    double  timeFinal = stepsMax*dt;
    double  af = 1.0;
    double  am = 1.0;
    double  gamma = 0.5+am;
    //double  beta  = 1.0; // 2nd-order accurate
    double  beta  = 28.0/27.0; // 2nd-order accurate
    //double  beta  = am + 1.0/12.0;                         // 3rd-order accurate
    double  alpha1 = 0.0;
    double  amDbetaDTT;

    VectorXd  Flocal1(30), Flocal2(10), rhsVecFa(nNode*ndof), matK1(nNode*ndof);
    //MatrixXd  Kup(30,10), Kpu(10,30), Kpp(4,4);
    MatrixXd  Kup, Kpu, Kpp;
    //SparseMatrixXd  matK1inv(totalDOF, totalDOF), matKpu(nElem, totalDOF);


    // Store variables
    SolnData.var1Prev        =  SolnData.var1;
    SolnData.var1DotPrev     =  SolnData.var1Dot;
    SolnData.var1DotDotPrev  =  SolnData.var1DotDot;
    SolnData.var3Prev        =  SolnData.var3;
    SolnData.var2Prev        =  SolnData.var2;

    cout << "Semi-Implicit scheme - Matrices and Vectors set " << endl;

    filecount=1;

    // calculate the mass matrix
    calcMassMatrixForExplicitDynamics();


    amDbetaDTT = am/(beta*dt*dt);

    matK1 = amDbetaDTT*MassVector;

    solverEigen->matKuuInv *= 0.0;
    for(ii=0; ii<dispDOF; ii++)
    {
      jj = assyForSoln[ii];

      solverEigen->matKuuInv.coeffRef(ii,ii) = 1.0/matK1(jj);
    }

    //printVector(MassVector);

    addExternalForces();

    //Time loop
    timeFinal = 100.0+1.0e-10;
    while( (stepsCompleted <= stepsMax ) && (timeNow <= timeFinal) )
    {
        /////////////////////////////
        //
        // Part 1: Mechanical problem
        //
        /////////////////////////////

        /*
        dtCrit=1.0e10;
        for(ee=0; ee<nElem; ee++)
          dtCrit = min(dtCrit, elems[ee]->calcCriticalTimeStep(true));

        dt = dtCrit*0.75;
        //cout << " dt = " << dt << endl;

        amDbetaDTT = am/(beta*dt*dt);

        matK1 = amDbetaDTT*MassVector;

        solverEigen->matKuuInv *= 0.0;
        for(ii=0; ii<dispDOF; ii++)
        {
          jj = assyForSoln[ii];

          solverEigen->matKuuInv.coeffRef(ii,ii) = 1.0/matK1(jj);
        }
        //matK1inv.makeCompressed();
        */

        timeNow += dt;

        myTime.dt = dt;
        myTime.cur = timeNow;

        timeFunction[0].update();

        //if(intVarFlag)
          //copyElemInternalVariables();

        //cout << " aaaaaaaaaaa " << endl;

        assignBoundaryConditions();

      //for(iter=0; iter<niter; ++iter)
      //{
        // loop over all the elements and calculate the residual vector
        //////////////////////////////////////////
        solverEigen->rhsVec.setZero();
        solverEigen->rhsVec2.setZero();

        solverEigen->matKup *= 0.0;
        solverEigen->matKpu *= 0.0;
        solverEigen->matKpp *= 0.0;

        for(ee=0;ee<nElem;ee++)
        {
          elems[ee]->calcStiffnessAndResidualMixed2(Kup, Kpu, Kpp, Flocal1, Flocal2);
          //elems[ee]->calcLoadVector(Flocal1);
          //printMatrix(Kup);
          //printMatrix(Kpu);
          //printMatrix(Kpp);

          elems[ee]->applyDirichletBCs2field(0, Kup, Kup, Kpu, Kpp, Flocal1, Flocal2, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres, SolnData.var1applied, SolnData.var2applied);

          solverEigen->assembleMatrixAndVectorUPexpl(0, Kup, Kpu, Kpp, Flocal1, Flocal2, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres);
        }

        // add specified nodal force
        addExternalForces();

        // compute Fa, due to the acceleration term
        timeFact = timeFunction[0].prop;
        for(ii=0; ii<dispDOF; ii++)
        {
          jj = assyForSoln[ii];

          fact  = matK1(jj)*(SolnData.var1(jj) - SolnData.var1Prev(jj));
          fact -= ((1.0/beta/dt)*MassVector(jj))*SolnData.var1DotPrev(jj) ;
          fact += ((1.0-am*0.5/beta)*MassVector(jj))*SolnData.var1DotDotPrev(jj);

          solverEigen->rhsVec(ii) += (ForceVectorExternal(jj) - fact);
        }

        //rNorm = solverEigen->rhsVec.norm();
        //printf("Iter= %d \t Norm = %12.6E \n", (iter+1), rNorm );
        //if( rNorm < 1.0e-6)
          //break;

        auto  time1 = chrono::steady_clock::now();

        //solverEigen->solverSchurExplicitDirect(10);
        solverEigen->solverSchurExplicitDirect(outputfreq);
        //printVector(solnPres);

        auto  time2 = chrono::steady_clock::now();
        auto  duration = chrono::duration_cast<chrono::microseconds>(time2-time1).count();
        //cout << " Time (microseconds)   = " << duration << endl;

        // update the solution
        for(ii=0; ii<dispDOF; ii++)
        {
          SolnData.var1(assyForSoln[ii]) += solverEigen->var1(ii);
        }
        for(ii=0; ii<presDOF; ii++)
        {
          SolnData.var2[assyForSolnPres[ii]] += solverEigen->var2[ii];
        }

        // add Dirichlet BCs
        for(ii=0; ii<DirichletBCs.size(); ii++)
        {
          nn  = (int) (DirichletBCs[ii][0]);
          dof = (int) (DirichletBCs[ii][1]);

          ind = nn*ndof+dof;

          SolnData.var1(ind) += SolnData.var1applied(ind);
        }

        // check if the deformation is too large
        if( (SolnData.var1.maxCoeff() > 1.0e10) ||  (abs(SolnData.var1.minCoeff()) > 1.0e10) )
        {
          cerr << "displacement too high or NAN in explicit simulation..." << endl;
          cerr << "Simulation aborted...\n\n\n" << endl;
          exit(1);
        }

        // update acceleration and velocity
        SolnData.var1DotDot = (SolnData.var1 - SolnData.var1Prev - dt*SolnData.var1DotPrev)/(beta*dt*dt) - (0.5/beta-1.0)*SolnData.var1DotDotPrev;
        SolnData.var1Dot    = SolnData.var1DotPrev + dt*( gamma*SolnData.var1DotDot + (1.0-gamma)*SolnData.var1DotDotPrev );

      //} // for(iter=0; iter<niter; ++iter)

        // update the coordinates
        updateGeometry();

        /////////////////////////////
        //
        // Part 2: Electric Field problem
        //
        /////////////////////////////

        bool dummyFlag = ((stepsCompleted-1)%outputfreq == 0);
        //bool dummyFlag = ((stepsCompleted-1)%10 == 0);

        solveElectricField(dummyFlag, 1);
        //solveElectricField(true, 1);

        // write the data for every "outputfreq" time steps
        if( (stepsCompleted%outputfreq == 0) )
        {
          printf("\n\n\n stepsCompleted = %5d ;  timeNow = %12.8f \n", stepsCompleted, timeNow);
          // write the output data
          writeNodalData();

          postProcess(4, 10, 1, false, 0, 0, resln);
          //cout <<  " aaaaaaaaaa " <<  endl;
          filecount++;
        }

        //store the variables
        SolnData.var1Prev        =  SolnData.var1;
        SolnData.var1DotPrev     =  SolnData.var1Dot;
        SolnData.var1DotDotPrev  =  SolnData.var1DotDot;
        SolnData.var3Prev        =  SolnData.var3;
        SolnData.var2Prev        =  SolnData.var2;

        ++stepsCompleted;
    }

    auto  time2 = chrono::steady_clock::now();
    auto  duration = chrono::duration_cast<chrono::milliseconds>(time2-time1).count();
    cout << " SemiImplicit Scheme took  " << duration << " milliseconds " <<  endl;

    return 0;
}



int ElectromechFEM::solveElectricField(bool update_matrix, int niter)
{
    //cout << "     ElectromechFEM: solveElectricField ...\n\n";
    //cout <<  "update_matrix = " <<  update_matrix <<  endl;

    int  ee, ii, jj, size1, r, c, iter, nn;
    double rNorm = 0.0, tol = 1.0e-6;
    MatrixXd  Kff;
    VectorXd  FlocalF, solnTemp(epotDOF);
    vector<int>  vecTemp;


    for(iter=0; iter<niter; iter++)
    {
        if(update_matrix)
          spmtxElecField *= 0.0;
        rhsElecField.setZero();

        for(ee=0;ee<nElem;ee++)
        {
            //cout << "       elem... : " << (ee+1) << endl;

            elems[ee]->calcStiffnessAndResidualElecField(Kff, FlocalF);
            //elems[ee]->calcLoadVector(FlocalU);

            if(iter == 0)
              elems[ee]->applyDirichletBCsElecField(Kff, FlocalF);

            //printMatrix(Kff);
            //printVector(FlocalF);

            vecTemp = elems[ee]->forAssyVecEpot;
            size1 = vecTemp.size();

            if(update_matrix)
            {
              for(ii=0; ii<size1; ii++)
              {
                r = vecTemp[ii];
                if( r != -1 )
                {
                  rhsElecField[r] += FlocalF(ii);

                  for(jj=0; jj<size1; jj++)
                  {
                    c = vecTemp[jj];
                    if( c != -1 )
                      spmtxElecField.coeffRef(r, c) += Kff(ii,jj);
                  }
                }
              }
            }
            else
            {
              for(ii=0; ii<size1; ii++)
              {
                r = vecTemp[ii];
                if( r != -1 )
                {
                  rhsElecField[r] += FlocalF(ii);
                }
              }
            }
        }

        //cout << " RHS " << endl;        printVector(rhsElecField);

        if(iter == 0)
        {
          for(ii=0; ii<DirichletBCs_Epot.size(); ii++)
          {
            nn  = (int) (DirichletBCs_Epot[ii][0]);

            SolnData.var3[nn] += SolnData.var3applied[nn] ;
          }
        }

        rNorm = rhsElecField.norm();

        //cout << " rNorm = " << iter << '\t' << rNorm << endl;

        if(rNorm < tol)
          break;

        if(update_matrix)
          solverElecField.compute(spmtxElecField);

        solnTemp = solverElecField.solve(rhsElecField);

        //printVector(assyForSolnEpot);
        for(ii=0; ii<epotDOF; ii++)
        {
          SolnData.var3[assyForSolnEpot[ii]] += solnTemp[ii];
        }
    }


    return 0;
}




int ElectromechFEM::calcMassMatrixForExplicitDynamics()
{
    cout << " calcMassMatrixForExplicitDynamics()  STARTED" << endl;

    int  ee, ii, jj, nn, n1, n2, n3, n4;
    int idd = SolnData.ElemProp[0]->id;

    MatrixXd  Mlocal, Mlocal2;
    VectorXd  Flocal(100);

    MassVector.resize(nNode*ndof);
    MassVector.setZero();

    ForceVectorExternal.resize(nNode*ndof);
    ForceVectorExternal.setZero();

    for(ee=0;ee<nElem;ee++)  // loop over all the elements
    {
        elems[ee]->calcMassMatrix(Mlocal, true);
        elems[ee]->calcLoadVector(Flocal);

        //printMatrix(Mlocal);

        //printVector(elems[ee]->forAssyVec);
        nn=elems[ee]->nodeNums.size();

        for(ii=0; ii<nn; ii++)
        {
          n1 = ndof*ii;
          n2 = ndof*elems[ee]->nodeNums[ii];

          for(jj=0; jj<ndof; jj++)
          {
            n3 = n1+jj;
            n4 = n2+jj;

            ForceVectorExternal(n4) += Flocal(n3);

            MassVector(n4) += Mlocal(n3, n3);
          }
        }
    }

    cout << " calcMassMatrixForExplicitDynamics() ENDED" << endl;

    solverEigen->rhsVec.resize(dispDOF);
    solverEigen->rhsVec.setZero();

    solverEigen->rhsVec2.resize(presDOF);
    solverEigen->rhsVec2.setZero();
    solverEigen->var2 = solverEigen->rhsVec2;

    return 0;
}


/*
int  ElectromechFEM::ModalAnalysis(int nModes, bool flag, double fact)
{
    //elementDiffStiffTest(0.01, 1, 6, 6, true);

    cout << " ElectromechFEM::ModalAnalysis() ... STARTED " << endl;
    flag = false;
    flag = true;
    cout << " flag = " << flag << endl;

    MatrixXd  Kglobal(totalDOF, totalDOF), Mglobal(totalDOF, totalDOF), eigen_vectors;
    MatrixXd  Kuu, Mlocal;
    VectorXd  eigen_values, eigvec, Flocal, vecTemp;

    int ee, ii, jj, nn, rr, cc;


    Kglobal.setZero();
    Mglobal.setZero();
    for(ee=0;ee<nElem;ee++)  // loop over all the elements
    {
        elems[ee]->calcStiffnessAndResidual(Kuu, Flocal, firstIter);

        elems[ee]->calcMassMatrix(Mlocal, flag);

        //printVector(elems[ee]->forAssyVec);
        nn=elems[ee]->forAssyVec.size();

        for(ii=0; ii<nn; ii++)
        {
          rr = elems[ee]->forAssyVec[ii];
          if( rr != -1 )
          {
            for(jj=0; jj<nn; jj++)
            {
              cc = elems[ee]->forAssyVec[jj];
              if( cc != -1 )
              {
                 Kglobal(rr, cc) += Kuu(ii,jj);
                 Mglobal(rr, cc) += Mlocal(ii,jj);
              }
            }
          }
        }
    }

    cout << "  Solving eigenvalue problem ... " << endl;

    //GeneralizedSelfAdjointEigenSolver<MatrixXd>  es(Kglobal, Mglobal, EigenvaluesOnly);
    GeneralizedSelfAdjointEigenSolver<MatrixXd>  es(Kglobal, Mglobal);
 
    eigen_values  = es.eigenvalues();
    eigen_vectors = es.eigenvectors();

    cout << " The first " << min(20, (int) eigen_values.rows()) << " natural frequencies are ... " << endl;

    for(ii=0;ii<min(20, (int) eigen_values.rows());ii++)
    //for(int ii=0;ii<eigen_values.rows();ii++)
      printf("\t %5d \t %12.10f \n", (ii+1), sqrt(abs(eigen_values(ii)))/2.0/PI);

    cout << "  Writing Mode shapes ... " << endl;
    int  resln[3]={1,1,1};
    vecTemp.resize(totalDOF);
    filecount=0;
    for(jj=0; jj<min(20, (int) eigen_values.rows()); jj++)
    {
      eigvec = eigen_vectors.col(jj);
      for(ii=0; ii<totalDOF; ii++)
        vecTemp[ii] = abs(vecTemp[ii]);

      fact = vecTemp.maxCoeff();
      eigvec /= fact;
      //cout << jj << '\t' << fact << '\t' << eigvec.maxCoeff() << endl;

      SolnData.var1.setZero();
      for(ii=0; ii<totalDOF; ii++)
      {
        SolnData.var1[assyForSoln[ii]] = eigvec[ii];
      }

      postProcess(4, 7, 1, false, 0, 0, resln);
      filecount++;
    }

    cout << " ElectromechFEM::ModalAnalysis() ... ENDED " << endl;

    return 0;
}
*/



//
int  ElectromechFEM::ModalAnalysis(int nModes, bool flag, double fact)
{
    // to compute inf-sup number

    int ee, ii, jj, nn, rr, cc, size1, size2;

    /*
    MatrixXd  Kglobal, Mglobal, eigen_vectors;
    MatrixXd  Kuu(12,12), Kup(12,1), Kpp(1,1);
    VectorXd  eigen_values, eigvec, Flocal, vecTemp;
    MatrixXd  Vmat, Qmat, Bmat;

    Vmat.resize(dispDOF, dispDOF);
    Bmat.resize(dispDOF, nElem);
    Qmat.resize(dispDOF, dispDOF);

    Kglobal.resize(dispDOF, dispDOF);
    Mglobal.resize(nElem, nElem);

    Kglobal.setZero();
    Mglobal.setZero();
    Vmat.setZero();
    Bmat.setZero();
    Qmat.setZero();

    /////////////////////////////////////////
    // compute and assemble matrices
    /////////////////////////////////////////

    for(ee=0;ee<nElem;ee++)  // loop over all the elements
    {
        //cout << "       elem... : " << (e+1) << endl;

        elems[ee]->toComputeInfSupCondition(Kuu, Kup, Kpp);

        //printVector(elems[ee]->forAssyVec);
        nn=elems[ee]->forAssyVec.size();

        // Kuu
        for(ii=0; ii<nn; ii++)
        {
          rr = elems[ee]->forAssyVec[ii];
          if( rr != -1 )
          {
            for(jj=0; jj<nn; jj++)
            {
              cc = elems[ee]->forAssyVec[jj];
              if( cc != -1 )
              {
                 Kglobal(rr, cc) += Kuu(ii,jj);
              }
            }
          }
        }

        // Kup
        for(ii=0; ii<nn; ii++)
        {
          rr = elems[ee]->forAssyVec[ii];
          if( rr != -1 )
          {
            Bmat(rr,ee) += Kup(ii,0);
          }
        }

        // Kpp
        Mglobal(ee,ee) += Kpp(0,0);
    }
    */


    //
    MatrixXf  Kglobal, Mglobal, eigen_vectors;
    MatrixXd  Kuu(12,12), Kup(12,4), Kpp(4,4);
    VectorXf  eigen_values, eigvec, Flocal, vecTemp;
    MatrixXf  Vmat, Qmat, Bmat;

    Vmat.resize(dispDOF, dispDOF);
    Bmat.resize(dispDOF, presDOF);
    Qmat.resize(dispDOF, dispDOF);

    Kglobal.resize(dispDOF, dispDOF);
    Mglobal.resize(presDOF, presDOF);

    Kglobal.setZero();
    Mglobal.setZero();
    Vmat.setZero();
    Bmat.setZero();
    Qmat.setZero();

    /////////////////////////////////////////
    // compute and assemble matrices
    /////////////////////////////////////////

    for(ee=0;ee<nElem;ee++)  // loop over all the elements
    {
        //cout << "       elem... : " << (ee+1) << endl;

        elems[ee]->toComputeInfSupCondition(Kuu, Kup, Kpp);

        //printMatrix(Kuu);
        //printf("\n\n");
        //printMatrix(Kup);
        //printf("\n\n");
        //printMatrix(Kpp);
        //printf("\n\n");

        //printVector(elems[ee]->forAssyVec);
        size1 = elems[ee]->forAssyVec.size();
        size2 = elems[ee]->forAssyVecPres.size();

        // Kuu
        for(ii=0; ii<size1; ii++)
        {
          rr = elems[ee]->forAssyVec[ii];
          if( rr != -1 )
          {
            for(jj=0; jj<size1; jj++)
            {
              cc = elems[ee]->forAssyVec[jj];
              if( cc != -1 )
              {
                Kglobal(rr, cc) += Kuu(ii,jj);
              }
            }

            for(jj=0; jj<size2; jj++)
            {
              cc = elems[ee]->forAssyVecPres[jj];

              if(cc != -1)
              {
                Bmat(rr, cc) += Kup(ii, jj);
              }
            }
          }
        }

        // Kpp
        for(ii=0; ii<size2; ii++)
        {
          rr = elems[ee]->forAssyVecPres[ii];

          if(rr != -1)
          {
            for(jj=0; jj<size2; jj++)
            {
              cc = elems[ee]->forAssyVecPres[jj];

              if(cc != -1)
              {
                Mglobal(rr, cc) += Kpp(ii,jj);
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
    }
    //

    cout << "  Computing the matrix ... " << endl;

    //SparseMatrixXf  Kglobal_sp = Kglobal.sparseView();
    //SparseMatrixXf  Mglobal_sp = Mglobal.sparseView();

    Kglobal = Kglobal.inverse();
    //Kglobal_sp = Kglobal_sp.inverse();
    //printMatrix(Vmat);
    //printf("\n\n");
    //printMatrix(globalK);
    //printf("\n\n");
    //printMatrix(globalM);
    //printf("\n\n");
    Kglobal = (Bmat.transpose()*Kglobal)*Bmat;

    cout << "  Solving eigenvalue problem ... " << endl;

    GeneralizedSelfAdjointEigenSolver<MatrixXf> es(Kglobal, Mglobal, EigenvaluesOnly);
    //GeneralizedSelfAdjointEigenSolver<SparseMatrixXf> es(Kglobal_sp, Mglobal_sp, EigenvaluesOnly);
    eigen_values = es.eigenvalues();
    //eigen_vectors = es.eigenvectors();

    //EigenSolver<MatrixXd>  es(Kglobal, EigenvaluesOnly);
    //cout << "The eigenvalues of A are:" << endl << es.eigenvalues() << endl;
    //eigen_values = es.eigenvalues().col(0);

    cout << "  Eigen analysis successfully completed ... " << endl;

    cout << " The first " << min(20, (int) eigen_values.rows()) << " eigenvalues are ... " << endl;

    for(int ii=0;ii<min(20, (int) eigen_values.rows());ii++)
      printf("\t %5d \t %12.10f \t %12.10f \n", (ii+1), eigen_values(ii), sqrt(abs(eigen_values(ii))));

    return 0;
}
//


void ElectromechFEM::elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt)
{
    //SolnData.var3[0] = 0.0;    SolnData.var3[1] = 0.0;    SolnData.var3[2] = 0.1;
    //SolnData.var3[3] = 0.1;    SolnData.var3[4] = 0.2;    SolnData.var3[5] = 0.2;


  if(elnum > nElem)
    cout << "  Error! Requested Element number is out of range " << endl;
  else
  {
     cout << endl;
     cout << endl;
     cout << "             ///////////////////////////////////// " << endl;
     cout << "             // " << endl;
     cout << "             //     diffStiffTest for Element # : " << elnum << endl;
     cout << "             // " << endl;
     cout << "             ///////////////////////////////////// " << endl;
     cout << endl;
     cout << endl;
     elems[elnum]->diffStiffTest(ddd, dig, dig2, 1);
     cout << endl;
     cout << endl;
  }

  return;
}




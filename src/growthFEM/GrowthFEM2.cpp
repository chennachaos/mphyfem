
#include "SolverPardisoEigen.h"
#include "GrowthFEM.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "ElementBase.h"
#include "SolutionData.h"
#include "Files.h"
#include <chrono>
#include "util.h"
#include "FunctionsProgram.h"
#include "NewMaterial.h"

#include "SolverPetsc.h"
#include "SolverPardisoPetsc.h"


extern MyTime myTime;
extern List<TimeFunction> timeFunction;
extern Files files;

//using namespace Spectra;




void GrowthFEM::setSolver()
{
    printf("\n     GrowthFEM::setSolver()  .... STARTED ...\n");

    if(solverEigen != NULL)
      delete solverEigen;
    solverEigen = NULL;

    if(solverPetsc != NULL)
      delete solverPetsc;
    solverPetsc = NULL;


    //Eigen::initParallel();
    Eigen::setNbThreads(0);

    switch(SOLVER_LIB_TYPE)
    {
        case  SOLVER_LIB_EIGEN: // SolverEigen ..........................

            solverEigen = (SolverEigen*) new SolverEigen;

            solverEigen->setAlgorithmType(1);

            prepareMatrixPattern();

            if(totalDOF > 0)
            {
              if(solverEigen->initialise(0,0,totalDOF) != 0)
                return;
            }
            //solver->setSolverAndParameters();

            solverEigen->printInfo();

        break;

        //case  5:                                           // PARDISO(sym) with Eigen
        case  SOLVER_LIB_EIGEN_PARDISO: // PARDISO(unsym) with Eigen

            solverEigen = (SolverEigen*) new SolverPardisoEigen;

            prepareMatrixPattern();

            //if(slv == 5)
            //{
              //if (solverEigen->initialise(1, PARDISO_STRUCT_SYM, totalDOF) != 0)
                //return;
            //}
            //if(slv == 6)
            //{
              //solverEigen->initialise(1, PARDISO_STRUCT_SYM, totalDOF);
              solverEigen->initialise(1, PARDISO_UNSYM, totalDOF);
            //}

        break;

        case  SOLVER_LIB_PETSC: // SolverPetsc ..........................

            printf("\n     PETSc Solver ...\n");
        
            solverPetsc = (SolverPetsc*) new SolverPetsc;

            prepareMatrixPattern();

            if(totalDOF > 0)
            {
              if(solverPetsc->initialise(0, 0, totalDOF) != 0)
                return;
            }

            solverPetsc->setSolverAndParameters();

        break;

        //case  9: // PARDISO(sym) with Petsc
        case  SOLVER_LIB_PETSC_PARDISO: // PARDISO(unsym) with Petsc

            solverPetsc = (SolverPetsc*) new SolverPardisoPetsc;

            prepareMatrixPattern();

            //if(slv == 9)
            //{
              //if(solverPetsc->initialise(numProc, PARDISO_STRUCT_SYM, totalDOF) != 0)
                //return;
            //}
            //if(slv == 10)
            //{
              //if(solverPetsc->initialise(1, PARDISO_UNSYM, totalDOF) != 0)
                //return;
            //}

            solverPetsc->initialise(1, PARDISO_UNSYM, totalDOF);

            solverPetsc->setSolverAndParameters();

            //solverPetsc->printInfo();

        break;

        default: // invalid slv ...................

             cout << " this solver has not been implemented yet!\n\n";

        break;
    }
  //}
  //else
    //prepareMatrixPattern();


    solverOK = true;

    if( (SolnData.tis > 0) )
      setInitialConditions();


    return;
}




int GrowthFEM::prepareMatrixPattern()
{
    printf("\n     GrowthFEM::prepareMatrixPattern()  .... STARTED ...\n");

    int  r, c, r1, c1, count=0, count1=0, count2=0, iii, e, ind, nsize;
    int *tt, *tt1, *tt2,  val1, val2, n1, n2, a, b, ll, pp, nnz;
    int  ind1, ind2, ee, ii, jj, kk, e1, e2;

    /////////////////////////////////////////////////////////////
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    IEN.resize(nElem);

    for(ee=0;ee<nElem;ee++)
    {
      IEN[ee] = elems[ee]->nodeNums;
    }

    cout << " nNode =" << nNode << endl;

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

    cout << " nNode =" << nNode << endl;

    dispDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndof;jj++)
      {
        if(!NodeType[ii][jj])
          ID[ii][jj] = dispDOF++;
      }
    }

    cout << " dispDOF =" << dispDOF << endl;

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

    cout << " Preparing element data " << endl;

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

    totalDOF = dispDOF + presDOF + mpotDOF + tempDOF;

    cout << " nElem    = " << nElem    << endl;
    cout << " nNode    = " << nNode    << endl;
    cout << " npElem   = " << npElem   << endl;
    cout << " ndof     = " << ndof     << endl;
    cout << " dispDOF  = " << dispDOF  << endl;
    cout << " presDOF  = " << presDOF  << endl;
    cout << " totalDOF = " << totalDOF << endl;

    printf("\n element DOF values initialised \n\n");
    printf("\n Preparing matrix pattern \n\n");


  // matrix pattern needs to be prepared only for the implicit solver
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

    printf("\n Preparing matrix pattern DONE \n\n");

        ndofs_local = totalDOF;

        PetscInt  *d_nnz, *o_nnz;

        ierr  = PetscMalloc1(ndofs_local,  &d_nnz);CHKERRQ(ierr);
        ierr  = PetscMalloc1(ndofs_local,  &o_nnz);CHKERRQ(ierr);

        PetscInt  nnz_max_row = 0, size1;
        nnz = 0;
        for(ii=0;ii<totalDOF;ii++)
        {
          findUnique(forAssyMat[ii]);

          size1 = forAssyMat[ii].size();

          nnz_max_row = max(nnz_max_row, size1);

          d_nnz[ii] = size1;

          nnz += size1;
        }

        cout << " nnz " << nnz << endl;


        //if(totalDOF > 0)
        //{
        ierr = MatCreate(PETSC_COMM_WORLD, &solverPetsc->mtx);CHKERRQ(ierr);

        ierr = MatSetSizes(solverPetsc->mtx, ndofs_local, ndofs_local, totalDOF, totalDOF);CHKERRQ(ierr);

        ierr = MatSetFromOptions(solverPetsc->mtx);CHKERRQ(ierr);

        //ierr = MatMPIAIJSetPreallocation(solverPetsc->mtx, 20, d_nnz, 20, o_nnz);CHKERRQ(ierr);
        ierr = MatSeqAIJSetPreallocation(solverPetsc->mtx, 20, d_nnz);CHKERRQ(ierr);

        //ierr = MatSetOption(solverPetsc->mtx, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
        
        MPI_Barrier(MPI_COMM_WORLD);

        cout << " jjjjjjjjjjjjjjjj " << endl;

        ierr  = PetscMalloc1(nnz_max_row,  &colTemp);CHKERRQ(ierr);
        ierr  = PetscMalloc1(nnz_max_row,  &arrayTemp);CHKERRQ(ierr);

        for(jj=0; jj<nnz_max_row; jj++)
          arrayTemp[jj] = 0.0;

        for(ii=0; ii<totalDOF; ii++)
        {
          ind2 = forAssyMat[ii].size();

          //ierr = MatSetValues(solverPetsc->mtx, 1, &ii, ind2, &(forAssyMat[ii][jj]), arrayTemp, INSERT_VALUES);

          for(jj=0;jj<ind2;jj++)
          {
              ierr = MatSetValue(solverPetsc->mtx, ii, forAssyMat[ii][jj], 0.0, INSERT_VALUES);
              //ierr = MatSetValues(solverPetsc->mtx, 1, &ii, 1, &(forAssyMat[ii][jj]), arrayTemp, INSERT_VALUES);
          }
        }
        //}
        cout << " jjjjjjjjjjjjjjjj " << endl;

        VecCreate(PETSC_COMM_WORLD, &solverPetsc->soln);
        //cout << " iiiiiiiiiiiii " << ndofs_local << '\t' << totalDOF << endl;
        ierr = VecSetSizes(solverPetsc->soln, ndofs_local, totalDOF); CHKERRQ(ierr);
        //ierr = VecCreateMPI(PETSC_COMM_WORLD, ndofs_local, totalDOF, &solverPetsc->soln); CHKERRQ(ierr);
        //cout << " bbbbbbbb " << this_mpi_proc << endl;
        ierr = VecSetFromOptions(solverPetsc->soln);CHKERRQ(ierr);
        ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->rhsVec);CHKERRQ(ierr);
        ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->solnPrev);CHKERRQ(ierr);
        ierr = VecDuplicate(solverPetsc->soln, &solverPetsc->reac);CHKERRQ(ierr);

  }//if(IMPLICIT_SOLVER)
  else
  {
    if( SolnData.tis == 400 )                              // fully-explicit scheme
    {
      if(MIXED_ELEMENT)
      {
        if( SolnData.TRULY_INCOMPRESSIBLE)
        {
          cerr << " Fully-explicit scheme is not valid for truly incompressible material models " << endl;
          exit(-1);
        }
      }

      solverEigen->rhsVec2.resize(presDOF);
      solverEigen->rhsVec2.setZero();
      solverEigen->var2     = solverEigen->rhsVec2;
      solverEigen->var2Prev = solverEigen->rhsVec2;
    }
    else  // semi-implicit scheme
    {
      if(MIXED_ELEMENT)
      {
        setSolverDataForTIC();
      }
      else
      {
        cerr << " Semi implicit scheme is not valid for displacement formulation " << endl;
        exit(-2);
      }
    }
  }

    // remove data objects

    solverPetsc->currentStatus = PATTERN_OK;

    // remove data objects
    nodePosData.clear();
    elemConn.clear();
    IEN.clear();
    LM.clear();
    ID.clear();
    forAssyMat.clear();

    printf("\n     GrowthFEM::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return 1;
}


// set the off-diagonal terms for the solver
void GrowthFEM::setSolverDataForTIC()
{
    int  idd = SolnData.ElemProp[0]->id;
    int  ee, ii, jj, size1, size2, row, col;
    vector<int>  vecIntTemp(10);

    if(IMPLICIT_SOLVER)
    {
      for(ee=0; ee<nElem; ee++)
      {
        size1 = elems[ee]->forAssyVec.size();
        size2 = elems[ee]->forAssyVecPres.size();

        //cout << size1 << '\t' << size2 << endl;
        //printVector(elems[ee]->forAssyVecPres);
        //printVector(elems[ee]->forAssyVec);

        for(ii=0; ii<size2; ii++)
        {
          row = elems[ee]->forAssyVecPres[ii];

          if(row != -1)
          {
            row += dispDOF ;
            for(jj=0; jj<size1; jj++)
            {
              col = elems[ee]->forAssyVec[jj];
              if(col != -1)
              {
                solverEigen->mtx.coeffRef(row, col) = 0.0;
                solverEigen->mtx.coeffRef(col, row) = 0.0;
              }
            }
            for(jj=0; jj<size2; jj++)
            {
              col = elems[ee]->forAssyVecPres[jj];
              if(col != -1)
              {
                col += dispDOF;

                solverEigen->mtx.coeffRef(row, col) = 0.0;
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
      } //for(ee=0;)
    }

    return;
}




int GrowthFEM::solve()
{
    if(totalDOF == 0)
      return 0;

    if ( (SolnData.tis > 0) || (NONLINEAR_SOLVER_TYPE == SOLVER_TYPE_STANDARD) )
      solveWithNewtonRaphson();
    else
      solveWithArclength();

    return 0;
}



/*
int GrowthFEM::solveStepDeflation(int iter_max, double tol1)
{
    cout << "GrowthFEM::solveStepDeflation " << endl;
    SolnData.var1Incr.setZero();

    double  factor = timeFunction[0].prop;
    cout << "factor = " << factor << endl;

    int solns_max = 2;
    double  deflationshift = 1.;
    int  num_solns = 0, solncount, i, j, ii, jj, iter, count, count2;
    vector<int>  soln_indices;
    vector<bool>  soln_conv_flag(solns_max, false);
    vector<VectorXd>  soln_vec(solns_max);
    VectorXd  solnVecTempI, solnVecTempJ;
    count = SolnData.var1.rows()+SolnData.var2.rows();
    VectorXd  dM(count), var1var2Full(count), dM1(totalDOF);
    var1var2Full.setZero();
    int  resln[3];


    double  rNorm = -1.0;
    double  rNormPrev = -1.0;
    double  p = 2.0, normtempj, temp;
    bool  convergenceflag = false;

    for(solncount = 0; solncount<solns_max; solncount++)
    {
        printf("Solution %d \n\n", solncount);

        if(solncount > 0)  SolnData.var1applied.setZero();

        for(iter=0; iter<iter_max; iter++)
        {
            if(iter == 0)
              firstIter = true;
            else
              firstIter = false;

            updateIterStep();


            count = SolnData.var1.rows();
            for(ii=0; ii<count; ii++)
              var1var2Full[ii] = SolnData.var1[ii];

            count2 = SolnData.var2.rows();
            for(ii=0; ii<count2; ii++)
              var1var2Full[count+ii] = SolnData.var2[ii];

            //cout << " aaaaaaaaaaa " << endl;
            calcStiffnessAndResidual();

            //cout << " bbbbbbbbbbb " << endl;
            for(int ii=0; ii<dispDOF; ii++)
              solverEigen->rhsVec[ii] += (factor*solverEigen->Fext[assyForSoln[ii]]);

            rNorm = solverEigen->rhsVec.norm();

            COUT << "GrowthFEM"; printf("    %d \t %11.4e\n", (iter+1), rNorm);

            if( converged() )
            {
              convergenceflag = true;
              SolnData.saveSolution();
              break;
            }

            //cout << " ccccccccccccc " << endl;
            solverEigen->factoriseAndSolve();


            double   M = 1.0;
            dM.setZero();

            for(j=0; j<=solncount-1; j++)
            {
              cout << j << '\t' << soln_conv_flag[j] << endl;

              if(soln_conv_flag[j])
              {
                solnVecTempJ = var1var2Full-soln_vec[j];

                //for(ii=0; ii<var1var2Full.rows(); ii++)
                  //cout << ii << '\t' << var1var2Full[ii] << '\t' << soln_vec[j][ii] << endl;

                normtempj = solnVecTempJ.norm();
                cout << "normtempj = " <<  normtempj  << endl;

                M *= (deflationshift + 1.0/pow(normtempj, p));

                temp = 1.0;
                for(i=0; i<=solncount-1; i++)
                {
                  if(soln_conv_flag[i])
                  {
                    if(j != i)
                    {
                      solnVecTempI = var1var2Full-soln_vec[i];

                      temp *= (deflationshift + 1.0/pow(solnVecTempI.norm(), p));
                    }
                  }
                }

                temp *= -p/pow(normtempj, p+2);

                dM += (solnVecTempJ*temp);
              }
            }

            if(solncount > 0)
            {
                for(ii=0; ii<dispDOF; ii++)
                {
                  dM1[ii] = dM[assyForSoln[ii]];
                }
                count = SolnData.var1.rows();
                for(ii=0; ii<presDOF; ii++)
                {
                  dM1[dispDOF+ii] = dM[count+assyForSolnPres[ii]];
                }

                cout << "M = " <<  M  << endl;

                double  fact = solverEigen->soln.dot(dM1);
                fact = -M/(M-fact);

                solverEigen->soln *= fact;
            }


            for(ii=0; ii<dispDOF; ii++)
            {
              SolnData.var1Incr[assyForSoln[ii]] = solverEigen->soln[ii];
            }
            SolnData.var1 += SolnData.var1Incr;

            if(MIXED_ELEMENT)
            {
              //printVector(assyForSolnPres);
              for(ii=0; ii<presDOF; ii++)
              {
                //cout << ii << '\t' << assyForSolnPres[ii] << '\t' << solverEigen->soln[dispDOF+ii] << endl;
                SolnData.var2[assyForSolnPres[ii]] += solverEigen->soln[dispDOF+ii];
              }
              //printVector(SolnData.var2);
            }
            //cout << " ddddddddddddd " << endl;
            //printVector(SolnData.var1);

            postProcess(4, 10, 0, true, 0.0, 1.0, resln);
        }


        if(convergenceflag)
        {
            soln_conv_flag[solncount] = true;
            num_solns++;
            soln_indices.push_back(solncount);
            soln_vec[solncount] = var1var2Full;

            //var1var2Full = var1var2Full*1.12345678 + VectorXd::Ones(var1var2Full.rows())*0.123456789;
            //var1var2Full = var1var2Full + VectorXd::Ones(var1var2Full.rows())*0.0123456789;

            //SolnData.var1 *= 1.123456789;
            //SolnData.var2 *= 1.123456789;
            //postProcess(4, 10, 0, true, 0.0, 1.0, resln);

            for(ii=0; ii<dispDOF; ii++)
            {
              jj = assyForSoln[ii];
              // + 0.00123456789;
              SolnData.var1[jj] = SolnData.var1[jj]*1.0123456789 + 0.123456789;
            }
            for(ii=0; ii<presDOF; ii++)
            {
              jj = assyForSolnPres[ii];
              // + 0.123456789;
              SolnData.var2[jj] = SolnData.var2[jj]*1.0123456789;
            }
            //setInitialConditions();

            //SolnData.var2 = SolnData.var2*1.123456789 + VectorXd::Ones(SolnData.var2.rows())*0.0123456789;
            //SolnData.var2 = SolnData.var2*1.123456789 + VectorXd::Ones(SolnData.var2.rows())*0.0123456789;

            //count = SolnData.var1.rows();
            //for(ii=0; ii<dispDOF; ii++)
              //SolnData.var1[assyForSoln[ii]] = var1var2Full[assyForSoln[ii]];

            //count2 = SolnData.var2.rows();
            //for(ii=0; ii<presDOF; ii++)
              //SolnData.var2[assyForSolnPres[ii]] = var1var2Full[count+assyForSolnPres[ii]];
        }
    }

    if( !convergenceflag )
      return -1;


    return 0;
}
*/



int GrowthFEM::solveStepDeflation()
{
    cout << "GrowthFEM::solveStepDeflation " << endl;
    SolnData.var1Incr.setZero();

    double  factor = timeFunction[0].prop, xx, yy;
    cout << "factor = " << factor << endl;

    int solns_max = 3;
    double  deflationShift = 0.000123456789, deflationPower = 20.0;
    int  num_solns = 0, solncount, i, j, ii, jj, iter, count, count2;
    vector<int>  soln_indices;
    vector<bool>  soln_conv_flag(solns_max, false);
    vector<VectorXd>  soln_vec(solns_max);
    VectorXd  solnVecTempI, solnVecTempJ;
    //count = SolnData.var1.rows()+SolnData.var2.rows();
    count = SolnData.var1.rows();
    VectorXd  dM(count), var1var2Full(count), dM1(totalDOF);
    var1var2Full.setZero();
    int  resln[3];


    double  rNorm = -1.0;
    double  rNormPrev = -1.0;
    double  normtempj, temp;
    bool  convergenceflag = false;


    for(solncount = 0; solncount<solns_max; solncount++)
    {
        printf("Solution %d \n\n", solncount);

        if(solncount > 0)
        {
            SolnData.var1.setZero();
            //SolnData.var1applied.setZero();
            //iter_max = 11;

            for(ii=0; ii<dispDOF; ii++)
            {
                jj = assyForSoln[ii];
                //SolnData.var1[jj] = SolnData.var1[jj]*1.000123456789;
                SolnData.var1[jj] = 0.00000123456789;
            }
        }

        for(iter=0; iter<max_iters; iter++)
        {
            if(iter == 0)
              firstIter = true;
            else
              firstIter = false;

            updateIterStep();

            //cout << " aaaaaaaaaaa " << endl;
            try
            {
              calcStiffnessAndResidual();
            }
            catch(runtime_error& err)
            {
              //SolnData.var1 = soln_vec[0];
              //updateIterStep();
              cout << err.what() << endl;
              break;
            }

            //cout << " bbbbbbbbbbb " << endl;
            for(int ii=0; ii<dispDOF; ii++)
              solverEigen->rhsVec[ii] += (factor*solverEigen->Fext[assyForSoln[ii]]);

            //printVector(solverEigen->rhsVec);

            rNorm = solverEigen->rhsVec.norm();

            printf(" GrowthFEM   %d \t %11.4e\n", (iter+1), rNorm);

            if( converged() )
            {
              convergenceflag = true;
              SolnData.saveSolution();
              break;
            }

            //cout << " ccccccccccccc " << endl;
            solverEigen->factoriseAndSolve();


            double   M = 1.0;
            dM.setZero();

            for(j=0; j<=solncount-1; j++)
            {
              cout << j << '\t' << soln_conv_flag[j] << endl;

              if(soln_conv_flag[j])
              {
                solnVecTempJ = SolnData.var1-soln_vec[j];

                //for(ii=0; ii<var1var2Full.rows(); ii++)
                  //cout << ii << '\t' << SolnData.var1[ii] << '\t' << soln_vec[j][ii] << endl;

                normtempj = solnVecTempJ.norm();
                cout << "normtempj = " <<  normtempj  << endl;

                M = M * (deflationShift + 1.0/pow(normtempj, deflationPower));

                temp = 1.0;
                for(i=0; i<=solncount-1; i++)
                {
                  if(soln_conv_flag[i])
                  {
                    if(j != i)
                    {
                      solnVecTempI = SolnData.var1-soln_vec[i];

                      temp = temp * (deflationShift + 1.0/pow(solnVecTempI.norm(), deflationPower));
                    }
                  }
                }

                temp = temp * -deflationPower/pow(normtempj, deflationPower+2);

                dM = dM + (solnVecTempJ*temp);
              }
            }

            if(solncount > 0)
            {
                for(ii=0; ii<dispDOF; ii++)
                {
                  dM1[ii] = dM[assyForSoln[ii]];
                }
                count = SolnData.var1.rows();
                for(ii=0; ii<presDOF; ii++)
                {
                  dM1[dispDOF+ii] = dM[count+assyForSolnPres[ii]];
                }

                double  fact = dM1.dot(solverEigen->soln);
                //double  fact = SolnData.var1.dot(dM);
                fact = 1.0/(1.0-fact/M);

                cout << "M = " <<  M  << '\t' << fact <<  endl;

                solverEigen->soln *= fact;
            }


            for(ii=0; ii<dispDOF; ii++)
            {
              SolnData.var1Incr[assyForSoln[ii]] = solverEigen->soln[ii];
            }
            SolnData.var1 += SolnData.var1Incr;

            postProcess();
        }


        if(convergenceflag)
        {
            soln_conv_flag[solncount] = true;
            num_solns++;
 
            soln_indices.push_back(solncount);
 
            soln_vec[solncount] = SolnData.var1;

            postProcess();
        }
    }

    //postProcess();

    if( !convergenceflag )
      return -1;


    return 0;
}



/*
int GrowthFEM::solveWithNewtonRaphson()
{
    cout << "GrowthFEM::solveWithNewtonRaphson " << endl;

    SolnData.var1Incr.setZero();

    double  loadFactor = timeFunction[0].prop;
    cout << "loadFactor = " << loadFactor << endl;

    double  DtFactor = 1.0;                                //dt/dtPrev;
    //double  DtFactor = dt/dtPrev;

    //SolnData.var1  = (1.0+DtFactor)*SolnData.var1Prev - DtFactor*SolnData.var1Prev2;
    //SolnData.var2  = (1.0+DtFactor)*SolnData.var2Prev - DtFactor*SolnData.var2Prev2;

    //printf("time step = %14.10f \t %14.10f \t %14.10f \n", dt, dtPrev, DtFactor);


    int idd = SolnData.ElemProp[0]->id;

    if( (idd == 6002) || (idd == 6005) || (idd == 6052) )
    {
        for(int ee=0; ee<nElem; ee++)
        {
            elems[ee]->presDOF = (1.0+DtFactor)*elems[ee]->presDOFprev - DtFactor*elems[ee]->presDOFprev2;
        }
    }

    rNormPrev = rNorm = -1.0;

    for(int iter=0; iter<max_iters; iter++)
    {
        firstIter = (iter == 0);

        // Compute the velocity and acceleration and the respective values at n+af, n+am
        updateIterStep();

        //cout << " aaaaaaaaaaa " << endl;
        // compute the global stiffness matrix and residual
        calcStiffnessAndResidual();
        //cout << " bbbbbbbbbbb " << endl;


        // compute contributions from the external loads
        addExternalForces();

        for(int ii=0; ii<dispDOF; ii++)
        solverEigen->rhsVec[ii] += (loadFactor*solverEigen->Fext[assyForSoln[ii]]);

        //cout << " RHS " << endl;                           //        printVector(solverEigen->rhsVec);
        //for(int ii=dispDOF; ii<totalDOF; ii++)
            //cout << ii << '\t' << solverEigen->rhsVec[ii] << endl;

        rNormPrev = rNorm;
        rNorm = solverEigen->rhsVec.norm();

        printf(" GrowthFEM ...  %d \t %11.4e \n", (iter+1), rNorm);
        //printf(" GrowthFEM ...  %d \t %11.4e \t %11.4e \n", (iter+1), rNormPrev, rNorm);

        // check for convergence and divergence of the iterations
        if( converged() )
        {
          printf(" GrowthFEM ...  Iterations CONVERGED \n");
          break;
        }
        else if( (rNorm/rNormPrev) > 10000.0 )
        {
          printf(" GrowthFEM ...  Iterations are diverging. NR loop is terminated. \n");
          break;
        }

        //cout << " ccccccccccccc " << endl;
        // solve the matrix system and update the unknown DOFs
        factoriseSolveAndUpdate();
        //cout << " ddddddddddddd " << endl;
        //printVector(SolnData.var1);
    }

    //dtPrev2 = dtPrev;
    //dtPrev  = dt;

    // if the residual is converged, then save the DOFs vectors
    if( converged() )
    {
        SolnData.saveSolution();

        if( (idd == 6002) || (idd == 6005) || (idd == 6052) )
        {
            for(int ee=0; ee<nElem; ee++)
            {
                elems[ee]->presDOFprev2 = elems[ee]->presDOFprev;
                elems[ee]->presDOFprev  = elems[ee]->presDOF;
            }
        }

        //dt = min(2.0*dt, dtMax);

        return 0;
    }
    else
    {
        //dt      = 0.5*dt;

        return -1;
    }
}
*/

int GrowthFEM::solveWithNewtonRaphson()
{
    cout << "\nGrowthFEM::solveWithNewtonRaphson \n\n\n" << endl;

    int  stepsCompleted=0, err = 0;

    postProcess();
    writeNodalData();

    convergedFlagPrev = convergedFlag = false;

    //setInitialConditions();

    //Time loop
    for(int tstep=1; tstep <= max_steps; tstep++)
    {
        myTime.update();


        printf(" ==================================================================== \n");
        printf(" Time step number = %d \n", tstep);
        printf(" Time step size   = %f \n", myTime.dt);
        printf(" Current time     = %f \n", myTime.cur);
        printf(" ==================================================================== \n");


        if(myTime.cur > (tfinal+EPSILON))
        {
            cout << "\n\n Simulation reached the specified final time ... \n\n\n" << endl;
            break;
        }


        // do a time update: reset variables and flags, and prediction step of fields
        timeUpdate();


        double  loadFactor = timeFunction[0].prop;
        cout << "loadFactor = " << loadFactor << endl;

        convergedFlagPrev = convergedFlag;
        convergedFlag = false;

        rNormPrev = rNorm = -1.0;

        cout << endl;    cout << endl;
        for(int iter=1; iter <= max_iters; iter++)
        {
            firstIter = (iter == 1);

            // Compute the velocity and acceleration and the respective values at n+af, n+am
            updateIterStep();

            //printVector(SolnData.var1DotDot);

            // compute the global stiffness matrix and residual

            try
            {
                calcStiffnessAndResidual();
            }
            catch(runtime_error& err)
            {
                cerr << err.what() << endl;
                break;
            }

            // compute contributions from the external loads
            addExternalForces();


            for(int ii=0; ii<dispDOF; ii++)
            {
              //solverEigen->rhsVec[ii] += (loadFactor*solverEigen->Fext[assyForSoln[ii]]);
              VecSetValue(solverPetsc->rhsVec, ii, loadFactor*solverPetsc->Fext[assyForSoln[ii]], ADD_VALUES);
            }

            //cout << " RHS " << endl;                           //        printVector(solverEigen->rhsVec);
            //for(int ii=dispDOF; ii<totalDOF; ii++)
                //cout << ii << '\t' << solverEigen->rhsVec[ii] << endl;


            rNormPrev = rNorm;

            VecAssemblyBegin(solverPetsc->rhsVec);
            VecAssemblyEnd(solverPetsc->rhsVec);

            VecNorm(solverPetsc->rhsVec, NORM_2, &rNorm);


            printf(" MagnetomechFEM ...  %3d \t %11.4e \n", (iter), rNorm);

            // check for convergence and divergence of the iterations
            if( converged() )
            {
              printf("\n MagnetomechFEM ...  Iterations CONVERGED \n\n\n");

              convergedFlag = true;

              break;
            }
            else if( (iter > 3) && diverging(1.0e7) )
            {
              printf(" MagnetomechFEM ...  Iterations are diverging. NR loop is terminated. \n\n\n");
              break;
            }


            // solve the matrix system and update the unknown DOFs
            factoriseSolveAndUpdate();
            //printVector(SolnData.var1);
            postProcess();
        }

        // if the residual is converged, then save the DOFs vectors
        if( convergedFlag )
        {
            stepsCompleted++;

            postProcess();

            writeNodalData();

            saveSolution();

            myTime.stck();
        }
        else
        {
            myTime.cut();

            //reset();
        }
    }

    return err;
}





int GrowthFEM::solveWithArclength()
{
    cout << "GrowthFEM::solveWithArclength " << endl;
    cout << "SolnData.timeStepCount = " << SolnData.timeStepCount << endl;

    SolnData.var1Incr.setZero();

    VectorXd  DuFull(SolnData.var1.rows()), du(SolnData.var1.rows());
    double  Dl, dl, DsFactor;
    int resln[3];

    if(SolnData.timeStepCount == 1)
    {
      loadfactor = myTime.dt;
    }


    if(SolnData.timeStepCount > 1)
    {
      DsFactor = arclenIncr/arclenIncrPrev;

      SolnData.var1  = (1.0+DsFactor)*SolnData.var1Prev - DsFactor*SolnData.var1Prev2;
      //SolnData.var2  = (1.0+DsFactor)*SolnData.var2Prev - DsFactor*SolnData.var2Prev2;
      loadfactor     = (1.0+DsFactor)*loadfactorPrev - DsFactor*loadfactorPrev2;

      printf("arclenIncr = %14.10f \t %14.10f \t %14.10f \n", arclenIncr, arclenIncrPrev, DsFactor);
      printf("loadfactor = %14.10f \t %14.10f \t %14.10f \n", loadfactor, loadfactorPrev, DsFactor);
    }

    cout << "arclenIncr = " << arclenIncr << '\t' << arclenIncrPrev << '\t' << DsFactor << endl;
    cout << "loadfactor = " << loadfactor << '\t' << loadfactor << endl;

    DuFull = SolnData.var1 - SolnData.var1Prev;
    Dl = loadfactor - loadfactorPrev;

    convergedFlagPrev = convergedFlag;
    convergedFlag = false;


    for(int iter=0; iter<max_iters; iter++)
    {
        firstIter = (iter == 0);

        updateIterStep();

        try
        {
            calcStiffnessAndResidual();
        }
        catch(runtime_error& err)
        {
            cout << err.what() << endl;
            break;
        }


        for(int ii=0; ii<dispDOF; ii++)
          solverEigen->rhsVec[ii] += (loadfactor*solverEigen->Fext[assyForSoln[ii]]);

        //printVector(solverEigen->rhsVec);

        rNorm = solverEigen->rhsVec.norm();

        printf("\t GrowthFEM ... %d \t %11.4E \n", (iter+1), rNorm);

        if( rNorm < tol )
        {
          convergedFlag = true;
          break;
        }

        //solverEigen->currentStatus = ASSEMBLY_OK;
        solverEigen->solveArclengthSystem(SolnData.timeStepCount, DuFull, Dl, arclenIncr, dl);


        for(int ii=0; ii<dispDOF; ii++)
        {
          SolnData.var1Incr[assyForSoln[ii]] = solverEigen->soln[ii];
          DuFull[assyForSoln[ii]] += solverEigen->soln[ii];
        }
        SolnData.var1 += SolnData.var1Incr;

        if(MIXED_ELEMENT)
        {
          //printVector(assyForSolnPres);
          for(int ii=0; ii<presDOF; ii++)
          {
            //cout << ii << '\t' << assyForSolnPres[ii] << '\t' << solverEigen->soln[dispDOF+ii] << endl;
            //SolnData.var2[assyForSolnPres[ii]] += solverEigen->soln[dispDOF+ii];
            SolnData.var2[ii] += solverEigen->soln[dispDOF+ii];
          }
          //printVector(SolnData.var2);
        }

        loadfactor += dl;
        Dl += dl;

        //postProcess(4, 10, 0, true, 0.0, 1.0, resln);

        //printVector(solverEigen->soln);
        //cout << endl;  cout << endl;  cout << endl;
        //cout << "dl         " << dl << endl;
        //cout << "loadfactor " << loadfactor << endl;
        //cout << endl;  cout << endl;  cout << endl;
    }

    if(convergedFlag)
    {
      //printVector(SolnData.var1);

      if(SolnData.timeStepCount == 1)
      {
        arclenIncr = sqrt(DuFull.dot(DuFull) + loadfactor*loadfactor*solverEigen->Fext.dot(solverEigen->Fext));

        arclenIncrMax = arclenIncr;
        arclenIncrMin = arclenIncr/1024.0;
      }

      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;

      SolnData.var1Prev2 = SolnData.var1Prev;
      SolnData.var1Prev  = SolnData.var1;

      SolnData.var2Prev2 = SolnData.var2Prev;
      SolnData.var2Prev  = SolnData.var2;

      arclenIncrPrev = arclenIncr;
      if(convergedFlagPrev)
        arclenIncr = min(max(2.0*arclenIncr, arclenIncrMin), arclenIncrMax);

      loadfactorVec.push_back(loadfactor);

      loadStepConverged = loadStepConverged + 1;

      postProcess();

      writeNodalData();
    }
    else
    {
      if(convergedFlagPrev)
        arclenIncr = max(arclenIncr*0.5, arclenIncrMin);
      else
        arclenIncr = max(arclenIncr*0.25, arclenIncrMin);
    }

    cout << " arclenIncr = " << arclenIncr << endl;

    return 0;
}




int GrowthFEM::solveStepArcLengthGrowthModel()
{
    cout << "GrowthFEM::solveStepArcLengthGrowthModel " << endl;
    cout << "SolnData.timeStepCount = " << SolnData.timeStepCount << endl;

    SolnData.var1Incr.setZero();

    VectorXd  DuFull(SolnData.var1.rows()), du(SolnData.var1.rows());
    double  Dl, dl, DsFactor;

    if(SolnData.timeStepCount == 1)
    {
      loadfactor = myTime.dt;
    }


    if(SolnData.timeStepCount > 1)
    {
      DsFactor = arclenIncr/arclenIncrPrev;

      SolnData.var1  = (1.0+DsFactor)*SolnData.var1Prev - DsFactor*SolnData.var1Prev2;
      loadfactor     = (1.0+DsFactor)*loadfactorPrev - DsFactor*loadfactorPrev2;
    }

    cout << "arclenIncr = " << arclenIncr << '\t' << arclenIncrPrev << '\t' << DsFactor << endl;
    cout << "loadfactor = " << loadfactor << '\t' << loadfactor << endl;

    DuFull = SolnData.var1 - SolnData.var1Prev;
    Dl = loadfactor - loadfactorPrev;

    convergedFlagPrev = convergedFlag;
    convergedFlag = false;


    for(int iter=0; iter<max_iters; iter++)
    {
        if(iter == 0)
          firstIter = true;
        else
          firstIter = false;

        updateIterStep();

        calcStiffnessAndResidual();

        for(int ii=0; ii<dispDOF; ii++)
          solverEigen->rhsVec[ii] += (loadfactor*solverEigen->Fext[assyForSoln[ii]]);

        computerInternalForceDerivativeOfGrowth();

        rNorm = solverEigen->rhsVec.norm();

        printf("  GrowthFEM  %d \t %11.4e\n", (iter+1), rNorm);

        if( rNorm < tol )
        {
          convergedFlag = true;
          break;
        }

        //solverEigen->currentStatus = ASSEMBLY_OK;
        solverEigen->solveArclengthSystemGrowthModel(SolnData.timeStepCount, DuFull, Dl, arclenIncr, dl);


        for(int ii=0; ii<dispDOF; ii++)
        {
          SolnData.var1Incr[assyForSoln[ii]] = solverEigen->soln[ii];
          DuFull[assyForSoln[ii]] += solverEigen->soln[ii];
        }
        SolnData.var1 += SolnData.var1Incr;

        if(MIXED_ELEMENT)
        {
          //printVector(assyForSolnPres);
          for(int ii=0; ii<presDOF; ii++)
          {
            //cout << ii << '\t' << assyForSolnPres[ii] << '\t' << solverEigen->soln[dispDOF+ii] << endl;
            SolnData.var2[assyForSolnPres[ii]] += solverEigen->soln[dispDOF+ii];
          }
          //printVector(SolnData.var2);
        }

        loadfactor += dl;
        Dl += dl;

        //printVector(solverEigen->soln);
        //cout << endl;  cout << endl;  cout << endl;
        //cout << "dl " << dl << endl;
        //cout << endl;  cout << endl;  cout << endl;
    }

    //cout << " SolnData.timeStepCount =  "  << SolnData.timeStepCount << endl;

    if(convergedFlag)
    {
      //printVector(SolnData.var1);

      if(SolnData.timeStepCount == 1)
      {
        arclenIncr = sqrt(DuFull.dot(DuFull) + loadfactor*loadfactor*solverEigen->Fext.dot(solverEigen->Fext));

        arclenIncrMax = arclenIncr;
        arclenIncrMin = arclenIncr/1024.0;
      }

      loadfactorPrev2 = loadfactorPrev;
      loadfactorPrev  = loadfactor;

      SolnData.var1Prev2 = SolnData.var1Prev;
      SolnData.var1Prev  = SolnData.var1;

      SolnData.var2Prev2 = SolnData.var2Prev;
      SolnData.var2Prev  = SolnData.var2;

      arclenIncrPrev = arclenIncr;
      if(convergedFlagPrev)
        arclenIncr = min(max(2.0*arclenIncr, arclenIncrMin), arclenIncrMax);

      loadfactorVec.push_back(loadfactor);

      loadStepConverged = loadStepConverged + 1;
    }
    else
    {
      if(convergedFlagPrev)
        arclenIncr = max(arclenIncr*0.5, arclenIncrMin);
      else
        arclenIncr = max(arclenIncr*0.25, arclenIncrMin);
    }

    cout << " arclenIncr = " << arclenIncr << endl;

    return 0;
}




int GrowthFEM::calcStiffnessAndResidual(int printRes, bool zeroMtx, bool zeroRes)
{
    if(DEBUG) {cout << "     GrowthFEM::calcStiffnessAndResidual ...STARTED \n\n";}

    int  ee, ii, jj, nsize;

    MatrixXd  Kuu, Kup, Kpu, Kpp;
    VectorXd  FlocalU, FlocalP;

    vector<int>  vecTemp;

    SolnData.reac.setZero();
    solverPetsc->zeroMtx();

    //printVector(solver->rhsVec);
    //printVector(SolnData.var1applied);
    //printVector(SolnData.var1);

    for(ee=0;ee<nElem;ee++)  // loop over all the elements
    {
        //cout << "       elem... : " << (ee+1) << endl;
        if(MIXED_ELEMENT)
        {
          try
          {
            elems[ee]->calcStiffnessAndResidualMixed(Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP, firstIter);
            //cout << " AAAAAAAAAAA " << endl;
          }
          catch(runtime_error& err)
          {
            cerr << err.what() << endl;
            throw runtime_error("Negative Jacobian encountered");
          }

          elems[ee]->calcLoadVector(FlocalU);
          //cout << " MMMMMMMMMMM " << endl;

          if(firstIter)
            elems[ee]->applyDirichletBCs2field(1, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres, SolnData.var1applied, SolnData.var2applied);
            //elems[ee]->applyDirichletBCsMixed(1, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP);
          //cout << " BBBBBBBBBBB " << endl;

          solverPetsc->assembleMatrixAndVector2field(dispDOF, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres);
          //cout << " MMMMMMMMMMM " << '\t' << ee << endl;
        }
        else
        {
          try
          {
            //cout << " MMMMMMMMMMM " << endl;
            elems[ee]->calcStiffnessAndResidual(Kuu, FlocalU, firstIter);

            //cout << " OOOOOOOOOOO " << endl;
          }
          catch(runtime_error& err)
          {
            cerr << err.what() << endl;
            throw runtime_error("Negative Jacobian encountered");
          }

          elems[ee]->calcLoadVector(FlocalU);
          //cout << " PPPPPPPPPPP " << endl;

          if(firstIter)
            elems[ee]->applyDirichletBCs(Kuu, FlocalU);

          //cout << " MMMMMMMMMMM " << endl;
          //printMatrix(Kuu);
          printVector(FlocalU);
          //printVector(elems[ee]->forAssyVec);
          solverPetsc->assembleMatrixAndVector(0, 0, elems[ee]->forAssyVec, elems[ee]->forAssyVec, Kuu, FlocalU);
        }

        // add up reaction forces
        vecTemp = elems[ee]->globalDOFnums;
        nsize = vecTemp.size();
        for(ii=0;ii<nsize;ii++)
        {
          SolnData.reac[vecTemp[ii]] += FlocalU[ii];
        }
    }

    solverPetsc->currentStatus = ASSEMBLY_OK;

    if(DEBUG) {cout << "     GrowthFEM::calcStiffnessAndResidual ... ENDED \n\n";}

    return 0;
}



int GrowthFEM::factoriseSolveAndUpdate()
{
    if(DEBUG) {cout << "     MagnetomechFEM::factoriseSolveAndUpdate ... STARTED \n\n";}

    time_t tstart, tend;

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec); printf("\n\n\n");
    //tstart = time(0);
    //for(int ii=dispDOF-100; ii<totalDOF; ii++)
      //cout << ii << '\t' << solverEigen->rhsVec[ii] << endl;

    SolnData.var1Incr.setZero();

    // add specified Dirichlet boundary conditions if first iteration
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

        if(mpotDOF != 0)
        {
            for(ii=0; ii<DirichletBCs_Mpot.size(); ii++)
            {
                nn  = (int) (DirichletBCs_Mpot[ii][0]);

                SolnData.var3[nn] += SolnData.var3applied[nn] ;
            }
        }
    }

    if(DEBUG) {cout << "     matrix solution STARTED \n\n";}

    solverPetsc->factoriseAndSolve();

    if(DEBUG) {cout << "     matrix solution DONE \n\n";}

    //tend = time(0);
    //printf("It took %8.4f second(s) \n ", difftime(tend, tstart) );

    //printVector(solverPetsc->soln);

    PetscScalar *arrayTemp;

    VecGetArray(solverPetsc->soln, &arrayTemp);

    // update solution vector
    for(int ii=0; ii<dispDOF; ii++)
    {
        SolnData.var1Incr[assyForSoln[ii]] = arrayTemp[ii];
    }
    SolnData.var1 += SolnData.var1Incr;


    printVector(SolnData.var1);

    
    for(int ii=0; ii<presDOF; ii++)
    {
        //SolnData.var2[assyForSolnPres[ii]] += arrayTemp[dispDOF+ii];
        SolnData.var2[ii] += arrayTemp[dispDOF+ii];
    }

    if(COUPLED_PROBLEM)
    {
        int mpot_offset = dispDOF + presDOF;
        for(int ii=0; ii<mpotDOF; ii++)
        {
            SolnData.var3[assyForSolnMpot[ii]] += arrayTemp[mpot_offset+ii];
        }
    }
    VecRestoreArray(solverPetsc->soln, &arrayTemp);

    // solve pressure variable in the mixed formulation
    // only for the constant pressure elements, Quad4/1, Hex8/1, TRIA6/1, TET10/1

    int  idd = SolnData.ElemProp[0]->id;
    if( (idd == 6002) || (idd == 6052) )
    {
      for(int ee=0;ee<nElem;ee++)  // loop over all the elements
      {
        elems[ee]->solveForPressure();
      }
    }

    if(DEBUG) {cout << "     MagnetomechFEM::factoriseSolveAndUpdate ... ENDED \n\n";}

    return 0;
}





void GrowthFEM::applyBoundaryConditions()
{
    //cout <<  " applying boundary conditions .... " << endl;
    //cout << " tis = " << tis << endl;

    //printVector(SolnData.td);

    int ii, jj, nn, dof, aa, ind;
    double  specVal, PENALTY=1.0e8, val1, val2;

    double  af = SolnData.td(2);

      /*
      for(aa=0;aa<DirichletBCs.size();aa++)
      {
        nn  = (int) (DirichletBCs[aa][0]);
        dof = (int) (DirichletBCs[aa][1]);
        specVal = DirichletBCs[aa][2];

        ind = nn*ndof+dof;
        //specVal = SolnData.var1applied[ind];

        //if( ElemProp[elemConn[0][1]].id == 7 )
          //specVal  -=  SolnData.var1Cur[ind];
        //else
          specVal  -=  SolnData.var1Cur[ind];

        //cout << start << '\t' << nn << '\t' << ind << '\t' << specVal << endl;

        val1 = PENALTY*af;
        val2 = PENALTY*specVal;

        //ind += start;

        solverEigen->mtx.coeffRef(ind, ind) += val1;
        solverEigen->rhsVec[ind]   += val2;
      }
      */

    return;
}



void GrowthFEM::addExternalForces()
{
    if(DEBUG) {cout <<  "    GrowthFEM::addExternalForces() ... STARTED \n\n";}

    int  ee, nn, dof, ii, jj, ind, ind1, ind2;
    double specVal=0.0, loadfactor=0.0, fact1=0.0;

    //MatrixXd  Klocal(30, 30);
    VectorXd  vecTemp(nNode*ndof), vecTempMpot(nNode), Flocal;
    vecTemp.setZero();
    vecTempMpot.setZero();

    ind1 = nNode*(ndof+1);
    if (solverPetsc->Fext.rows() != ind1)
        solverPetsc->Fext.resize(ind1);
    solverPetsc->Fext.setZero();

    // forces due to face/edge loads
    for(ee=0; ee<nElemFaces; ee++)
    {
      //cout << " jjjjjjjjjj " << '\t' << ee << endl;

      elemsFaces[ee]->calcLoadVector(Flocal);
      //elemsFaces[ee]->calcStiffnessAndResidual(Klocal, Flocal);

      //printVector(elemsFaces[ee]->forAssyVec);
      //printMatrix(Klocal);
      //printVector(Flocal);

      //solverPetsc->assembleMatrixAndVectorElecField(elemsFaces[ee]->forAssyVec, Klocal, Flocal);

      //
      for(ii=0; ii<elemsFaces[ee]->nodeNums.size(); ii++)
      {
        ind1 = elemsFaces[ee]->nodeNums[ii]*ndof;
        ind2 = ii*ndof;
        for(jj=0; jj<ndof; jj++)
        {
          solverPetsc->Fext(ind1+jj) += Flocal(ind2+jj);
          //vecTemp(ind1+jj) += Flocal(ind2+jj);
        }
      }
      //
    }
    //cout << "AAAAAAAAAA" << endl;

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

            solverPetsc->Fext[ind] += specVal;
            //vecTemp[ind] += specVal;
        }
        else
        {
            vecTempMpot[nn] += specVal;
        }
    }
    //printVector(vecTemp);

    if(DEBUG) {cout <<  " GrowthFEM::addExternalForces() ... ENDED \n\n" <<  endl;}

    return;
}



void  GrowthFEM::computerInternalForceDerivativeOfGrowth()
{
    int  ee, aa, ii, jj, size1;

    if (solverEigen->FintDerGrowth.rows() != totalDOF)
        solverEigen->FintDerGrowth.resize(totalDOF);
    solverEigen->FintDerGrowth.setZero();

    vector<int>  vecIntTemp;
    VectorXd  Flocal1, Flocal2;

    // forces due to face/edge loads
    for(ee=0; ee<nElem; ee++)
    {
      elems[ee]->calcForceVectorGrowthModel(Flocal1, Flocal2);

      //printVector(elemsFaces[ee]->nodeNums);
      //printVector(Flocal);

      vecIntTemp = elems[ee]->forAssyVec;
      size1 = vecIntTemp.size();

      for(ii=0; ii<size1; ii++)
      {
        aa = vecIntTemp[ii];

        if( aa != -1 )
        {
          solverEigen->FintDerGrowth[aa] += Flocal1[ii];
        }
      }

      if(MIXED_ELEMENT)
      {
        vecIntTemp = elems[ee]->forAssyVecPres;
        size1 = vecIntTemp.size();

        for(ii=0; ii<size1; ii++)
        {
          aa = vecIntTemp[ii];

          if( aa != -1 )
          {
            solverEigen->FintDerGrowth[dispDOF+aa] += Flocal2[ii];
          }
        }
      }
    }

    return;
}



void  GrowthFEM::computeElementErrors(int index)
{
    cout << " GrowthFEM::computeElementErrors " << endl;

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




void GrowthFEM::elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt)
{
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


//
int  GrowthFEM::ModalAnalysis(int nModes, bool flag, double fact)
{
    cout << " GrowthFEM::ModalAnalysis() ... STARTED " << endl;
    cout << " nModes = " << nModes << endl;
    cout << " flag = " << flag << endl;

    //MatrixXd  Kglobal(totalDOF, totalDOF), Mglobal(totalDOF, totalDOF);
    MatrixXcd  eigen_vectors;
    VectorXcd  eigvec, eigen_values, vecTemp;

    int ee, k, ii, jj, nn, rr, cc;

    calcStiffnessAndResidual();

    cout << "  Solving eigenvalue problem ... " << endl;

    /*
    for(k=0; k<solverEigen->mtx.outerSize(); ++k)
    {
      for(SparseMatrixXd::InnerIterator it(solverEigen->mtx,k); it; ++it)
      {
        ii = it.row();
        if(ii >= dispDOF)
        {
          jj = it.col();
          //cout << ii << '\t' << jj << '\t' << it.value() << endl;

          solverEigen->mtx.coeffRef(ii, jj) *= -1.0;
        }
      }
    }
    */
/*
    // Construct matrix operation object using the wrapper class SparseSymMatProd
    SparseGenMatProd<double>  op(solverEigen->mtx);
    //SparseSymMatProd<double>  op(solverEigen->mtx);

    cout << "  Solving eigenvalue problem ... " << endl;

    // Construct eigen solver object, requesting the largest three eigenvalues
    //SymEigsSolver<double, LARGEST_ALGE, DenseSymMatProd<double> > eigs(&op, 3, 6);

    // Construct eigen solver object, requesting the largest three eigenvalues
    //GenEigsSolver<double, SMALLEST_REAL, SparseGenMatProd<double> > eigs(&op, nModes, 2*nModes+5);
    GenEigsSolver<double, SMALLEST_MAGN, SparseGenMatProd<double> > eigs(&op, nModes, 2*nModes+5);
    //SymEigsSolver<double, SMALLEST_MAGN, SparseSymMatProd<double> > eigs(&op, nModes, 2*nModes+5);

    // Initialize and compute
    eigs.init();

    cout << "  Solving eigenvalue problem ... " << endl;

    int nconv = eigs.compute();
    cout << nconv << " eigenvalues converged.\n" << endl;

    // Retrieve results
    //if(eigs.info() == SUCCESSFUL)
    //{
      cout << "  Eigen solver SUCCESSFUL ... " << endl;

      eigen_values = eigs.eigenvalues();
      eigen_vectors = eigs.eigenvectors();

      cout << " The first " << min(nModes, (int) eigen_values.rows()) << " natural frequencies are ... " << endl;

      for(ii=0;ii<min(nModes, (int) eigen_values.rows());ii++)
        cout << ii+1 << '\t' << eigen_values(ii) << '\t' << abs(eigen_values(ii)) << '\t' << sqrt(abs(eigen_values(ii))) << endl;
      //for(int ii=0;ii<eigen_values.rows();ii++)
        //printf("\t %5d \t %12.10f \t %12.10f \n", (ii+1), eigen_values(ii), sqrt(abs(eigen_values(ii))));

      //
      cout << "  Writing Mode shapes ... " << endl;
      int  resln[3]={1,1,1};
      vecTemp.resize(totalDOF);
      filecount=0;
      for(jj=0; jj<min(nModes, (int) eigen_values.rows()); jj++)
      {
        eigvec = eigen_vectors.col(jj);
        for(ii=0; ii<totalDOF; ii++)
        {
          //cout << jj << '\t' << fact << '\t' << eigvec[ii] << endl;
          vecTemp[ii] = abs(vecTemp[ii]);
        }

        //fact = vecTemp.maxCoeff();
        //eigvec /= fact;

        SolnData.var1.setZero();
        for(ii=0; ii<dispDOF; ii++)
        {
          SolnData.var1[assyForSoln[ii]] = eigvec[ii].real();
          //SolnData.var1[assyForSoln[ii]] = eigvec[ii];
        }

        postProcess();
        filecount++;
      }
    //}
    //
*/
    cout << " GrowthFEM::ModalAnalysis() ... ENDED " << endl;

    return 0;
}
//


/*
int  GrowthFEM::ModalAnalysis(int nModes, bool flag, double fact)
{
    cout << " GrowthFEM::ModalAnalysis() ... STARTED " << endl;
    flag = false;
    flag = true;
    cout << " flag = " << flag << endl;
    //totalDOF = 3;

    MatrixXd  Kglobal(totalDOF, totalDOF);
    MatrixXd  eigen_vectors;
    VectorXd  eigvec, eigen_values, vecTemp;

    int ee, ii, jj, nn, rr, cc;

    calcStiffnessAndResidual();

    cout << "  Solving eigenvalue problem ... " << endl;

    Kglobal = solverEigen->mtx;
    SelfAdjointEigenSolver<MatrixXd>  eigs(Kglobal);

    // Retrieve results

    eigen_values = eigs.eigenvalues();
    eigen_vectors = eigs.eigenvectors();

    cout << " The first " << min(nModes, (int) eigen_values.rows()) << " natural frequencies are ... " << endl;

    for(ii=0;ii<min(nModes, (int) eigen_values.rows());ii++)
    //for(int ii=0;ii<eigen_values.rows();ii++)
      printf("\t %5d \t %12.10f \n", (ii+1), sqrt(abs(eigen_values(ii))));

    //
    cout << "  Writing Mode shapes ... " << endl;
    int  resln[3]={1,1,1};
    vecTemp.resize(totalDOF);
    filecount=0;
    for(jj=0; jj<min(nModes, (int) eigen_values.rows()); jj++)
    {
      eigvec = eigen_vectors.col(jj);
      for(ii=0; ii<totalDOF; ii++)
      {
        //cout << jj << '\t' << fact << '\t' << eigvec[ii] << endl;
        vecTemp[ii] = abs(vecTemp[ii]);
      }

      //fact = vecTemp.maxCoeff();
      //eigvec /= fact;

      SolnData.var1.setZero();
      for(ii=0; ii<totalDOF; ii++)
      {
        SolnData.var1[assyForSoln[ii]] = eigvec[ii];
      }

      postProcess(4, 7, 1, false, 0, 0, resln);
      filecount++;
    }

    //

    cout << " GrowthFEM::ModalAnalysis() ... ENDED " << endl;

    return 0;
}
*/


/*
int  GrowthFEM::ModalAnalysis(int nModes, bool flag, double fact)
{
    cout << " GrowthFEM::ModalAnalysis() ... STARTED " << endl;
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

    cout << " GrowthFEM::ModalAnalysis() ... ENDED " << endl;

    return 0;
}
*/


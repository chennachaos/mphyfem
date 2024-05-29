
#include "SolverPardisoEigen.h"
#include "SolverPardisoPetsc.h"
#include "femSolidmechanics.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "ElementBase.h"
#include "SolutionData.h"
#include <chrono>
#include "util.h"

extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime                myTime;
extern bool debug;

using namespace std;



int femSolidmechanics::setSolver(int slv)
{
    if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femSolidmechanics::setSolver ... STARTED \n");

    solverPetsc.reset();

    //Eigen::initParallel();
    Eigen::setNbThreads(0);

    int numProc=1;

    switch(slv)
    {
        case  SOLVER_LIB_EIGEN: // SolverEigen ..........................

            solverEigen = make_unique<SolverEigen>();

            solverEigen->setAlgorithmType(1);

            prepareMatrixPattern();

            if(ntotdofs_global > 0)
            {
              if(solverEigen->initialise(0,0,ntotdofs_global) != 0)
                return -1;
            }

            solverEigen->printInfo();

        break;

        case  SOLVER_LIB_EIGEN_PARDISO: // PARDISO with Eigen

            solverEigen = make_unique<SolverPardisoEigen>();

            numProc = 1;

            prepareMatrixPattern();

            //if(slv == 5)
            //{
              //if (solverEigen->initialise(numProc, PARDISO_STRUCT_SYM, //ntotdofs_global) != 0)
              //  return -2;
            //}
            //if(slv == 6)
            //{
              solverEigen->initialise(numProc, PARDISO_UNSYM, ntotdofs_global);
            //}

        break;

        case  SOLVER_LIB_PETSC: // SolverPetsc ..........................

            printf("\n     PETSc Solver ...\n");

            solverPetsc = make_unique<SolverPetsc>();

            prepareMatrixPattern();
            printf("\n     PETSc Solver ...\n");

            solverPetsc->setSolverAndParameters();
            printf("\n     PETSc Solver ...\n");

        break;

        case  SOLVER_LIB_PETSC_PARDISO: // PARDISO with Petsc

            solverPetsc = make_unique<SolverPardisoPetsc>();

            prepareMatrixPattern();

            solverEigen->initialise(numProc, PARDISO_UNSYM, ntotdofs_global);

            solverPetsc->setSolverAndParameters();

        break;

        default: // invalid slv ...................

             cout << " this solver has not been implemented yet!\n\n";

        break;
    }

    solverOK = true;

    if( convert2TISSOLID(SolnData.timeIntegrationScheme) != TISSOLID::STATIC )
      setInitialConditions();

    solverPetsc->Fext.resize(nNode_global*ndof);

    setSolverDataForFullyImplicit();

    if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femSolidmechanics::setSolver ... ENDED \n");

    return 0;
}




int femSolidmechanics::prepareMatrixPattern()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n     femSolidmechanics::prepareMatrixPattern()  .... STARTED ...\n\n");

    int  r, c, r1, c1, count=0, count1=0, count2=0, iii, e, ind, nsize;
    int  npElem, val1, val2, n1, n2, a, b, ll, pp, nnz;
    int  ind1, ind2, ee, ii, jj, kk, e1, e2, nn, dof, size1, size2;
    int  side, start1, start2, nr1, nr2, count_diag, count_offdiag, tempInt;

    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////


    vector<vector<int> >  NodeDofArray(nNode_global, vector<int>(ndim, -1)), LM;
    vector<vector<bool> >  NodeDofType(nNode_global, vector<bool>(ndim, false));

    setSpecifiedDOFs_Displacement(NodeDofType);


    dispDOF = 0;
    for(ii=0;ii<nNode_global;++ii)
    {
      for(jj=0;jj<ndof;++jj)
      {
        if(!NodeDofType[ii][jj])
        {
          NodeDofArray[ii][jj] = dispDOF++;
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    assyForSoln.resize(dispDOF);
    count = 0;
    for(ii=0;ii<nNode_global;ii++)
    {
      //printVector(NodeDofArray[ii]);
      for(jj=0;jj<ndof;jj++)
      {
        if( NodeDofArray[ii][jj] != -1)
          assyForSoln[count++] = ii*ndof + jj;
      }
    }

    //cout << " assyForSoln " << endl;
    //printVector(assyForSoln);

    // process DOF data for pressure
    //
    prepareDataForMixedElements();


    ntotdofs_global = dispDOF + presDOF;

    if(ARC_LENGTH)  ntotdofs_global += 1;


    PetscPrintf(MPI_COMM_WORLD, " Mesh statistics ...\n ");
    PetscPrintf(MPI_COMM_WORLD, " nElem_global     = %d \n", nElem_global);
    PetscPrintf(MPI_COMM_WORLD, " nNode_global     = %d \n", nNode_global);
    PetscPrintf(MPI_COMM_WORLD, " ndof             = %d \n", ndof);
    PetscPrintf(MPI_COMM_WORLD, " dispDOF          = %d \n", dispDOF);
    PetscPrintf(MPI_COMM_WORLD, " presDOF          = %d \n", presDOF);
    PetscPrintf(MPI_COMM_WORLD, " ntotdofs_global  = %d \n", ntotdofs_global);


    nNode_owned    = nNode_global;
    ntotdofs_local = ntotdofs_global;
    row_start      = 0;
    row_end        = ntotdofs_global-1;

    cout << "ntotdofs_global = " << ntotdofs_global << endl;

    elem_proc_id.resize(nElem_global);
    for(ee=0; ee<nElem_global; ee++)
      elem_proc_id[ee] = 0;

    node_proc_id.resize(nNode_global);
    for(ee=0; ee<nNode_global; ee++)
      node_proc_id[ee] = 0;

    if(n_mpi_procs > 1)
    {
        // compute first and last row indices of the rows owned by the local processor
        row_start  =  1e9;
        row_end    = -1e9;
        ntotdofs_local = 0;
        for(ii=node_start; ii<=node_end; ii++)
        {
            for(jj=0; jj<ndof; jj++)
            {
              if(NodeDofType[ii][jj] == false)
              {
                ind = NodeDofArray[ii][jj];
                row_start  = min(row_start, ind);
                row_end    = max(row_end,   ind);
                ntotdofs_local++;
              }
            }
        }

        cout << "this_mpi_proc   = " << this_mpi_proc   << '\t'
             << "ntotdofs_local  = " << ntotdofs_local  << '\t'
             << "ntotdofs_global = " << ntotdofs_global << '\t'
             << "row_start       = " << row_start       << '\t'
             << "row_end         = " << row_end << endl;

        // check if the sum of local problem sizes is equal to that of global problem size
        ind=0;
        errpetsc = MPI_Allreduce(&ntotdofs_local, &ind, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

        if(ind != ntotdofs_global)
        {
          cerr << "Sum of local problem sizes is not equal to global size" << endl;
          cout << this_mpi_proc << '\t' << "ind = " << ind << '\t' << ntotdofs_global << endl;
        }
    }
    errpetsc = MPI_Barrier(MPI_COMM_WORLD);
    PetscPrintf(MPI_COMM_WORLD, "\n\n Calculating forAssyVec arrays \n\n");

    vector<vector<int> >  forAssyMat;

    //printVector(elem_proc_id);

    vector<int> nodeNums, forAssyVec, globalDOFnums;
    for(ee=0; ee<nElem_global; ee++)
    {
      //if(elem_proc_id[ee] == this_mpi_proc)
      //{
        nodeNums = elems[ee]->nodeNums;
        npElem = nodeNums.size();
        //printVector(nodeNums);

        nsize = ndof*npElem;

        forAssyVec.resize(nsize);
        globalDOFnums.resize(nsize);

        for(ii=0; ii<npElem; ii++)
        {
          n1 = ndof*ii;
          n2 = ndof*nodeNums[ii];

          kk = nodeNums[ii];

          for(dof=0; dof<ndof; dof++)
          {
            globalDOFnums[n1+dof] = n2+dof;

            forAssyVec[n1+dof] = NodeDofArray[kk][dof];
          }

        }
        elems[ee]->forAssyVec = forAssyVec;
        //printVector(forAssyVec);
        //printVector(globalDOFnums);
      //}
    }
    errpetsc = MPI_Barrier(MPI_COMM_WORLD);
    PetscPrintf(MPI_COMM_WORLD, "\n\n Prepared forAssyVec arrays \n\n");

    PetscPrintf(MPI_COMM_WORLD, "\n\n Element DOF values initialised \n\n");
    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    PetscPrintf(MPI_COMM_WORLD, "\n element DOF values initialised \n\n");
    PetscPrintf(MPI_COMM_WORLD, "\n Preparing matrix pattern \n\n");

    // matrix pattern needs to be prepared only for the implicit solver

    vector<int>::iterator  it;

    forAssyMat.resize(ntotdofs_global);
    //for(ii=row_start; ii<=row_end; ii++)
      //forAssyMat[ii].reserve(500);

    for(ee=0; ee<nElem_global; ee++)
    {
      if(elem_proc_id[ee] == this_mpi_proc)
      {
        forAssyVec = elems[ee]->forAssyVec;
        nsize = forAssyVec.size();

        for(ii=0;ii<nsize;ii++)
        {
          r = forAssyVec[ii];

          if(r != -1)
          {
            //if(r >= row_start && r <= row_end)
            //{
              for(jj=0;jj<nsize;jj++)
              {
                if(forAssyVec[jj] != -1)
                {
                  forAssyMat[r].push_back(forAssyVec[jj]);
                }
              }
            //}
          }
        }
      }
    }
    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    if(presDOF > 0)
    {
          cout << "Adding matrix entries for the displacement-pressure coupling " << endl;

          int  *ttP,  *ttU;
          int  sizeP, sizeU, row, col, offset=dispDOF;

          for(ee=0;ee<nElem_global;ee++)
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


    PetscPrintf(MPI_COMM_WORLD, "\n\n Preparing matrix pattern DONE \n\n");


    PetscInt  *diag_nnz, *offdiag_nnz;
    PetscInt  *colTemp;
    PetscScalar  *arrayTemp;

    errpetsc  = PetscMalloc1(ntotdofs_local,  &diag_nnz);CHKERRQ(errpetsc);
    errpetsc  = PetscMalloc1(ntotdofs_local,  &offdiag_nnz);CHKERRQ(errpetsc);

    PetscInt  nnz_max_row = 0;
    kk = 0;
    for(ii=row_start; ii<=row_end; ii++)
    {
      findUnique(forAssyMat[ii]);
      size1 = forAssyMat[ii].size();

      nnz_max_row = max(nnz_max_row, size1);


      count_diag=0, count_offdiag=0;
      for(it=forAssyMat[ii].begin(); it!=forAssyMat[ii].end(); ++it)
      {
        tempInt = *it;

        if(tempInt >= row_start && tempInt <= row_end)
          count_diag++;
        else
          count_offdiag++;
      }

      //cout << " count_diag ..." << ii << '\t' << count_diag << '\t' << count_offdiag << endl;

      diag_nnz[kk]    = count_diag;
      offdiag_nnz[kk] = count_offdiag;
      kk++;
    }

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);
    PetscPrintf(MPI_COMM_WORLD, "\n\n Initialising petsc solver \n\n");

    // Initialize the petsc solver
    solverPetsc->initialise(ntotdofs_local, ntotdofs_global, diag_nnz, offdiag_nnz);
    errpetsc = MPI_Barrier(MPI_COMM_WORLD);


    PetscPrintf(MPI_COMM_WORLD, " Initialise the Matrix pattern \n", errpetsc);

    errpetsc  = PetscMalloc1(nnz_max_row,  &colTemp);CHKERRQ(errpetsc);
    errpetsc  = PetscMalloc1(nnz_max_row,  &arrayTemp);CHKERRQ(errpetsc);

    for(jj=0; jj<nnz_max_row; jj++)
      arrayTemp[jj] = 0.0;

    PetscScalar  Klocal[nnz_max_row];
    for(ii=0; ii<nnz_max_row; ii++)  Klocal[ii] = 0.0;

    PetscInt  rows[1];
    size1 = 1;
    for(ii=row_start; ii<=row_end; ii++)
    {
      rows[0] = ii;

      forAssyVec = forAssyMat[ii];
      size2 = forAssyVec.size();

      errpetsc = MatSetValues(solverPetsc->mtx, size1, rows, size2, &forAssyVec[0], Klocal, INSERT_VALUES);

      //errpetsc = MatSetValue(solverPetsc->mtx, ii, forAssyMat[ii][jj], 0.0, INSERT_VALUES);
      //errpetsc = MatSetValues(solverPetsc->mtx, 1, &ii, 1, &(forAssyMat[ii][jj]), arrayTemp, INSERT_VALUES);
    }
    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    // Create reaction vector
    errpetsc = VecCreate(PETSC_COMM_WORLD, &(solverPetsc->reacVec));
    CHKERRQ(errpetsc);

    ind1 = nNode_owned*ndof;
    ind2 = nNode_global*ndof;

    errpetsc = VecSetSizes(solverPetsc->reacVec, ind1, ind2);
    CHKERRQ(errpetsc);

    errpetsc = VecSetFromOptions(solverPetsc->reacVec);
    CHKERRQ(errpetsc);

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);


    solverPetsc->currentStatus = PATTERN_OK;


    errpetsc  = PetscFree(diag_nnz);   CHKERRQ(errpetsc);
    errpetsc  = PetscFree(offdiag_nnz);   CHKERRQ(errpetsc);
    errpetsc  = PetscFree(colTemp);   CHKERRQ(errpetsc);
    errpetsc  = PetscFree(arrayTemp);   CHKERRQ(errpetsc);

    PetscPrintf(MPI_COMM_WORLD, "\n     femSolids::prepareMatrixPattern()  .... FINISHED ...\n\n");

    printf("\n     femSolidmechanics::prepareMatrixPattern()  .... FINISHED ...\n\n");

    return 0;
}





// set the off-diagonal terms for the solver
int femSolidmechanics::setSolverDataForFullyImplicit()
{
    cout <<  " femSolidmechanics::setSolverDataForFullyImplicit() ... STARTED " << endl;
/*
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
        size1 = elems[ee]->forAssyVec.size();

        //printVector(elems[ee]->forAssyVec);

        for(ii=0; ii<size1; ii++)
        {
          row = elems[ee]->forAssyVec[ii];

          if(row != -1)
          {
            for(jj=0; jj<size1; jj++)
            {
              col = elems[ee]->forAssyVec[jj];

              //cout << ii << '\t' << jj << '\t' << row << '\t' << col << endl;

              if(col != -1)
              {
                matK.coeffRef(row, col) = 0.0;
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
    } //for(ee=0;)

    matK.makeCompressed();

    //initialise_pardiso();
    solver.analyzePattern(matK);
    solver.factorize(matK);
    //solver.compute(matK);
*/
    cout <<  " femSolidmechanics::setSolverDataForFullyImplicit() ... ENDED " << endl;

    return 0;
}





// set the off-diagonal terms for the solver
int femSolidmechanics::setSolverDataForTIC()
{
    int  idd = 0;//SolnData.ElemProp[0]->id;
    int  ee, ii, jj, size1, size2, row, col;
    vector<int>  vecIntTemp(10);

    if(IMPLICIT_SOLVER)
    {
      for(ee=0; ee<nElem_global; ee++)
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
    else
    {
      cout << " dispDOF = " << dispDOF << endl;
      cout << " presDOF = " << presDOF << endl;

      // inverse of A matrix
      solverEigen->matAinv.setZero();
      solverEigen->matAinv.resize(dispDOF, dispDOF);
      solverEigen->matAinv.reserve(dispDOF);
      for(ii=0; ii<dispDOF; ii++)
      {
        solverEigen->matAinv.coeffRef(ii,ii) = 0.0;
      }

      VectorXi  nnzVec(dispDOF), nnzVec2(presDOF);

      // B matrix
      solverEigen->matB.setZero();
      solverEigen->matB.resize(dispDOF, presDOF);
      solverEigen->matB.reserve(ceil(presDOF*dispDOF*0.1));
      jj = ceil(presDOF*0.2);
      for(ii=0; ii<dispDOF; ii++)
        nnzVec(ii) = jj;
      solverEigen->matB.reserve(nnzVec);

      // C matrix
      solverEigen->matC.setZero();
      solverEigen->matC.resize(presDOF, dispDOF);
      solverEigen->matC.reserve(ceil(presDOF*dispDOF*0.1));
      jj = ceil(dispDOF*0.2);
      for(ii=0; ii<presDOF; ii++)
        nnzVec2(ii) = jj;
      solverEigen->matC.reserve(nnzVec2);

      // D matrix
      solverEigen->matD.setZero();
      solverEigen->matD.resize(presDOF, presDOF);
      solverEigen->matD.reserve(ceil(presDOF*presDOF*0.1));
      jj = ceil(presDOF*0.2);
      for(ii=0; ii<presDOF; ii++)
        nnzVec2(ii) = jj;
      solverEigen->matD.reserve(nnzVec2);

      for(ee=0; ee<nElem_global; ee++)
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
            for(jj=0; jj<size1; jj++)
            {
              col = elems[ee]->forAssyVec[jj];

              if(col != -1)
              {
                solverEigen->matC.coeffRef(row, col) = 0.0;
                solverEigen->matB.coeffRef(col, row) = 0.0;
              }
            }

            for(jj=0; jj<size2; jj++)
            {
              col = elems[ee]->forAssyVecPres[jj];

              if(col != -1)
              {
                solverEigen->matD.coeffRef(row, col) = 0.0;
              }
            }
          }//if(row != -1)
        } //for(ii=0;)
      } //for(ee=0;)

      //cout << " size1  size2 " << endl;

      solverEigen->matAinv.makeCompressed();
      solverEigen->matB.makeCompressed();
      solverEigen->matC.makeCompressed();
      solverEigen->matD.makeCompressed();

      solverEigen->rhsVec.resize(dispDOF);
      solverEigen->rhsVec.setZero();
      //solverEigen->disp     = solverEigen->rhsVec;
      //solverEigen->dispPrev = solverEigen->rhsVec;

      solverEigen->rhsVec2.resize(presDOF);
      solverEigen->rhsVec2.setZero();
      //solverEigen->pres     = solverEigen->rhsVec2;
      //solverEigen->presPrev = solverEigen->rhsVec2;

      solverEigen->currentStatus = PATTERN_OK;
    }

    return 0;
}




int femSolidmechanics::solveFullyImplicit()
{
    if(ntotdofs_global == 0)
      return 0;

    if( (SolnData.timeIntegrationScheme != "STEADY") || (NONLINEAR_SOLVER_TYPE == SOLVER_TYPE_NEWTONRAPHSON) )
      solveWithNewtonRaphson();
    else if(NONLINEAR_SOLVER_TYPE == SOLVER_TYPE_ARCLENGTH)
      solveWithArclength();
    else
    {
        cout << " femSolidmechanics::solveFullyImplicit() ... Algorithm type is not available " << endl;
    }

    return 0;
}





int femSolidmechanics::solveWithNewtonRaphson()
{
    cout << "\femSolidmechanics::solveWithNewtonRaphson \n\n\n" << endl;

    int  stepsCompleted=1, err = 0;

    solverPetsc->Fext.resize(nNode_global*ndof);

    setInitialConditions();
    postProcess();
    writeNodalData();

    convergedFlagPrev = convergedFlag = false;

    //Time loop
    while( (myTime.cur <= (timeFinal-EPSILON)) && (stepsCompleted <= stepsMax) )
    {
        // do a time update: reset variables and flags, and prediction step of fields
        timeUpdate();

        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n");
        PetscPrintf(MPI_COMM_WORLD, " Time step number     =  %d  \n", stepsCompleted);
        PetscPrintf(MPI_COMM_WORLD, " Time step size       =  %f  \n", myTime.dt);
        PetscPrintf(MPI_COMM_WORLD, " Current time         =  %f  \n", myTime.cur);
        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n");

        convergedFlagPrev = convergedFlag;
        convergedFlag = false;

        rhsNormPrev = rhsNorm = -1.0;

        PetscPrintf(MPI_COMM_WORLD, "\n\n");
        for(int iter=1; iter <= iterationsMax; iter++)
        {
            firstIteration = (iter == 1);

            // Compute the velocity and acceleration and the respective values at n+af, n+am
            updateIterStep();
            errpetsc = MPI_Barrier(MPI_COMM_WORLD);

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
            errpetsc = MPI_Barrier(MPI_COMM_WORLD);

            for(int ii=0; ii<dispDOF; ii++)
            {
              //solverEigen->rhsVec[ii] += solverEigen->Fext[assyForSoln[ii]];
              VecSetValue(solverPetsc->rhsVec, ii, solverPetsc->Fext[assyForSoln[ii]], ADD_VALUES);
            }
            errpetsc = MPI_Barrier(MPI_COMM_WORLD);
            //cout << " RHS " << endl;                           //        printVector(solverEigen->rhsVec);
            //for(int ii=dispDOF; ii<totalDOF; ii++)
                //cout << ii << '\t' << solverEigen->rhsVec[ii] << endl;


            rhsNormPrev = rhsNorm;

            VecAssemblyBegin(solverPetsc->rhsVec);
            VecAssemblyEnd(solverPetsc->rhsVec);

            VecNorm(solverPetsc->rhsVec, NORM_2, &rhsNorm);

            PetscPrintf(MPI_COMM_WORLD, " femSolidmechanics ...  %3d \t %11.4e \n", (iter), rhsNorm);

            // check for convergence and divergence of the iterations
            if( converged() )
            {
              PetscPrintf(MPI_COMM_WORLD, "\n femSolidmechanics ...  Iterations CONVERGED \n\n\n");

              convergedFlag = true;

              break;
            }
            else if( (iter > 3) && diverging(1.0e7) )
            {
              PetscPrintf(MPI_COMM_WORLD, " femSolidmechanics ...  Iterations are diverging. NR loop is terminated. \n\n\n");
              break;
            }


            // solve the matrix system and update the unknown DOFs
            factoriseSolveAndUpdate();
            //printVector(SolnData.var1);
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

            reset();
        }
    }

    PetscPrintf(MPI_COMM_WORLD, "\n\n\n Simulation reached the specified final time or maximum steps specified ... \n\n\n");

    return err;
}



int femSolidmechanics::solveWithArclength()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n femSolids::solveWithArclength \n\n");

    int stepsCompleted=1, err = 0;
    VectorXd  DuFull(SolnData.disp.rows()), du(SolnData.disp.rows());
    double  Dl, dl, DsFactor, value;
    int resln[3];

    postProcess();
    writeNodalData();

    // compute contributions from the external loads
    addExternalForces();

    solverPetsc->dispDOF = dispDOF;
    solverPetsc->presDOF = presDOF;
    solverPetsc->assyForSoln = assyForSoln;

    convergedFlagPrev = convergedFlag = false;

    loadFactor = myTime.dt;

    //Time loop
    while( (myTime.cur <= (timeFinal-EPSILON)) && (stepsCompleted <= stepsMax) )
    {
        // do a time update: reset variables and flags, and prediction step of fields
        timeUpdate();

        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n");
        PetscPrintf(MPI_COMM_WORLD, " Time step number     =  %d  \n", stepsCompleted);
        PetscPrintf(MPI_COMM_WORLD, " Load factor          =  %f  \n", loadFactor);
        PetscPrintf(MPI_COMM_WORLD, " ==================================================================== \n\n\n");


        convergedFlagPrev = convergedFlag;
        convergedFlag = false;

        rhsNormPrev = rhsNorm = -1.0;


        SolnData.dispIncr.setZero();


        if(stepsCompleted > 1)
        {
          DsFactor = arclenIncr/arclenIncrPrev;

          SolnData.disp  = (1.0+DsFactor)*SolnData.dispPrev - DsFactor*SolnData.dispPrev2;
          //SolnData.pres  = (1.0+DsFactor)*SolnData.presPrev - DsFactor*SolnData.presPrev2;
          loadFactor     = (1.0+DsFactor)*loadFactorPrev - DsFactor*loadFactorPrev2;
        }

        PetscPrintf(MPI_COMM_WORLD, "arclenIncr = %14.10f \t %14.10f \t %14.10f \n", arclenIncr, arclenIncrPrev, DsFactor);
        PetscPrintf(MPI_COMM_WORLD, "loadFactor = %14.10f \t %14.10f \t %14.10f \n", loadFactor, loadFactorPrev, DsFactor);

        DuFull = SolnData.disp - SolnData.dispPrev;
        Dl = loadFactor - loadFactorPrev;

        convergedFlagPrev = convergedFlag;
        convergedFlag = false;


        for(int iter=1; iter<iterationsMax; iter++)
        {
            firstIteration = (iter == 1);

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
            {
              //solverEigen->rhsVec[ii] += (loadFactor*solverEigen->Fext[assyForSoln[ii]]);
              value = loadFactor*solverPetsc->Fext[assyForSoln[ii]];
              VecSetValue(solverPetsc->rhsVec, ii, value, ADD_VALUES);
            }

            //printVector(solverEigen->rhsVec);

            rhsNormPrev = rhsNorm;

            VecAssemblyBegin(solverPetsc->rhsVec);
            VecAssemblyEnd(solverPetsc->rhsVec);

            VecNorm(solverPetsc->rhsVec, NORM_2, &rhsNorm);

            PetscPrintf(MPI_COMM_WORLD, "\t femSolids ... %d \t %11.4E \n", iter, rhsNorm);

            if( rhsNorm < conv_tol )
            {
              convergedFlag = true;
              break;
            }

            //solverEigen->currentStatus = ASSEMBLY_OK;
            solverPetsc->solveArclengthSystem(stepsCompleted, DuFull, Dl, arclenIncr, dl);

            PetscScalar *arrayTemp;

            VecGetArray(solverPetsc->solnVec, &arrayTemp);

            for(int ii=0; ii<dispDOF; ii++)
            {
              SolnData.dispIncr[assyForSoln[ii]] = arrayTemp[ii];
              DuFull[assyForSoln[ii]] += arrayTemp[ii];
            }

            VecRestoreArray(solverPetsc->solnVec, &arrayTemp);

            SolnData.disp += SolnData.dispIncr;

            if(MIXED_ELEMENT)
            {
              //printVector(assyForSolnPres);
              for(int ii=0; ii<presDOF; ii++)
              {
                //cout << ii << '\t' << assyForSolnPres[ii] << '\t' << solverEigen->soln[dispDOF+ii] << endl;
                //SolnData.var2[assyForSolnPres[ii]] += solverEigen->soln[dispDOF+ii];
                //SolnData.pres[ii] += solverEigen->soln[dispDOF+ii];
              }
              //printVector(SolnData.var2);
            }

            loadFactor += dl;
            Dl += dl;
        }

        if(convergedFlag)
        {
          //printVector(SolnData.var1);

          if(stepsCompleted == 1)
          {
            arclenIncr = sqrt(DuFull.dot(DuFull) + loadFactor*loadFactor*solverPetsc->Fext.dot(solverPetsc->Fext));

            arclenIncrMax = arclenIncr;
            arclenIncrMin = arclenIncr/1024.0;
          }

          loadFactorPrev2 = loadFactorPrev;
          loadFactorPrev  = loadFactor;

          SolnData.dispPrev2 = SolnData.dispPrev;
          SolnData.dispPrev  = SolnData.disp;

          SolnData.presPrev2 = SolnData.presPrev;
          SolnData.presPrev  = SolnData.pres;

          arclenIncrPrev = arclenIncr;
          if(convergedFlagPrev)
            arclenIncr = min(max(2.0*arclenIncr, arclenIncrMin), arclenIncrMax);

          loadFactorVec.push_back(loadFactor);

          loadStepConverged = loadStepConverged + 1;

          postProcess();

          writeNodalData();

          stepsCompleted++;
        }
        else
        {
          if(convergedFlagPrev)
            arclenIncr = max(arclenIncr*0.5, arclenIncrMin);
          else
            arclenIncr = max(arclenIncr*0.25, arclenIncrMin);
        }

        PetscPrintf(MPI_COMM_WORLD, " arclenIncr = %f \n", arclenIncr);
    }

    return 0;
}







int femSolidmechanics::calcStiffnessAndResidual()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femSolidmechanics::calcStiffnessAndResidual ...STARTED \n\n");}

    int  ee, ii, jj, nsize;

    MatrixXd  Kuu, Kup, Kpu, Kpp;
    VectorXd  FlocalU, FlocalP;

    vector<int>  vecTemp;

    SolnData.reac.setZero();
    solverPetsc->zeroMtx();

    for(ee=0;ee<nElem_global;ee++)  // loop over all the elements
    {
    if(elem_proc_id[ee] == this_mpi_proc)
    {
        if(MIXED_ELEMENT)
        {
          try
          {
            elems[ee]->calcStiffnessAndResidualMixed(Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP);
          }
          catch(runtime_error& err)
          {
            cerr << err.what() << endl;
            throw runtime_error("Negative Jacobian encountered");
          }

          if(firstIteration)
            elems[ee]->applyDirichletBCs2field(1, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres, SolnData.dispApplied, SolnData.presApplied);
            //elems[ee]->applyDirichletBCsMixed(1, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP);
          //cout << " BBBBBBBBBBB " << endl;

          solverPetsc->assembleMatrixAndVector2field(dispDOF, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres);
        }
        else
        {
          try
          {
            elems[ee]->calcStiffnessAndResidual(Kuu, FlocalU);
          }
          catch(runtime_error& err)
          {
            cerr << err.what() << endl;
            throw runtime_error("Negative Jacobian encountered");
          }

          if(firstIteration)
            elems[ee]->applyDirichletBCs(Kuu, FlocalU);

          //printMatrix(Kuu);
          //printVector(FlocalU);
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
    }

    solverPetsc->currentStatus = ASSEMBLY_OK;

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femSolidmechanics::calcStiffnessAndResidual ... ENDED \n\n");}

    return 0;
}



int femSolidmechanics::factoriseSolveAndUpdate()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femSolidmechanics::factoriseSolveAndUpdate ... STARTED \n\n");}

    time_t tstart, tend;

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec); printf("\n\n\n");
    //for(int ii=dispDOF-100; ii<totalDOF; ii++)
      //cout << ii << '\t' << solverEigen->rhsVec[ii] << endl;

    SolnData.dispIncr.setZero();

    // add specified Dirichlet boundary conditions if first iteration
    if(firstIteration)
    {
        int ii, dof;
        for(ii=0; ii<dofs_specified_disp.size(); ii++)
        {
            dof = dofs_specified_disp[ii];
            SolnData.dispIncr[dof] += SolnData.dispApplied[dof] ;
        }

        for(ii=0; ii<dofs_specified_pres.size(); ii++)
        {
            dof = dofs_specified_pres[ii];
            SolnData.pres[dof] += SolnData.presApplied[dof] ;
        }
    }

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     matrix solution STARTED \n\n");}

    tstart = time(0);
    if( solverPetsc->factoriseAndSolve() )
    {
        PetscPrintf(MPI_COMM_WORLD, " PETSc solver not converged. \n\n");
        return -1;
    }

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     matrix solution DONE \n\n");}

    tend = time(0);
    PetscPrintf(MPI_COMM_WORLD, "It took %8.4f second(s) \n ", difftime(tend, tstart) );

    //printVector(solverPetsc->soln);

    PetscScalar *arrayTempSoln;
    Vec            vecseq;
    VecScatter     ctx;

    VecScatterCreateToAll(solverPetsc->solnVec, &ctx, &vecseq);


    VecScatterBegin(ctx, solverPetsc->solnVec, vecseq, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctx,   solverPetsc->solnVec, vecseq, INSERT_VALUES, SCATTER_FORWARD);

    VecGetArray(vecseq, &arrayTempSoln);

    // update solution vector
    for(int ii=0; ii<dispDOF; ii++)
    {
        SolnData.dispIncr[assyForSoln[ii]] = arrayTempSoln[ii];
    }
    SolnData.disp += SolnData.dispIncr;


    //if(debug) printVector(SolnData.disp);


    for(int ii=0; ii<presDOF; ii++)
    {
        SolnData.pres[assyForSolnPres[ii]] += arrayTempSoln[dispDOF+ii];
        //SolnData.pres[ii] += arrayTempSoln[dispDOF+ii];
    }

    VecRestoreArray(vecseq, &arrayTempSoln);

    VecScatterDestroy(&ctx);
    VecDestroy(&vecseq);

    // solve pressure variable in the mixed formulation
    // only for the constant pressure elements, Quad4/1, Hex8/1, TRIA6/1, TET10/1

    if( MIXED_ELEMENT_P0 )
    {
      for(int ee=0;ee<nElem_global;ee++)  // loop over all the elements
      {
        elems[ee]->solveForPressure();
      }
    }

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femSolidmechanics::factoriseSolveAndUpdate ... ENDED \n\n");}

    return 0;
}







int  femSolidmechanics::computeElementErrors(int ind)
{
    int  ii, ee, count=0, dd, domTemp;

  for(int index=0; index<3; index++)
  {
    //cout << " index = " << index << endl;

    if(index < 4) // L2 or H1 norm based errors
    {
      if(ndim == 2)
      {
        totalError = 0.0;
        for(int ee=0;ee<nElem_global;ee++)  // loop over all the elements
        {
          elems[ee]->calcError2D(index);

          totalError += elems[ee]->getError();
        }
      }
      else
      {
        totalError = 0.0;
        for(int ee=0;ee<nElem_global;ee++)  // loop over all the elements
        {
          //cout << " ee = " << ee << endl;
          //cout << elems[ee]->getElementTypeNumber() << endl;
          //printVector(elems[ee]->nodeNums);
          //printVector(elems[ee]->forAssyVec);

          elems[ee]->calcError3D(index);
          //cout << " ee = " << ee << endl;

          totalError += elems[ee]->getError();
        }
      }


      totalError = sqrt(totalError);

        if(index == 0)
          printf(" \n\n \t Displacement Error = %12.6E \n\n " , totalError);
        else if(index == 1)
          printf(" \n\n \t Pressure Error     = %12.6E \n\n " , totalError);
        else if(index == 2)
          printf(" \n\n \t Stress Error       = %12.6E \n\n " , totalError);
        else
          printf(" \n\n \t H1 Error = %12.6E \n\n " , totalError);
    }
    else if(index == 10) // total energy
    {
      VectorXd  energyElem(3), energyGlobal(3);
      energyGlobal.setZero();

      for(int ee=0;ee<nElem_global;ee++)  // loop over all the elements
      {
        elems[ee]->computeEnergy(0, 0, energyElem);

        energyGlobal += energyElem;
      }
      energyGlobal[2] = energyGlobal[0] + energyGlobal[1];

      char        tmp[200];
      string    tmpStr;

      sprintf(tmp," \t %12.6E \t %12.6E \t %12.6E", energyGlobal[0], energyGlobal[1], energyGlobal[2]);
    }
    else if(index == 11) // total momentum
    {
      VectorXd  momElem(6), momGlobal(6);
      momGlobal.setZero();

      for(int ee=0;ee<nElem_global;ee++)  // loop over all the elements
      {
        //elems[ee]->computeMomentum(0, 0, momElem);

        momGlobal += momElem;
      }

      char        tmp[200];
      string    tmpStr;

      //sprintf(tmp," \t %12.6E \t %12.6E \t %12.6E", momGlobal[0], momGlobal[1], momGlobal[5]);
    }
    else
    {
      cout << " Compute element errors not defined for thie index value " << endl;
    }
    //printf(" \n\n \t totalError = %12.6E \n\n " , totalError);
  }

    return 0;
}





/*
int femSolidmechanics::solveExplicitStep(int stepsMax)
{
    // using the stiffness matrix
    // Subroutine for CH-alpha scheme

    // just in case this function is called by mistake while using the implicit solver
    if(IMPLICIT_SOLVER)
      return 0;

    //cout << stepsMax << endl;

    MatrixXd  Kuu;
    VectorXd  Flocal, dispTemp, rhsVec2, rhsVec(nNode_global*ndof);

    cout << " Calculating stiffness matrix " << endl;

    solverEigen->zeroMtx();
    for(int ee=0;ee<nElem_global;ee++)  // loop over all the elements
    {
        //cout << "       elem... : " << (ee+1) << endl;
        elems[ee]->calcStiffnessAndResidual(Kuu, Flocal, firstIter);

        solverEigen->assembleMatrixAndVector(0, 0, elems[ee]->forAssyVec, Kuu, Flocal);
    }
    cout << " Calculating stiffness matrix " << endl;


    int  iter, bb, ee, ii, jj, nn, n1, n2, n3, rr, cc;
    int  dof, ind, ind1, ind2;
    int  resln[3]={1,1,1};
    double specVal=0.0, fact=0.0, fact1=0.0, vec[10];

    //int     stepsMax = 10000;
    int     stepsCompleted=1;
    double  dt   = myTime.dt;
    double  timeNow=dt;
    double  timeFinal = stepsMax*dt;
    double  DTT  = dt*dt;
    double  IDTT = 1.0/DTT;
    double  af = 1.0;
    double  am = 1.0;
    double  gamma = 0.5+am;
    //double  beta  = 1.0; // 2nd-order accurate
    //double  beta  = 28.0/27.0; // 2nd-order accurate
    double  beta  = am + 1.0/12.0; // 3rd-order accurate

    char fname[200];
    sprintf(fname,"%s%s", files.Tfile.asCharArray(),"-explicit.dat");

    ofstream fout(fname);

    if(fout.fail())
    {
       cout << " Could not open the Output file" << endl;
       exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(8);

    //printVector(SolnData.dispDot);

    // Store variables
    SolnData.dispPrev2       =  SolnData.dispPrev;
    SolnData.dispPrev        =  SolnData.disp;
    SolnData.dispDotPrev     =  SolnData.dispDot;
    SolnData.dispDotDotPrev  =  SolnData.dispDotDot;

    dispTemp.resize(totalDOF);
    rhsVec2.resize(totalDOF);

    filecount=1;

    //Time loop
    while( (stepsCompleted <= (stepsMax+1) ) && (timeNow <= timeFinal) )
    {
        for(ii=0; ii<totalDOF; ii++)
        {
          jj = assyForSoln[ii];

          dispTemp(ii) = SolnData.disp(jj);
        }

        // calculate the residual vector
        rhsVec2 = solverEigen->mtx * dispTemp;

        // update the global residual vector
        rhsVec.setZero();
        for(ii=0; ii<totalDOF; ii++)
        {
          jj = assyForSoln[ii];

          rhsVec(jj) = -rhsVec2(ii);
        }

        //cout << " aaaaaaaaaaaaa " << endl;

        //Add specified nodal force 
        //////////////////////////////////////////

        //if(timeNow <= 4.0e-5)
        if(timeNow <= 20.0)
          fact = timeFunction[0].prop;
        else
          fact = 0.0;

        //cout << " fact=" << timeNow << '\t' << fact << endl;

        for(ii=0; ii<nodeForcesData.size(); ii++)
        {
          nn  = (int) (nodeForcesData[ii][0] - 1);
          dof = (int) (nodeForcesData[ii][1] - 1);
          specVal = nodeForcesData[ii][2];

          rhsVec(nn*ndof+dof) += (specVal*fact);
        }
        //cout << " bbbbbbbbbbb " << endl;
        //rNorm = rhsVec.norm();
        //printVector(rhsVec);

        //G-alpha scheme
        for(ii=0; ii<totalDOF; ii++)
        {
          jj = assyForSoln[ii];

          //Solver for acceleration
          SolnData.dispDotDot(jj) = (rhsVec(jj)-(1.0-am)*MassVector(jj)*SolnData.dispDotDotPrev(jj))/(am*MassVector(jj));

          //Compute displacement and velocity
          SolnData.disp(jj)    = SolnData.dispPrev(jj) + dt*SolnData.dispDotPrev(jj) + dt*dt*( beta*SolnData.dispDotDot(jj) + (0.5-beta)*SolnData.dispDotDotPrev(jj) );

	  SolnData.dispDot(jj) = SolnData.dispDotPrev(jj) + dt*( gamma*SolnData.dispDotDot(jj) + (1.0-gamma)*SolnData.dispDotDotPrev(jj) );

          if( abs(SolnData.disp(jj)) > 1.0e40 )
          {
            cerr << "displacement too high in explicit simulation..." << endl;
            cerr << "Simulation aborted...\n\n\n" << endl;
            exit(1);
          }
        }

        // Apply Dirichlet BCs
        //////////////////////////////////////////

        //cout << " applying BCs " << endl;

        fact = timeNow;

        for(ii=0; ii<DirichletBCs.size(); ii++)
        {
          nn  = (int) (DirichletBCs[ii][0]);
          dof = (int) (DirichletBCs[ii][1]);
          specVal = DirichletBCs[ii][2];

          jj = nn*ndof+dof;

          //cout << nn << '\t' << dof << '\t' << specVal << '\t' << timeNow << endl;
          SolnData.disp(jj) = specVal*fact;

          // update acceleration and velocity
          SolnData.dispDotDot(jj) = (SolnData.disp(jj)-SolnData.dispPrev(jj)-dt*SolnData.dispDotPrev(jj))/(beta*dt*dt) - (0.5/beta-1.0)*SolnData.dispDotDotPrev(jj);

          SolnData.dispDot(jj)    = SolnData.dispDotPrev(jj) + dt*( gamma*SolnData.dispDotDot(jj) + (1.0-gamma)*SolnData.dispDotDotPrev(jj) );
        }
        //printVector(SolnData.dispDot);


        SolnData.dispIncr = SolnData.disp;
        int  idd = SolnData.ElemProp[0]->id;
        if( (idd == 204) || (idd == 208) )
        {
          for(ee=0;ee<nElem_global;ee++)  // loop over all the elements
          {
            elems[ee]->solveForPressure();
          }
        }
        //cout << " aaaaaaaaaaaaa " << endl;
        //update variables at (n+af)
        updateGeometry();

        // write the data for every N time steps
        if( (stepsCompleted%1000 == 0) )
        {
          printf("\n\n\n stepsCompleted = %5d ;  timeNow = %12.8f \n", stepsCompleted, timeNow);
          writeNodalData();

          postProcess(4, 0, 1, false, 0, 0, resln);
          filecount++;
        }

        //store the variables
        SolnData.dispPrev2       =  SolnData.dispPrev;
        SolnData.dispPrev        =  SolnData.disp;
        SolnData.dispDotPrev     =  SolnData.dispDot;
        SolnData.dispDotDotPrev  =  SolnData.dispDotDot;

        stepsCompleted++;

        timeNow += dt;

        //updateIterStep();
    }

    fout.close();

    //if( !converged() )
      //return 1;

    return 0;
}
*/



int femSolidmechanics::solveExplicitStep(int* inputdata)
{
    // just in case this function is called by mistake while using the implicit solver
    if(IMPLICIT_SOLVER)
      return 0;

    // calculate the mass matrix
    calcMassMatrixForExplicitDynamics();

    auto  time1 = chrono::steady_clock::now();

    // find the solution
    //if(tis == 500)
      //solveExplicitStepFIC(inputdata);
    //else
      //solveExplicitStepNIC(inputdata);

    auto  time2 = chrono::steady_clock::now();
    auto  duration = chrono::duration_cast<chrono::milliseconds>(time2-time1).count();
    cout << " Time (milliseconds)   = " << duration << endl;

    return 0;
}










int femSolidmechanics::calcMassMatrixForExplicitDynamics()
{
    cout << " calcMassMatrixForExplicitDynamics()  STARTED" << endl;

    int  ee, ii, jj, nn, n1, n2, n3, n4;
    int idd = 0;// SolnData.ElemProp[0]->id;

    MatrixXd  Mlocal, Mlocal2;
    VectorXd  Flocal(100);

    MassVector.resize(nNode_global*ndof);
    MassVector.setZero();

    ForceVectorExternal.resize(nNode_global*ndof);
    ForceVectorExternal.setZero();

    for(ee=0;ee<nElem_global;ee++)  // loop over all the elements
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

    if( (idd == 211) || (idd == 212) || (idd == 213) || (idd == 214) )
    {
        MassVectorPres.resize(nNode_global);
        MassVectorPres.setZero();

        if(
          (idd == 211) || // Bezier TRIA6/1
          (idd == 213)    // Bezier TET10/1
          )
        {
          for(ee=0;ee<nElem_global;ee++)
          {
            elems[ee]->calcMassMatrixPressure(Mlocal2, true);

            MassVectorPres(ee) += Mlocal2(0, 0);
          }
        }
        else if( (idd == 212) || (idd == 214) )
        {
          if(idd == 212) // Bezier TRIA6/3
            nn = 3;
          else if(idd == 214) // Bezier TET10/4
            nn = 4;

          for(ee=0;ee<nElem_global;ee++)
          {
            elems[ee]->calcMassMatrixPressure(Mlocal2, true);

            //printMatrix(Mlocal2);

            for(ii=0; ii<nn; ii++)
            {
              MassVectorPres(elems[ee]->nodeNums[ii]) += Mlocal2(ii, ii);
            }
          }
        }
        else
        {
          cerr << " Wrong element type in 'femSolidmechanics::calcMassMatrixForExplicitDynamics()' " << endl;
          //exit(999);
        }
    }

    cout << " calcMassMatrixForExplicitDynamics() ENDED" << endl;

    //printVector(MassVector);
    //printVector(MassVectorPres);

    return 0;
}


/*
int  femSolidmechanics::ModalAnalysis(int nModes, bool flag, double fact)
{
    //elementDiffStiffTest(0.01, 1, 6, 6, true);

    cout << " femSolidmechanics::ModalAnalysis() ... STARTED " << endl;
    flag = false;
    flag = true;
    cout << " flag = " << flag << endl;

    MatrixXd  Kglobal(totalDOF, totalDOF), Mglobal(totalDOF, totalDOF), eigen_vectors;
    MatrixXd  Kuu, Mlocal;
    VectorXd  eigen_values, eigvec, Flocal, vecTemp;

    int ee, ii, jj, nn, rr, cc;


    Kglobal.setZero();
    Mglobal.setZero();
    for(ee=0;ee<nElem_global;ee++)  // loop over all the elements
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

      SolnData.disp.setZero();
      for(ii=0; ii<totalDOF; ii++)
      {
        SolnData.disp[assyForSoln[ii]] = eigvec[ii];
      }

      postProcess(4, 7, 1, false, 0, 0, resln);
      filecount++;
    }

    cout << " femSolidmechanics::ModalAnalysis() ... ENDED " << endl;

    return 0;
}
*/



//
int  femSolidmechanics::ModalAnalysis(int nModes, bool flag, double fact)
{
    cout << " femSolidmechanics::ModalAnalysis() ... STARTED " << endl;

    // to compute inf-sup number

    int ee, ii, jj, nn, rr, cc, size1, size2;

    MatrixXf  Kglobal, Mglobal, eigen_vectors;
    MatrixXd  Kuu(12,12), Kup(12,1), Kpp(1,1);
    VectorXf  eigen_values, eigvec, Flocal, vecTemp;
    MatrixXf  Vmat, Qmat, Bmat;

    int idd = 0;//SolnData.ElemProp[0]->id;
    cout <<  "idd = " <<  idd <<  endl;

    if( (idd == 204) || (idd == 208) )
    {
      Vmat.resize(dispDOF, dispDOF);
      Bmat.resize(dispDOF, nElem_global);
      Qmat.resize(dispDOF, dispDOF);

      Kglobal.resize(dispDOF, dispDOF);
      Mglobal.resize(nElem_global, nElem_global);

      Kglobal.setZero();
      Mglobal.setZero();
      Vmat.setZero();
      Bmat.setZero();
      Qmat.setZero();

      /////////////////////////////////////////
      // compute and assemble matrices
      /////////////////////////////////////////

      for(ee=0;ee<nElem_global;ee++)                               // loop over all the elements
      {
        //cout << "       elem... : " << (e+1) << endl;

        elems[ee]->toComputeInfSupCondition(Kuu, Kup, Kpp);

        //printVector(elems[ee]->forAssyVec);
        size1=elems[ee]->forAssyVec.size();

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
          }
        }

        // Kup
        for(ii=0; ii<size1; ii++)
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
    }
    else if( (idd == 215) || (idd == 216) )
    {
      presDOF = 2*nElem_global;

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

      for(ee=0;ee<nElem_global;ee++)                               // loop over all the elements
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
        //size2 = elems[ee]->forAssyVecPres.size();
        size2 = 2;

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

            Bmat(rr, 2*ee)   += Kup(ii, 0);
            Bmat(rr, 2*ee+1) += Kup(ii, 1);
          }
        }

        // Kpp
        ii = 2*ee;

        Mglobal(ii,   ii)   += Kpp(0,0);
        Mglobal(ii,   ii+1) += Kpp(0,1);
        Mglobal(ii+1, ii)   += Kpp(1,0);
        Mglobal(ii+1, ii+1) += Kpp(1,1);
      }
    }
    else
    {
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

      for(ee=0;ee<nElem_global;ee++)                               // loop over all the elements
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
        //size2 = elems[ee]->forAssyVecPres.size();
        size2 = 2;

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
    }


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


int femSolidmechanics::elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt)
{
  if(elnum > nElem_global)
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

  return 0;
}




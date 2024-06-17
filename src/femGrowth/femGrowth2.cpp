
#include "SolverPardisoEigen.h"
#include "femGrowth.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "ElementBase.h"
#include "SolutionData.h"
#include "FunctionsProgram.h"
#include <chrono>

#include "SolverPetsc.h"
#include "SolverPardisoPetsc.h"

extern MyTime myTime;
extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern bool debug;



int femGrowth::setSolver(int slv)
{
    if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femGrowth::setSolver ... STARTED \n");

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

    //setSolverDataForFullyImplicit();

    if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femGrowth::setSolver ... ENDED \n");

    return 0;
}




int femGrowth::prepareMatrixPattern()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n     femGrowth::prepareMatrixPattern()  .... STARTED ...\n\n");

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

    PetscPrintf(MPI_COMM_WORLD, "\n     femGrowth::prepareMatrixPattern()  .... FINISHED ...\n\n");

/*
    prepareDataForPressure();
    prepareDataForMagneticPotential();
    //prepareDataForTemperature();

    ntotdofs_global = dispDOF + presDOF + mpotDOF + tempDOF;

    cout << " nElem_global    = " << nElem_global    << endl;
    cout << " nNode    = " << nNode    << endl;
    cout << " npElem   = " << npElem   << endl;
    cout << " ndof     = " << ndof     << endl;
    cout << " dispDOF  = " << dispDOF  << endl;
    cout << " presDOF  = " << presDOF  << endl;
    cout << " mpotDOF  = " << mpotDOF  << endl;
    cout << " tempDOF  = " << tempDOF  << endl;
    cout << " ntotdofs_global = " << ntotdofs_global << endl;

        if(mpotDOF > 0)
        {
          cout << "Adding matrix entries for the displacement-electricpotential coupling " << endl;

          int  *ttF,  *ttU;
          int  sizeF, sizeU, row, col, offset=dispDOF+presDOF;

          for(ee=0;ee<nElem_global;ee++)
          {
            sizeF = elems[ee]->forAssyVecMpot.size();
            ttF   = &(elems[ee]->forAssyVecMpot[0]);

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
          int  sizeT, sizeU, row, col, offset=dispDOF+presDOF+mpotDOF;

          for(ee=0;ee<nElem_global;ee++)
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
*/

    return 0;
}



// set the off-diagonal terms for the solver
int femGrowth::setSolverDataForSemiImplicit()
{
/*
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


    for(ee=0; ee<nElem_global; ee++)
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

    rhsElecField.resize(mpotDOF);
    spmtxElecField.setZero();
    spmtxElecField.resize(mpotDOF, mpotDOF);
    spmtxElecField.reserve(ceil(mpotDOF*mpotDOF*0.1));
    jj = ceil(mpotDOF*0.2);
    cout << " jj = " << jj << endl;

    for(ii=0; ii<mpotDOF; ii++)
      nnzVec(ii) = jj;
    spmtxElecField.reserve(nnzVec);


    int  *ttF, *ttU, sizeF;

    for(ee=0;ee<nElem_global;ee++)
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
*/
    return 0;
}


int femGrowth::solveFullyImplicit()
{
    if(ntotdofs_global == 0)
      return 0;

    if( (SolnData.timeIntegrationScheme != "STEADY") || (NONLINEAR_SOLVER_TYPE == SOLVER_TYPE_NEWTONRAPHSON) )
      solveWithNewtonRaphson();
    else if(NONLINEAR_SOLVER_TYPE == SOLVER_TYPE_ARCLENGTH)
      solveWithArclength();
    else
    {
        cout << " femGrowth::solveFullyImplicit() ... Algorithm type is not available " << endl;
    }

    return 0;
}





int femGrowth::solveWithNewtonRaphson()
{
    cout << "\femGrowth::solveWithNewtonRaphson \n\n\n" << endl;

    int  stepsCompleted=1, err = 0;

    solverPetsc->Fext.resize(nNode_global*ndof);

    //setInitialConditions();
    //readResult();
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
        PetscPrintf(MPI_COMM_WORLD, " Maximum steps        =  %d  \n", stepsMax);
        PetscPrintf(MPI_COMM_WORLD, " Time step size       =  %f  \n", myTime.dt);
        PetscPrintf(MPI_COMM_WORLD, " Current time         =  %f  \n", myTime.cur);
        PetscPrintf(MPI_COMM_WORLD, " Final time           =  %f  \n", timeFinal);
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

            PetscPrintf(MPI_COMM_WORLD, " femGrowth ...  %3d \t %11.4e \n", (iter), rhsNorm);

            // check for convergence and divergence of the iterations
            if( converged() )
            {
              PetscPrintf(MPI_COMM_WORLD, "\n femGrowth ...  Iterations CONVERGED \n\n\n");

              convergedFlag = true;

              break;
            }
            else if( (iter > 3) && diverging(1.0e7) )
            {
              PetscPrintf(MPI_COMM_WORLD, " femGrowth ...  Iterations are diverging. NR loop is terminated. \n\n\n");
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

    writeResult();

    PetscPrintf(MPI_COMM_WORLD, "\n\n\n Simulation reached the specified final time or maximum steps specified ... \n\n\n");

    return err;
}




int femGrowth::solveWithArclength()
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






int femGrowth::calcStiffnessAndResidual()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femGrowth::calcStiffnessAndResidual ...STARTED \n\n");}

    //int  ee, ii, jj, nsize, row;

    MatrixXd  Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp;
    VectorXd  FlocalU, FlocalF, FlocalP;
    //MatrixXd  Kuu(30, 30), Kuf(30, 10), Kfu(10,30), Kff(10,10), Kup, Kpu, Kpp;
    //VectorXd  FlocalU(30), FlocalF(10), FlocalP(10);


    SolnData.reac.setZero();

    solverPetsc->zeroMtx();

    if(COUPLED_PROBLEM)
    {
      int mpot_offset = dispDOF+presDOF;
      // loop over all the elements
      for(int ee=0;ee<nElem_global;ee++)
      {
        //cout << "       elem... : " << (ee+1) << endl;

        try
        {
            elems[ee]->calcStiffnessAndResidualMixed(Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp, FlocalU, FlocalF, FlocalP);
            //cout << " AAAAAAAAAAA " << endl;
        }
        catch(runtime_error& err)
        {
            cerr << err.what() << endl;
            throw runtime_error("Negative Jacobian encountered");
        }

        elems[ee]->calcLoadVector(FlocalU);

        //if(firstIteration)
          //elems[ee]->applyDirichletBCs3field(1, Kuu, Kuf, Kfu, Kff, Kup, Kpu, Kpp, FlocalU, FlocalF, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecMpot, elems[ee]->forAssyVecPres, SolnData.dispAapplied, SolnData.mpotApplied, SolnData.presApplied);

        solverPetsc->assembleMatrixAndVector3field(dispDOF, mpot_offset, Kuu, Kup, Kpu, Kpp, Kuf, Kfu, Kff, FlocalU, FlocalP, FlocalF, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres, elems[ee]->forAssyVecMpot);
      }
    }
    else
    {
      if(MIXED_ELEMENT)
      {
        // loop over all the elements
        for(int ee=0;ee<nElem_global;ee++)
        {
          //cout << "       elem... : " << (ee+1) << endl;
          elems[ee]->calcStiffnessAndResidualMixed(Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP);

          if(firstIteration)
            elems[ee]->applyDirichletBCs2field(1, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres, SolnData.dispApplied, SolnData.presApplied);

          solverPetsc->assembleMatrixAndVector2field(dispDOF, Kuu, Kup, Kpu, Kpp, FlocalU, FlocalP, elems[ee]->forAssyVec, elems[ee]->forAssyVecPres);
        }
      }
      else
      {
        // loop over all the elements
        for(int ee=0;ee<nElem_global;ee++)
        {
          //cout << "       elem... : " << (ee+1) << endl;
          elems[ee]->calcStiffnessAndResidual(Kuu, FlocalU);
          elems[ee]->calcLoadVector(FlocalU);

          if(firstIteration)
            elems[ee]->applyDirichletBCs(Kuu, FlocalU);

          solverPetsc->assembleMatrixAndVector(0, 0, elems[ee]->forAssyVec, elems[ee]->forAssyVec, Kuu, FlocalU);
        }
      }
    }

    //cout << " BBBBBBBBBBBBB " << endl;
    //printf("\n solverPetsc->rhsVec norm = %12.6E \n", solverPetsc->rhsVec.norm());
    //cout << " RHS " << endl;        printVector(solverPetsc->rhsVec); printf("\n\n\n");

    solverPetsc->currentStatus = ASSEMBLY_OK;

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femGrowth::calcStiffnessAndResidual ... ENDED \n\n");}

    return 0;
}



int femGrowth::factoriseSolveAndUpdate()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femGrowth::factoriseSolveAndUpdate ... STARTED \n\n");}

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

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "     femGrowth::factoriseSolveAndUpdate ... ENDED \n\n");}

/*
    if(debug) {cout << "     femGrowth::factoriseSolveAndUpdate ... STARTED \n\n";}

    time_t tstart, tend;

    //cout << " RHS " << endl;        printVector(solverEigen->rhsVec); printf("\n\n\n");
    //tstart = time(0);
    //for(int ii=dispDOF-100; ii<ntotdofs_global; ii++)
      //cout << ii << '\t' << solverEigen->rhsVec[ii] << endl;

    SolnData.dispIncr.setZero();

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

    if(debug) {cout << "     matrix solution STARTED \n\n";}

    solverPetsc->factoriseAndSolve();

    if(debug) {cout << "     matrix solution DONE \n\n";}

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
    // only for the constant pressure elements, Quad4/1, Hex8/1

    int  idd = SolnData.ElemProp[0]->id;
    if( (idd == 2010) || (idd == 2060) )
    {
      for(int ee=0;ee<nElem_global;ee++)  // loop over all the elements
      {
        elems[ee]->solveForPressure();
      }
    }
*/

    return 0;
}















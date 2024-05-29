
#include "femINSstabilised.h"
#include "util.h"
#include "ElementBase.h"
#include <chrono>
#include "KimMoinFlow.h"
#include "metis.h"
#include "QuadratureUtil.h"
#include "BasisFunctionsLagrange.h"
#include "stabilisationRoutines.h"
#include "MyTime.h"
#include "TimeFunction.h"



extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime                 myTime;



int femINSstabilised::setSolver(int slv, int *parm, bool cIO)
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n     femINSstabilised::setSolver()  .... STARTED ...\n\n");

    Eigen::setNbThreads(0);

    solverPetsc = make_unique<SolverPetsc>();

    computerTimePattern = MPI_Wtime();
    prepareMatrixPattern();
    computerTimePattern = MPI_Wtime() - computerTimePattern;

    PetscPrintf(MPI_COMM_WORLD, "\n\n     femINSstabilised::setSolver()  .... ENDED ...\n\n");

    return 0;
}




int femINSstabilised::prepareMatrixPattern()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n     femINSstabilised::prepareMatrixPattern()  .... STARTED ...\n\n");

    int  size1, size2, row, col, npElem, row_start, row_end;
    int  tempDOF, domTemp, ind;
    int  r, c, r1, c1, count=0, count1=0, count2=0, ii, jj, ee, dd, ind1, ind2, nsize;
    int  val1, val2, nnz, nnz_max_row, n1, n2, kk, e1, e2, ll, aa, bb, nn, dof;
    int  side, start1, start2, nr1, nr2, count_diag, count_offdiag, tempInt;

    PetscInt  *colTemp;
    PetscScalar  *arrayTemp;
    double  tstart, tend;

    int ndomains = n_mpi_procs, subdomain=0;


    /////////////////////////////////////////////////////////////
    //
    // prepare the matrix pattern
    /////////////////////////////////////////////////////////////

    vector<vector<int> >     NodeDofArray;
    vector<vector<bool> >    NodeType;

    // set sizes of some data arrays
    vector<bool>  vecBoolTempFalse(ndof, false);
    NodeType.resize(mesh->nNode_global, vecBoolTempFalse);

    vector<int>  vecIntTempM1(ndof, -1);
    NodeDofArray.resize(mesh->nNode_global, vecIntTempM1);

    // fix the specified Dirichlet BCs
    dofs_specified.clear();
    setSpecifiedDOFs(NodeType);


    ntotdofs_global = 0;
    for(ii=0;ii<mesh->nNode_global;++ii)
    {
      for(jj=0;jj<ndof;++jj)
      {
        if(NodeType[ii][jj] == false)
        {
          NodeDofArray[ii][jj] = ntotdofs_global++;
        }
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    ntotdofs_local = ntotdofs_global;
    row_start      =  0;
    row_end        = ntotdofs_global-1;

    if(n_mpi_procs > 1)
    {
      // compute first and last row indices of the rows owned by the local processor
      row_start  =  1e9;
      row_end    = -1e9;
      ntotdofs_local = 0;
      for(ii=mesh->node_start; ii<=mesh->node_end; ii++)
      {
        for(jj=0; jj<ndof; jj++)
        {
          if(NodeType[ii][jj] == false)
          {
            ind = NodeDofArray[ii][jj];
            row_start  = min(row_start, ind);
            row_end    = max(row_end,   ind);
            ntotdofs_local++;
          }
        }
      }

      cout << "ntotdofs_local = " << ntotdofs_local << '\t' << ntotdofs_global << '\t' << this_mpi_proc << endl;

      cout << "row_start  = " << row_start  << '\t' << row_end     << '\t' << this_mpi_proc << endl;

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

    forAssyVecAll.resize(mesh->nElem_global);
    globalDOFnumsAll.resize(mesh->nElem_global);

    vector<int> nodeNums, forAssyVec, globalDOFnums;
    for(ee=0; ee<mesh->nElem_global; ++ee)
    {
      //if(elems[ee]->getSubdomainId() == this_mpi_proc)
      if(elem_proc_id[ee] == this_mpi_proc)
      {
        //elems[ee]->prepareElemData(nodeCoordsOrig);
        //npElem = elems[ee]->nodeNums.size();

        nodeNums = mesh->elemConn[ee];
        npElem = nodeNums.size();

        nsize = ndof*npElem;

        //elems[ee]->forAssyVec.resize(nsize);
        forAssyVec.resize(nsize);
        globalDOFnums.resize(nsize);

        for(ii=0; ii<npElem; ++ii)
        {
          n1 = ndof*ii;
          n2 = ndof*nodeNums[ii];

          //kk = elems[ee]->nodeNums[ii];
          kk = nodeNums[ii];

          for(dof=0; dof<ndof; ++dof)
          {
            globalDOFnums[n1+dof] = n2+dof;

            //elems[ee]->forAssyVec[n1+dof] = NodeDofArrayNew[kk][dof];
            forAssyVec[n1+dof] = NodeDofArray[kk][dof];
          }
        }

        forAssyVecAll[ee]    = forAssyVec;
        globalDOFnumsAll[ee] = globalDOFnums;
      }
    }

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    assyForSoln.resize(ntotdofs_global);
    ind = 0;
    for(ii=0;ii<mesh->nNode_global;++ii)
    {
      for(jj=0;jj<ndof;++jj)
      {
        if(NodeDofArray[ii][jj] != -1)
        {
          assyForSoln[ind++] = ii*ndof+jj;
        }
      }
    }

    PetscPrintf(MPI_COMM_WORLD, "\n\n Element DOF values initialised \n\n");
    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    PetscPrintf(MPI_COMM_WORLD, " Total DOF   = %d \t %d \n", ntotdofs_local, ntotdofs_global);


    vector<vector<int> > forAssyMat;
    vector<int>::iterator it;
    //vector<set<int> > forAssyMat;
    //set<int>::iterator it;

    forAssyMat.resize(ntotdofs_global);

    //for(ii=row_start; ii<=row_end; ii++)
      //forAssyMat[ii].reserve(500);

    for(ee=0; ee<mesh->nElem_global; ee++)
    {
      if(elem_proc_id[ee] == this_mpi_proc)
      {
        forAssyVec = forAssyVecAll[ee];
        nsize = forAssyVec.size();

        for(ii=0;ii<nsize;ii++)
        {
            r = forAssyVec[ii];

            if(r != -1)
            {
              if(r >= row_start && r <= row_end)
              {
              for(jj=0;jj<nsize;jj++)
              {
                if(forAssyVec[jj] != -1)
                {
                  //printf("ii.... %5d \t %5d \t %5d \t %5d \n",ii, jj, r, forAssyVec[jj]);
                    forAssyMat[r].push_back(forAssyVec[jj]);
                }
              }
              }
            }
        }
      }
    }

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);
    PetscPrintf(MPI_COMM_WORLD, "\n\n Preparing matrix pattern DONE \n\n");


    PetscInt  *diag_nnz, *offdiag_nnz;

    errpetsc  = PetscMalloc1(ntotdofs_local,  &diag_nnz);CHKERRQ(errpetsc);
    errpetsc  = PetscMalloc1(ntotdofs_local,  &offdiag_nnz);CHKERRQ(errpetsc);


    kk = 0;
    nnz_max_row = 0;
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

    //Create parallel matrix, specifying only its global dimensions.
    //When using MatCreate(), the matrix format can be specified at
    //runtime. Also, the parallel partitioning of the matrix is
    //determined by Petsc at runtime.
    //Performance tuning note: For problems of substantial size,
    //preallocation of matrix memory is crucial for attaining good
    //performance. See the matrix chapter of the users manual for details.

    PetscPrintf(MPI_COMM_WORLD, " Initialise the Matrix pattern \n", errpetsc);


    PetscScalar  Klocal[nnz_max_row];
    PetscInt  rows[1];
    size1 = 1;
    for(ii=row_start; ii<=row_end; ii++)
    {
      rows[0] = ii;

      forAssyVec = forAssyMat[ii];
      size2 = forAssyVec.size();

      errpetsc = MatSetValues(solverPetsc->mtx, size1, rows, size2, &forAssyVec[0], Klocal, INSERT_VALUES);
    }

    /*
    npElem = 27;
    ind = npElem*ndof;
    ind = ind*ind;
    PetscScalar  Klocal[ind];
    for(ii=0; ii<ind; ii++)  Klocal[ii] = 0.0;
    for(ee=0; ee<mesh->nElem_global; ee++)
    {
      if(elem_proc_id[ee] == this_mpi_proc)
      {
        forAssyVec = forAssyVecAll[ee];
        size1 = forAssyVec.size();

        errpetsc = MatSetValues(solverPetsc->mtx, size1, &forAssyVec[0], size1, &forAssyVec[0], Klocal, INSERT_VALUES);
      }
    } //for(ee=0;)
    */

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    //MatAssemblyBegin(solverPetsc->mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
    //MatAssemblyEnd(solverPetsc->mtx, MAT_FINAL_ASSEMBLY);//CHKERRQ(ierr);
    //MatZeroEntries(solverPetsc->mtx);
    //errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    // Create reaction vector
    errpetsc = VecCreate(PETSC_COMM_WORLD, &(solverPetsc->reacVec));
    CHKERRQ(errpetsc);

    ind1 = mesh->nNode_owned*ndof;
    ind2 = mesh->nNode_global*ndof;

    errpetsc = VecSetSizes(solverPetsc->reacVec, ind1, ind2);
    CHKERRQ(errpetsc);

    errpetsc = VecSetFromOptions(solverPetsc->reacVec);
    CHKERRQ(errpetsc);

    //VecGetSize(solverPetsc->reacVec, &size1);
    //cout << " size = " << size1 << endl;

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);


    solverPetsc->currentStatus = PATTERN_OK;

    errpetsc  = PetscFree(diag_nnz);   CHKERRQ(errpetsc);
    errpetsc  = PetscFree(offdiag_nnz);   CHKERRQ(errpetsc);

    for(ii=0; ii<NodeDofArray.size(); ii++)
    {
      NodeDofArray[ii].clear();
    }
    NodeDofArray.clear();

    PetscPrintf(MPI_COMM_WORLD, " femINSstabilised::prepareMatrixPattern()  .... FINISHED. Took %f  milliseconds \n", (tend-tstart)*1000);

    return 0;
}







int  femINSstabilised::solveFullyImplicit()
{
    PetscPrintf(MPI_COMM_WORLD, " Solving with the Fully-Implicit Scheme \n\n\n");

    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////


    int  stepsCompleted=1, ii, jj;

    PetscScalar *arrayTempSoln;
    Vec            vecseq;
    VecScatter     ctx;

    VecScatterCreateToAll(solverPetsc->solnVec, &ctx, &vecseq);


    setInitialConditions();
    postProcess();
    writeOutputData();

    convergedFlagPrev = convergedFlag = false;

    computerTimeTimeLoop = MPI_Wtime();

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


        MPI_Barrier(MPI_COMM_WORLD);

        PetscPrintf(MPI_COMM_WORLD, " Iteration loop \n");
        for(int iter=1; iter<=iterationsMax; iter++)
        {
            firstIteration = (iter == 1);

            // Compute the velocity and acceleration and the respective values at n+af, n+am
            updateIterStep();

            PetscPrintf(MPI_COMM_WORLD, "\n\n Zeroing matrix entries \n\n");
            solverPetsc->zeroMtx();

            MPI_Barrier(MPI_COMM_WORLD);

            timerVal = MPI_Wtime();

            calcStiffnessAndResidual(iter);

            timerVal = MPI_Wtime() - timerVal;
            computerTimeAssembly += timerVal;
            PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time for matrix assembly = %f seconds \n\n", timerVal );

            MPI_Barrier(MPI_COMM_WORLD);

            //PetscPrintf(MPI_COMM_WORLD, "\n Adding boundary conditions \n");

            // add boundary conditions
            if(firstIteration)
            {
              int nDBC = dofs_specified.size();
              for(ii=0; ii<nDBC; ++ii)
              {
                jj = dofs_specified[ii];

                soln[jj]  += solnApplied[jj];
              }
            }

            MPI_Barrier(MPI_COMM_WORLD);

            rhsNormPrev = rhsNorm;

            VecAssemblyBegin(solverPetsc->rhsVec);
            VecAssemblyEnd(solverPetsc->rhsVec);

            VecNorm(solverPetsc->rhsVec, NORM_2, &rhsNorm);
            solverPetsc->currentStatus = ASSEMBLY_OK;

            MPI_Barrier(MPI_COMM_WORLD);

            //VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

            PetscPrintf(MPI_COMM_WORLD, " Iteration = %d  \t RHS norm = %E \n", iter, rhsNorm);

            if(rhsNorm < conv_tol)
            {
                PetscPrintf(MPI_COMM_WORLD, " Solution converged below the specified tolerance. \n\n");

                convergedFlag = true;

                break;
            }
            else
            {
                PetscPrintf(MPI_COMM_WORLD, "Assembly done. Solving the matrix system. \n");

                timerVal = MPI_Wtime();
                if( solverPetsc->factoriseAndSolve() )
                {
                  PetscPrintf(MPI_COMM_WORLD, " PETSc solver not converged. \n\n");
                  return -1;
                }
                timerVal = MPI_Wtime() - timerVal;
                computerTimeSolver += timerVal;
                PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time for PETSc solver = %f seconds \n\n", timerVal );

                /////////////////////////////////////////////////////////////////////////////
                // get the solution vector onto all the processors
                /////////////////////////////////////////////////////////////////////////////

                VecScatterBegin(ctx, solverPetsc->solnVec, vecseq, INSERT_VALUES, SCATTER_FORWARD);
                VecScatterEnd(ctx,   solverPetsc->solnVec, vecseq, INSERT_VALUES, SCATTER_FORWARD);

                VecGetArray(vecseq, &arrayTempSoln);

                // update solution vector
                for(ii=0; ii<ntotdofs_global; ii++)
                {
                  soln[assyForSoln[ii]]   +=  arrayTempSoln[ii];
                }

                VecRestoreArray(vecseq, &arrayTempSoln);
            }
        } //Iteration Loop


        // if the residual is converged, then save the DOFs vectors
        if( convergedFlag )
        {
            PetscPrintf(MPI_COMM_WORLD, " Postprocessing... \n\n");

            timerVal = MPI_Wtime();

            if( stepsCompleted%outputFreq == 0 )
              postProcess();

            MPI_Barrier(MPI_COMM_WORLD);

            writeOutputData();

            computerTimePostprocess += (MPI_Wtime() - timerVal);

            saveSolution();

            myTime.stck();

            stepsCompleted++;
        }
        else
        {
            myTime.cut();

            reset();
        }
 
    } //Time loop


    VecScatterDestroy(&ctx);
    VecDestroy(&vecseq);

    computerTimeTimeLoop = MPI_Wtime() - computerTimeTimeLoop;

    PetscPrintf(MPI_COMM_WORLD, "\n\n Simulation reached the specified final time or maximum steps specified ... \n\n\n");

    return 0;
}





int femINSstabilised::solveTimeStep()
{
    int  ii, jj;


    PetscScalar *arrayTempSoln;
    Vec            vecseq;
    VecScatter     ctx;

    VecScatterCreateToAll(solverPetsc->solnVec, &ctx, &vecseq);


        convergedFlagPrev = convergedFlag;
        convergedFlag = false;

        rhsNormPrev = rhsNorm = -1.0;


        MPI_Barrier(MPI_COMM_WORLD);

        PetscPrintf(MPI_COMM_WORLD, " Iteration loop \n");
        for(int iter=1; iter<=iterationsMax; iter++)
        {
            firstIteration = (iter == 1);

            // Compute the velocity and acceleration and the respective values at n+af, n+am
            updateIterStep();

            PetscPrintf(MPI_COMM_WORLD, "\n solverPetsc->zeroMtx() \n");

            solverPetsc->zeroMtx();

            MPI_Barrier(MPI_COMM_WORLD);

            timerVal = MPI_Wtime();

            calcStiffnessAndResidual(iter);

            timerVal = MPI_Wtime() - timerVal;
            computerTimeAssembly += timerVal;
            PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time for matrix assembly = %f seconds \n\n", timerVal );

            MPI_Barrier(MPI_COMM_WORLD);

            //PetscPrintf(MPI_COMM_WORLD, "\n Adding boundary conditions \n");

            // add boundary conditions
            if(firstIteration)
            {
              int nDBC = dofs_specified.size();
              for(ii=0; ii<nDBC; ++ii)
              {
                jj = dofs_specified[ii];

                soln[jj]  += solnApplied[jj];
              }
            }

            MPI_Barrier(MPI_COMM_WORLD);

            rhsNormPrev = rhsNorm;

            VecAssemblyBegin(solverPetsc->rhsVec);
            VecAssemblyEnd(solverPetsc->rhsVec);

            VecNorm(solverPetsc->rhsVec, NORM_2, &rhsNorm);
            solverPetsc->currentStatus = ASSEMBLY_OK;

            MPI_Barrier(MPI_COMM_WORLD);

            //VecView(solverPetsc->rhsVec, PETSC_VIEWER_STDOUT_WORLD);

            PetscPrintf(MPI_COMM_WORLD, " Iteration = %d  \t RHS norm = %E \n", iter, rhsNorm);

            if(rhsNorm < conv_tol)
            {
                PetscPrintf(MPI_COMM_WORLD, " Navier-Stokes : Solution converged below the specified tolerance. \n\n");

                convergedFlag = true;

                break;
            }
            else
            {
                PetscPrintf(MPI_COMM_WORLD, "Assembly done. Solving the matrix system. \n");

                timerVal = MPI_Wtime();
                if( solverPetsc->factoriseAndSolve() )
                {
                  PetscPrintf(MPI_COMM_WORLD, " PETSc solver not converged. \n\n");
                  return -1;
                }
                timerVal = MPI_Wtime() - timerVal;
                computerTimeSolver += timerVal;
                PetscPrintf(MPI_COMM_WORLD, "\n\n Elapsed time for PETSc solver = %f seconds \n\n", timerVal );

                /////////////////////////////////////////////////////////////////////////////
                // get the solution vector onto all the processors
                /////////////////////////////////////////////////////////////////////////////

                VecScatterBegin(ctx, solverPetsc->solnVec, vecseq, INSERT_VALUES, SCATTER_FORWARD);
                VecScatterEnd(ctx,   solverPetsc->solnVec, vecseq, INSERT_VALUES, SCATTER_FORWARD);

                VecGetArray(vecseq, &arrayTempSoln);

                // update solution vector
                for(ii=0; ii<ntotdofs_global; ii++)
                {
                  soln[assyForSoln[ii]]   +=  arrayTempSoln[ii];
                }

                VecRestoreArray(vecseq, &arrayTempSoln);
            }
        } //Iteration Loop

    VecScatterDestroy(&ctx);
    VecDestroy(&vecseq);

    return 0;
}









int femINSstabilised::calcForcesOnBoundaries()
{


    return 0;
}





int femINSstabilised::calcStiffnessAndResidual(int iter)
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n femINSstabilised::calcStiffnessAndResidual ... STARTED \n\n");

    if(SCHEME_TYPE == 1)
    {
        if(mesh->getNDIM() == 2)
          return calcStiffnessAndResidualImplicit2D(iter);
        else
          return calcStiffnessAndResidualImplicit3D(iter);
    }
    else if(SCHEME_TYPE == 2)
    {
        if(mesh->getNDIM() == 2)
          return calcStiffnessAndResidualLinearised2D(iter);
        else
          return calcStiffnessAndResidualLinearised3D(iter);
    }
    else if(SCHEME_TYPE == 3)
    {
        if(mesh->getNDIM() == 2)
          return calcStiffnessAndResidualExtrapolated2D(iter);
        else
          return calcStiffnessAndResidualExtrapolated2D(iter);
    }

    PetscPrintf(MPI_COMM_WORLD, "\n\n femINSstabilised::calcStiffnessAndResidual ... ENDED \n\n");

    return 0;
}





int femINSstabilised::calcStiffnessAndResidualImplicit2D(int iter)
{

    return 0;
}




int femINSstabilised::calcStiffnessAndResidualExtrapolated2D(int iter)
{

    return 0;
}



int femINSstabilised::calcStiffnessAndResidualImplicit3D(int iter)
{

    return 0;
}




int femINSstabilised::calcStiffnessAndResidualExtrapolated3D(int iter)
{

    return 0;
}





int femINSstabilised::calcStiffnessAndResidualLinearised2D(int iter)
{
    int  aa, ee, ii, jj, row, col, ind, size1, size2, npElem;
    int  gp, nGP, TI, TIp1, TIp2, TJ, TJp1, TJp2, degree;
    ElementShape  elShape;

    double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, Jac;
    double  pres, Da, Db, rad, urdr, urdr2, tau[3], CI=4.0;
    double  fact, fact1, fact2, elemVol, charlen, param[2];
    double  xNode[50], yNode[50];

    VectorXd  res(3), res2(2), dp(2), Du(2), vel(2), velPrev(2), velDot(2), force(2), gradTvel(2), rStab(3);
    MatrixXd  Dj(2,3), grad(2,2), gradN(2,2), stress(2,2), gradPrev(2,2);
    Dj.setZero();
    VectorXd  velTemp(3);


    double  rho = fluidProperties[0];
    double  mu  = fluidProperties[1];
    double  af  = td[2];
    double  acceFact = td[8];

    double  gamma = 1.96, eps = 0.02, WellHeight=0.25, phi=0.0, eta=0.0, dphi[2];
    double  kappa = 3.0*gamma/(4.0*sqrt(2*WellHeight)*eps);

    bool axsy = false;

    vector<int> nodeNums, forAssyVec, globalDOFnums;
    vector<double>  gpts1, gpts2, gwts, elemVolGP;

    //VectorXd  Flocal(ind);
    //MatrixXd  Klocal(ind, ind);
    MatrixXd  matG(3,3);    matG.setZero();
    vector<VectorXd>  Nv, dNvdx, dNvdy;


    // loop over elements and compute matrix and residual
    PetscPrintf(MPI_COMM_WORLD, "\n Element loop \n\n");
    for(ee=0; ee<mesh->getnElemGlobal(); ee++)
    {
    if(elem_proc_id[ee] == this_mpi_proc)
    {
        //cout << " ee = " << ee << endl;

        // calculate stiffness and residual
        //
        //
        nodeNums      = mesh->elemConn[ee];
        forAssyVec    = forAssyVecAll[ee];
        globalDOFnums = globalDOFnumsAll[ee];
        
        //printVector(nodeNums);
        //printVector(globalDOFnums);
        //printVector(forAssyVec);

        npElem = nodeNums.size();

        for(ii=0;ii<npElem;ii++)
        {
          xNode[ii] = mesh->nodeCoordsOrig[node_map_get_old[nodeNums[ii]]][0];
          yNode[ii] = mesh->nodeCoordsOrig[node_map_get_old[nodeNums[ii]]][1];
        }

        if(npElem == 3)
        {
            nGP     = 1;
            elShape = ELEM_SHAPE_TRIA;
            degree  = 1;
            getGaussPointsTriangle(nGP, gpts1, gpts2, gwts);
        }
        else if(npElem == 6)
        {
            nGP     = 3;
            elShape = ELEM_SHAPE_TRIA;
            degree  = 2;
            getGaussPointsTriangle(nGP, gpts1, gpts2, gwts);
        }
        else if(npElem == 4)
        {
            nGP     = 4;
            elShape = ELEM_SHAPE_QUAD;
            degree  = 1;
            getGaussPointsQuad(nGP, gpts1, gpts2, gwts);
        }
        else if(npElem == 9)
        {
            nGP     = 9;
            elShape = ELEM_SHAPE_QUAD;
            degree  = 2;
            getGaussPointsQuad(nGP, gpts1, gpts2, gwts);
        }


        elemVolGP.resize(nGP);

        Nv.resize(nGP);
        dNvdx.resize(nGP);
        dNvdy.resize(nGP);

        elemVol=0.0;
        for(gp=0; gp<nGP; gp++)
        {
            Nv[gp].resize(npElem);
            dNvdx[gp].resize(npElem);
            dNvdy[gp].resize(npElem);

            param[0] = gpts1[gp];
            param[1] = gpts2[gp];

            computeBasisFunctions2D(false, elShape, degree, param, xNode, yNode, &Nv[gp][0], &dNvdx[gp][0], &dNvdy[gp][0], Jac);

            if(Jac < 0.0)
            {
              cout << " Jac = " << Jac << endl;
              cout << " Negative Jacobian in 'femINSstabilised::calcStiffnessAndResidualLinearised2D' " << endl;
              exit(1);
            }

            dvol = gwts[gp]*Jac;
            elemVolGP[gp] = dvol;

            elemVol += dvol;
        }//gp

        //cout << "volume = " << elemVol << endl;
        charlen = sqrt(4.0*abs(elemVol)/PI);

        fact = 4.0/charlen/charlen;
        matG(0,0) = fact;
        matG(1,1) = fact;

        ind = npElem*ndof;

        VectorXd  Flocal(ind);
        MatrixXd  Klocal(ind, ind);

        Klocal.setZero();
        Flocal.setZero();

        force.setZero();
        if(bodyForceTimeFunction > -1)
        {
          force[0] = bodyForce[0]*timeFunction[bodyForceTimeFunction]->getFactor();
          force[1] = bodyForce[1]*timeFunction[bodyForceTimeFunction]->getFactor();
        }

        for(gp=0; gp<nGP; gp++)
        {
            // compute the gradient of velocity
            xx = 0.0; yy = 0.0;
            vel[0] = vel[1] = 0.0;
            velPrev[0] = velPrev[1] = 0.0;
            velDot[0] = velDot[1] = 0.0;
            grad.setZero();
            gradPrev.setZero();
            Du.setZero();
            pres = 0.0;
            dp.setZero();
            phi=0.0; eta=0.0;
            dphi[0] = 0.0; dphi[1] = 0.0;

            for(ii=0; ii<npElem; ii++)
            {
                xx += xNode[ii]*Nv[gp][ii];
                yy += yNode[ii]*Nv[gp][ii];

                TI   = nodeNums[ii]*ndof;
                TIp1 = TI+1;
                TIp2 = TI+2;

                b1 = solnPrev[TI];
                b2 = solnPrev[TIp1];

                velPrev[0] += b1*Nv[gp][ii];
                velPrev[1] += b2*Nv[gp][ii];

                gradPrev(0,0) += b1*dNvdx[gp][ii];
                gradPrev(0,1) += b1*dNvdy[gp][ii];
                gradPrev(1,0) += b2*dNvdx[gp][ii];
                gradPrev(1,1) += b2*dNvdy[gp][ii];

                b1 = solnCur[TI];
                b2 = solnCur[TIp1];

                vel[0]     += b1*Nv[gp][ii];
                vel[1]     += b2*Nv[gp][ii];

                grad(0,0)  += b1*dNvdx[gp][ii];
                grad(0,1)  += b1*dNvdy[gp][ii];
                grad(1,0)  += b2*dNvdx[gp][ii];
                grad(1,1)  += b2*dNvdy[gp][ii];

                b4 = solnCur[TIp2];

                pres       += b4*Nv[gp][ii];
                dp[0]      += b4*dNvdx[gp][ii];
                dp[1]      += b4*dNvdy[gp][ii];

                velDot[0] += solnDotCur[TI]  *Nv[gp][ii];
                velDot[1] += solnDotCur[TIp1]*Nv[gp][ii];
            }


            if(MULTIPHASEFLOW)
            {
              rho = 0.0;
              mu  = 0.0;
              phi = 0.0;
              eta = 0.0;
            dphi[0] = 0.0; dphi[1] = 0.0;

              for(ii=0; ii<npElem; ii++)
              {
                jj = nodeNums[ii];
                b1 = Nv[gp][ii];

                rho += rhoNodal[jj]*b1;
                mu  += muNodal[jj]*b1;
                eta += etaNodal[jj]*b1;

                phi += phiNodal[jj]*b1;

                dphi[0]  += phiNodal[jj]*dNvdx[gp][ii];
                dphi[1]  += phiNodal[jj]*dNvdy[gp][ii];
              }
            }

            // this is pseudo-stress
            stress = mu*grad;
            stress(0,0) -= pres;
            stress(1,1) -= pres;


            gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

            res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
            res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;

            rStab(0) = res2(0) - mu*Du(0) + dp(0) - kappa*eta*dphi[0];
            rStab(1) = res2(1) - mu*Du(1) + dp(1) - kappa*eta*dphi[1] ;

            dvol = elemVolGP[gp];
            if(axsy)
            {
                rad = xx;

                urdr  = vel(0)/rad;
                urdr2 = urdr/rad;
                dvol *= (2.0*PI*rad);

                rStab(0) -= mu*(grad(0,0)/rad - urdr2 );
                rStab(1) -= mu*(grad(1,0)/rad );
            }

            velTemp(0) = velPrev(0);
            velTemp(1) = velPrev(1);
            velTemp(2) = 0.0;

            //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, myTime.dt,  beta, tau);
            //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, myTime.dt,  beta, tau);
            evaluateStabParams_algo3(velTemp, matG, myTime.dt, rho, mu, CI, tau);

            tau[0] *= stabilisationFactors[0];                             // SUPG
            tau[1] *= stabilisationFactors[1];                             // PSPG
            tau[2] *= stabilisationFactors[2];                             // LSIC

            for(ii=0;ii<npElem;ii++)
            {
                TI   = ndof*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
    
                b1 = dNvdx[gp][ii]*dvol;
                b2 = dNvdy[gp][ii]*dvol;
                b4 = Nv[gp][ii]*dvol;
    
                b5 = mu*af*b1;
                b6 = mu*af*b2;
                b8 = af*b4;

                Da = rho*(velPrev(0)*b1 + velPrev(1)*b2)*tau[0];

                for(jj=0;jj<npElem;jj++)
                {
                    TJ   = ndof*jj;
                    TJp1 = TJ+1;
                    TJp2 = TJ+2;
    
                    fact2 = rho*acceFact*Nv[gp][jj];
    
                    // time acceleration term
                    fact = b4*fact2 ;
    
                    // diffusion term
                    fact += ( b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj] );
    
                    Klocal(TI,   TJ)   += fact;
                    Klocal(TIp1, TJp1) += fact;

                    // convection term

                    gradN = gradPrev*(rho*Nv[gp][jj]);

                    Db = rho*(velPrev(0)*dNvdx[gp][jj] + velPrev(1)*dNvdy[gp][jj]);

                    gradN(0,0) += Db;
                    gradN(1,1) += Db;

                    Klocal(TI,   TJ)   += (b8*gradN(0,0));
                    Klocal(TI,   TJp1) += (b8*gradN(0,1));
                    Klocal(TIp1, TJ)   += (b8*gradN(1,0));
                    Klocal(TIp1, TJp1) += (b8*gradN(1,1));
    
                    // pressure term
                    Klocal(TI,   TJp2) -= (b1*af*Nv[gp][jj]);
                    Klocal(TIp1, TJp2) -= (b2*af*Nv[gp][jj]);
    
                    // continuity equation
                    Klocal(TIp2, TJ)   += (b8*dNvdx[gp][jj]);
                    Klocal(TIp2, TJp1) += (b8*dNvdy[gp][jj]);
    
                    // SUPG and PSPG stabilisation terms
    
                    gradN *= af;
    
                    Dj(0,0) = gradN(0,0) + fact2;
                    Dj(0,1) = gradN(0,1);
                    Dj(0,2) = af*dNvdx[gp][jj];
    
                    Dj(1,0) = gradN(1,0);
                    Dj(1,1) = gradN(1,1) + fact2;
                    Dj(1,2) = af*dNvdy[gp][jj];
    
                    if(axsy)
                    {
                      Dj(0,0) -= mu*af*(dNvdx[gp][jj]/rad - Nv[gp][jj]/rad/rad);
                      Dj(1,1) -= mu*af*(dNvdx[gp][jj]/rad);
                    }

                    // SUPG
                    Klocal(TI, TJ)     += Da*Dj(0,0);
                    Klocal(TI, TJp1)   += Da*Dj(0,1);
                    Klocal(TI, TJp2)   += Da*Dj(0,2);
    
                    Klocal(TIp1, TJ)   += Da*Dj(1,0);
                    Klocal(TIp1, TJp1) += Da*Dj(1,1);
                    Klocal(TIp1, TJp2) += Da*Dj(1,2);

                    // PSPG stabilisation
                    Klocal(TIp2, TJ)   += (b1*Dj(0,0) + b2*Dj(1,0))*tau[1];
                    Klocal(TIp2, TJp1) += (b1*Dj(0,1) + b2*Dj(1,1))*tau[1];
                    Klocal(TIp2, TJp2) += (b1*Dj(0,2) + b2*Dj(1,2))*tau[1];

                    // LSIC stabilisation

                    //fact = af*rho*tau[2];

                    //Klocal(TI,   TJ)   += (b1*fact*dN_dx[jj]);
                    //Klocal(TI,   TJp1) += (b1*fact*dN_dy[jj]);

                    //Klocal(TIp1, TJ)   += (b2*fact*dN_dx[jj]);
                    //Klocal(TIp1, TJp1) += (b2*fact*dN_dy[jj]);

                    if(axsy)
                    {
                      // diffusion term
                      Klocal(TI, TJ)     += (b4 * (mu/rad/rad) * (af*Nv[gp][jj]) );
                      Klocal(TI, TJp2)   -= (b4 * Nv[gp][jj]/rad);
    
                      // continuity equation
                      Klocal(TIp2, TJ)   += (b4 * af*Nv[gp][jj]/rad);
                    }
                }

                Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b4*kappa*eta*dphi[0]);
                Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b4*kappa*eta*dphi[1]);
                Flocal(TIp2) -= (b4*grad.trace());

                // SUPG stabilisation terms
                Flocal(TI)   -= Da*rStab(0);
                Flocal(TIp1) -= Da*rStab(1);

                // PSPG stabilisation terms
                Flocal(TIp2) -= (tau[1]*(b1*rStab(0)+b2*rStab(1)));

                // LSIC stabilisation terms
                //fact2 = tau[2]*rho*grad.trace();
                //Flocal(TI)   -= b1*fact2;
                //Flocal(TIp1) -= b2*fact2;

                if(axsy)
                {
                    Flocal(TI)   -= (b4 * (mu/rad/rad) * vel(0) );
                    Flocal(TI)   += (b4 * pres/rad);
                    Flocal(TIp2) -= (b4 * vel(0)/rad);
                }
            }
        } //gp

        //printMatrix(Klocal); printf("\n\n\n");
        //printVector(Flocal); printf("\n\n\n");

        // add contributions from boundary conditions

        size1 = forAssyVec.size();

        //PetscPrintf(MPI_COMM_WORLD, "\n Applying boundary conditions %d \n ", ee);
        if(firstIteration)
        {
            // apply boundary conditions
            for(ii=0; ii<size1; ii++)
            {
                aa = forAssyVec[ii];

                if(aa == -1) // this DOF has a prescribed value
                {
                    fact = solnApplied[globalDOFnums[ii]];

                    // check if fact is zero. We don't need to
                    // execute the for loop if fact is zero.
                    if( abs(fact) > 1.0e-10)
                    {
                        for(jj=0; jj<size1; jj++)
                        {
                            if( forAssyVec[jj] != -1 )
                            {
                                Flocal(jj) -= Klocal(jj, ii) * fact;
                            }
                        }
                    }
                }
            }
        } // if(firstIteration)

        // assemble matrices and vectors
        //PetscPrintf(MPI_COMM_WORLD, "\n Assembling matrices and vectors \n");
        solverPetsc->assembleMatrixAndVectorSerial(forAssyVec, Klocal, Flocal);

        VecSetValues(solverPetsc->reacVec, size1, &globalDOFnums[0], &Flocal[0], ADD_VALUES);

        // add to reaction vector
        //for(ii=0; ii<size1; ii++)
          //reacVec[globalDOFnums[ii]] += Flocal[ii];

    } // if(elem_proc_id[ee] == this_proc_id)
    } //Element Loop


    return 0;
}





int femINSstabilised::calcStiffnessAndResidualLinearised3D(int iter)
{
    int  aa, ee, ii, jj, row, col, ind, size1, size2, npElem;
    int  gp, nGP, TI, TIp1, TIp2, TIp3, TJ, TJp1, TJp2, TJp3, degree;
    ElementShape  elShape;

    double  dvol, b1, b2, b3, b4, b5, b6, b7, b8, xx, yy, zz, Jac;
    double  pres, Da, Db, rad, urdr, urdr2, tau[3], CI=10.0;
    double  fact, fact1, fact2, elemVol, charlen, param[3];
    double  xNode[50], yNode[50], zNode[50];

    VectorXd  res(4), res2(3), dp(3), Du(3), vel(3), velPrev(3), velDot(3), force(3), gradTvel(3), rStab(4);
    MatrixXd  Dj(3,4), grad(3,3), gradN(3,3), stress(3,3), gradPrev(3,3);
    Dj.setZero();
    VectorXd  velTemp(3);


    double  rho = fluidProperties[0];
    double  mu  = fluidProperties[1];
    double  am  = td[1];
    double  af  = td[2];
    double  acceFact = td[8];
    double  muTaf = mu*af;

    bool axsy = false;

    vector<int> nodeNums, forAssyVec, globalDOFnums;
    vector<double>  gpts1, gpts2, gpts3, gwts, elemVolGP;

    MatrixXd  matG(3,3);    matG.setZero();
    vector<VectorXd>  Nv, dNvdx, dNvdy, dNvdz;


    // loop over elements and compute matrix and residual
    PetscPrintf(MPI_COMM_WORLD, "\n Element loop \n\n");
    for(ee=0; ee<mesh->getnElemGlobal(); ee++)
    {
    if(elem_proc_id[ee] == this_mpi_proc)
    {
        //cout << " ee = " << ee << endl;

        // calculate stiffness and residual
        //
        //
        nodeNums      = mesh->elemConn[ee];
        forAssyVec    = forAssyVecAll[ee];
        globalDOFnums = globalDOFnumsAll[ee];
        
        //printVector(nodeNums);
        //printVector(globalDOFnums);
        //printVector(forAssyVec);

        npElem = nodeNums.size();

        for(ii=0;ii<npElem;ii++)
        {
          xNode[ii] = mesh->nodeCoordsOrig[node_map_get_old[nodeNums[ii]]][0];
          yNode[ii] = mesh->nodeCoordsOrig[node_map_get_old[nodeNums[ii]]][1];
          zNode[ii] = mesh->nodeCoordsOrig[node_map_get_old[nodeNums[ii]]][2];
        }

        if(npElem == 4)
        {
            nGP     = 1;
            elShape = ELEM_SHAPE_TETRA;
            degree  = 1;
            getGaussPointsTetra(nGP, gpts1, gpts2, gpts3, gwts);
        }
        else if(npElem == 10)
        {
            nGP     = 4;
            elShape = ELEM_SHAPE_TETRA;
            degree  = 2;
            getGaussPointsTetra(nGP, gpts1, gpts2, gpts3, gwts);
        }
        else if(npElem == 6)
        {
            nGP     = 2;
            elShape = ELEM_SHAPE_WEDGE;
            degree  = 1;
            getGaussPointsWedge(nGP, gpts1, gpts2, gpts3, gwts);
        }
        else if(npElem == 18)
        {
            nGP     = 9;
            elShape = ELEM_SHAPE_WEDGE;
            degree  = 2;
            getGaussPointsWedge(nGP, gpts1, gpts2, gpts3, gwts);
        }
        else if(npElem == 8)
        {
            nGP     = 8;
            elShape = ELEM_SHAPE_HEXA;
            degree  = 1;
            getGaussPointsHexa(nGP, gpts1, gpts2, gpts3, gwts);
        }
        else if(npElem == 27)
        {
            nGP     = 27;
            elShape = ELEM_SHAPE_HEXA;
            degree  = 2;
            getGaussPointsHexa(nGP, gpts1, gpts2, gpts3, gwts);
        }


        elemVolGP.resize(nGP);

        Nv.resize(nGP);
        dNvdx.resize(nGP);
        dNvdy.resize(nGP);
        dNvdz.resize(nGP);

        elemVol=0.0;
        for(gp=0; gp<nGP; gp++)
        {
            Nv[gp].resize(npElem);
            dNvdx[gp].resize(npElem);
            dNvdy[gp].resize(npElem);
            dNvdz[gp].resize(npElem);

            param[0] = gpts1[gp];
            param[1] = gpts2[gp];
            param[2] = gpts3[gp];

            computeBasisFunctions3D(false, elShape, degree, param, xNode, yNode, zNode, &Nv[gp][0], &dNvdx[gp][0], &dNvdy[gp][0], &dNvdz[gp][0], Jac);

            if(Jac < 0.0)
            {
              cout << " Jac = " << Jac << endl;
              cout << " Negative Jacobian in 'femINSstabilised::calcStiffnessAndResidualLinearised3D' " << endl;
              exit(1);
            }

            dvol = gwts[gp]*Jac;
            elemVolGP[gp] = dvol;

            elemVol += dvol;
        }//gp

        //cout << "volume = " << elemVol << endl;
        charlen = pow(6.0*abs(elemVol)/PI, 1.0/3.0);

        fact = 4.0/charlen/charlen;
        matG(0,0) = fact;
        matG(1,1) = fact;
        matG(2,2) = fact;

        ind = npElem*ndof;

        VectorXd  Flocal(ind);
        MatrixXd  Klocal(ind, ind);

        Klocal.setZero();
        Flocal.setZero();

        force.setZero();
        if(bodyForceTimeFunction > -1)
        {
          force[0] = bodyForce[0]*timeFunction[bodyForceTimeFunction]->getFactor();
          force[1] = bodyForce[1]*timeFunction[bodyForceTimeFunction]->getFactor();
          force[2] = bodyForce[2]*timeFunction[bodyForceTimeFunction]->getFactor();
        }

        for(gp=0; gp<nGP; gp++)
        {
            // compute the gradient of velocity
            xx = yy = zz = 0.0;
            vel[0] = vel[1] = vel[2] = 0.0;
            velPrev[0] = velPrev[1] = velPrev[2] = 0.0;
            velDot[0] = velDot[1] = velDot[2] = 0.0;
            grad.setZero();
            gradPrev.setZero();
            Du.setZero();
            pres = 0.0;
            dp.setZero();

            for(ii=0; ii<npElem; ii++)
            {
                xx += xNode[ii]*Nv[gp][ii];
                yy += yNode[ii]*Nv[gp][ii];
                zz += zNode[ii]*Nv[gp][ii];

                TI   = nodeNums[ii]*ndof;
                TIp1 = TI+1;
                TIp2 = TI+2;
                TIp3 = TI+3;

                b1 = solnPrev[TI];
                b2 = solnPrev[TIp1];
                b3 = solnPrev[TIp2];

                velPrev[0] += b1*Nv[gp][ii];
                velPrev[1] += b2*Nv[gp][ii];
                velPrev[2] += b3*Nv[gp][ii];

                gradPrev(0,0) += b1*dNvdx[gp][ii];
                gradPrev(0,1) += b1*dNvdy[gp][ii];
                gradPrev(0,2) += b1*dNvdz[gp][ii];

                gradPrev(1,0) += b2*dNvdx[gp][ii];
                gradPrev(1,1) += b2*dNvdy[gp][ii];
                gradPrev(1,2) += b2*dNvdz[gp][ii];

                gradPrev(2,0) += b3*dNvdx[gp][ii];
                gradPrev(2,1) += b3*dNvdy[gp][ii];
                gradPrev(2,2) += b3*dNvdz[gp][ii];

                b1 = solnCur[TI];
                b2 = solnCur[TIp1];
                b3 = solnCur[TIp2];

                vel[0]     += b1*Nv[gp][ii];
                vel[1]     += b2*Nv[gp][ii];
                vel[2]     += b3*Nv[gp][ii];

                grad(0,0)  += b1*dNvdx[gp][ii];
                grad(0,1)  += b1*dNvdy[gp][ii];
                grad(0,2)  += b1*dNvdz[gp][ii];

                grad(1,0)  += b2*dNvdx[gp][ii];
                grad(1,1)  += b2*dNvdy[gp][ii];
                grad(1,2)  += b2*dNvdz[gp][ii];

                grad(2,0)  += b3*dNvdx[gp][ii];
                grad(2,1)  += b3*dNvdy[gp][ii];
                grad(2,2)  += b3*dNvdz[gp][ii];

                b4 = solnCur[TIp3];

                pres      += b4*Nv[gp][ii];
                dp[0]     += b4*dNvdx[gp][ii];
                dp[1]     += b4*dNvdy[gp][ii];
                dp[2]     += b4*dNvdz[gp][ii];

                velDot[0] += solnDotCur[TI]  *Nv[gp][ii];
                velDot[1] += solnDotCur[TIp1]*Nv[gp][ii];
                velDot[2] += solnDotCur[TIp2]*Nv[gp][ii];
            }

            // this is pseudo-stress
            stress = mu*grad;
            stress(0,0) -= pres;
            stress(1,1) -= pres;
            stress(2,2) -= pres;


            gradTvel = gradPrev*vel + grad*velPrev - gradPrev*velPrev;

            res2(0) = rho*(velDot(0) + gradTvel(0) - force(0)) ;
            res2(1) = rho*(velDot(1) + gradTvel(1) - force(1)) ;
            res2(2) = rho*(velDot(2) + gradTvel(2) - force(2)) ;

            rStab(0) = res2(0) - mu*Du(0) + dp(0) ;
            rStab(1) = res2(1) - mu*Du(1) + dp(1) ;
            rStab(2) = res2(2) - mu*Du(2) + dp(2) ;


            //evaluateStabParams_algo1(&velTemp(0), h, rho, mu, myTime.dt,  beta, tau);
            //evaluateStabParams_algo2(&velTemp(0), h, rho, mu, myTime.dt,  beta, tau);
            evaluateStabParams_algo3(velPrev, matG, myTime.dt, rho, mu, CI, tau);

            tau[0] *= stabilisationFactors[0];                             // SUPG
            tau[1] *= stabilisationFactors[1];                             // PSPG
            tau[2] *= stabilisationFactors[2];                             // LSIC


            dvol = elemVolGP[gp];

            for(ii=0;ii<npElem;ii++)
            {
                TI   = 4*ii;
                TIp1 = TI+1;
                TIp2 = TI+2;
                TIp3 = TI+3;

                b1 = dNvdx[gp][ii]*dvol;
                b2 = dNvdy[gp][ii]*dvol;
                b3 = dNvdz[gp][ii]*dvol;
                b4 = Nv[gp][ii]*dvol;

                b5 = muTaf*b1;
                b6 = muTaf*b2;
                b7 = muTaf*b3;
                b8 = af*b4;

                Da = rho*(velPrev(0)*b1 + velPrev(1)*b2 + velPrev(2)*b3)*tau[0];

                for(jj=0;jj<npElem;jj++)
                {
                    TJ   = 4*jj;
                    TJp1 = TJ+1;
                    TJp2 = TJ+2;
                    TJp3 = TJ+3;

                    fact2 = rho*acceFact*Nv[gp][jj];

                    // time acceleration term
                    fact = b4*fact2 ;

                    // diffusion term
                    fact += (b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj]+b7*dNvdz[gp][jj]);

                    Klocal(TI,   TJ)   += fact;
                    Klocal(TIp1, TJp1) += fact;
                    Klocal(TIp2, TJp2) += fact;

                    // convection term - semi-implicit type B

                    gradN = gradPrev*(rho*Nv[gp][jj]);

                    Db = rho*(velPrev(0)*dNvdx[gp][jj] + velPrev(1)*dNvdy[gp][jj] + velPrev(2)*dNvdz[gp][jj]);

                    gradN(0,0) += Db;
                    gradN(1,1) += Db;
                    gradN(2,2) += Db;

                    Klocal(TI,   TJ)   += (b8*gradN(0,0));
                    Klocal(TI,   TJp1) += (b8*gradN(0,1));
                    Klocal(TI,   TJp2) += (b8*gradN(0,2));

                    Klocal(TIp1, TJ)   += (b8*gradN(1,0));
                    Klocal(TIp1, TJp1) += (b8*gradN(1,1));
                    Klocal(TIp1, TJp2) += (b8*gradN(1,2));

                    Klocal(TIp2, TJ)   += (b8*gradN(2,0));
                    Klocal(TIp2, TJp1) += (b8*gradN(2,1));
                    Klocal(TIp2, TJp2) += (b8*gradN(2,2));

                    // pressure term
                    Klocal(TI,   TJp3) -= (b1*af*Nv[gp][jj]);
                    Klocal(TIp1, TJp3) -= (b2*af*Nv[gp][jj]);
                    Klocal(TIp2, TJp3) -= (b3*af*Nv[gp][jj]);

                    // continuity equation
                    Klocal(TIp3, TJ)   -= (b8*dNvdx[gp][jj]);
                    Klocal(TIp3, TJp1) -= (b8*dNvdy[gp][jj]);
                    Klocal(TIp3, TJp2) -= (b8*dNvdz[gp][jj]);

                    // SUPG and PSPG stabilisation terms
                    //fact2 -= mu*d2N(jj);

                    gradN *= af;

                    Dj(0,0) = gradN(0,0) + fact2;
                    Dj(0,1) = gradN(0,1);
                    Dj(0,2) = gradN(0,2);
                    Dj(0,3) = af*dNvdx[gp][jj];

                    Dj(1,0) = gradN(1,0);
                    Dj(1,1) = gradN(1,1) + fact2;
                    Dj(1,2) = gradN(1,2);
                    Dj(1,3) = af*dNvdy[gp][jj];

                    Dj(2,0) = gradN(2,0);
                    Dj(2,1) = gradN(2,1);
                    Dj(2,2) = gradN(2,2) + fact2;
                    Dj(2,3) = af*dNvdz[gp][jj];

                    // SUPG
                    Klocal(TI, TJ)     += Da*Dj(0,0);
                    Klocal(TI, TJp1)   += Da*Dj(0,1);
                    Klocal(TI, TJp2)   += Da*Dj(0,2);
                    Klocal(TI, TJp3)   += Da*Dj(0,3);

                    Klocal(TIp1, TJ)   += Da*Dj(1,0);
                    Klocal(TIp1, TJp1) += Da*Dj(1,1);
                    Klocal(TIp1, TJp2) += Da*Dj(1,2);
                    Klocal(TIp1, TJp3) += Da*Dj(1,3);

                    Klocal(TIp2, TJ)   += Da*Dj(2,0);
                    Klocal(TIp2, TJp1) += Da*Dj(2,1);
                    Klocal(TIp2, TJp2) += Da*Dj(2,2);
                    Klocal(TIp2, TJp3) += Da*Dj(2,3);

                    // PSPG
                    Klocal(TIp3, TJ)   -= (b1*Dj(0,0) + b2*Dj(1,0) + b3*Dj(2,0))*tau[1];
                    Klocal(TIp3, TJp1) -= (b1*Dj(0,1) + b2*Dj(1,1) + b3*Dj(2,1))*tau[1];
                    Klocal(TIp3, TJp2) -= (b1*Dj(0,2) + b2*Dj(1,2) + b3*Dj(2,2))*tau[1];
                    Klocal(TIp3, TJp3) -= (b1*Dj(0,3) + b2*Dj(1,3) + b3*Dj(2,3))*tau[1];

                    // LSIC stabilisation

                    fact2 = rho*af*tau[2];

                    //Klocal(TI,   TJ)   += (b1*dNvdx[gp][jj])*fact2;
                    //Klocal(TI,   TJp1) += (b1*dNvdy[gp][jj])*fact2;
                    //Klocal(TI,   TJp2) += (b1*dNvdz[gp][jj])*fact2;

                    //Klocal(TIp1, TJ)   += (b2*dNvdx[gp][jj])*fact2;
                    //Klocal(TIp1, TJp1) += (b2*dNvdy[gp][jj])*fact2;
                    //Klocal(TIp1, TJp2) += (b2*dNvdz[gp][jj])*fact2;

                    //Klocal(TIp2, TJ)   += (b3*dNvdx[gp][jj])*fact2;
                    //Klocal(TIp2, TJp1) += (b3*dNvdy[gp][jj])*fact2;
                    //Klocal(TIp2, TJp2) += (b3*dNvdz[gp][jj])*fact2;
                }

                Flocal(TI)   -= (b4*res2(0) + b1*stress(0,0) + b2*stress(0,1) + b3*stress(0,2) );
                Flocal(TIp1) -= (b4*res2(1) + b1*stress(1,0) + b2*stress(1,1) + b3*stress(1,2) );
                Flocal(TIp2) -= (b4*res2(2) + b1*stress(2,0) + b2*stress(2,1) + b3*stress(2,2) );
                Flocal(TIp3) += (b4*grad.trace());

                // SUPG stabilisation terms
                Flocal(TI)   -= Da*rStab(0);
                Flocal(TIp1) -= Da*rStab(1);
                Flocal(TIp2) -= Da*rStab(2);

                // PSPG stabilisation terms
                Flocal(TIp3) += (tau[1]*(b1*rStab(0)+b2*rStab(1)+b3*rStab(2)));

                // LSIC stabilisation terms
                //fact2 = tau[2]*rho*grad.trace();

                //Flocal(TI)   -= b1*fact2;
                //Flocal(TIp1) -= b2*fact2;
                //Flocal(TIp2) -= b3*fact2;
            }
        } //gp

        //printMatrix(Klocal); printf("\n\n\n");
        //printVector(Flocal); printf("\n\n\n");

        // add contributions from boundary conditions

        size1 = forAssyVec.size();

        //PetscPrintf(MPI_COMM_WORLD, "\n Applying boundary conditions %d \n ", ee);
        if(firstIteration)
        {
            // apply boundary conditions
            for(ii=0; ii<size1; ii++)
            {
                aa = forAssyVec[ii];

                if(aa == -1) // this DOF has a prescribed value
                {
                    fact = solnApplied[globalDOFnums[ii]];

                    // check if fact is zero. We don't need to
                    // execute the for loop if fact is zero.
                    if( abs(fact) > 1.0e-10)
                    {
                        for(jj=0; jj<size1; jj++)
                        {
                            if( forAssyVec[jj] != -1 )
                            {
                                Flocal(jj) -= Klocal(jj, ii) * fact;
                            }
                        }
                    }
                }
            }
        } // if(firstIteration)

        // assemble matrices and vectors
        //PetscPrintf(MPI_COMM_WORLD, "\n Assembling matrices and vectors \n");
        solverPetsc->assembleMatrixAndVectorSerial(forAssyVec, Klocal, Flocal);

        VecSetValues(solverPetsc->reacVec, size1, &globalDOFnums[0], &Flocal[0], ADD_VALUES);

    } // if(elem_proc_id[ee] == this_proc_id)
    } //Element Loop


    return 0;
}









int femINSstabilised::addExternalForces(double loadFact)
{
    int  nn, dof, ii, ind;
    double specVal=0.0;

    VectorXd  vecTemp, Flocal;
    vecTemp.setZero();
/*
    // specified nodal forces
    for(ii=0;ii<nodeForcesData.size();++ii)
    {
      nn  = (int) (nodeForcesData[ii][0] - 1);
      dof = (int) (nodeForcesData[ii][1] - 1);
      specVal = nodeForcesData[ii][2];

      ind = nn*ndof+dof;

      vecTemp[ind] += specVal;
    }
    //printVector(vecTemp);
*/
    return 0;
}




int femINSstabilised::computeElementErrors(int ind)
{
    cout << " Computing errors \n " << endl;

    double totalError = 0.0, timeNow = 5.0;
    //cout << " index = " << index << endl;
    for(int index=0; index<4; index++)
    {
      totalError = 0.0;
      for(int ee=0; ee<mesh->getnElemGlobal(); ee++)
      {
        //Compute the element force vector, including residual force
        //totalError += elems[ee]->CalculateError(nodeCoordsOrig, elemData, timeData, soln, solnDot, timeNow, index);
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




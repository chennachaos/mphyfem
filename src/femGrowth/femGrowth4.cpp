
#include "SolverPardisoEigen.h"
#include "femGrowth.h"
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


extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime myTime;
extern bool debug;

//using namespace Spectra;





/*
int femGrowth::solveStepDeflation(int iter_max, double tol1)
{
    cout << "femGrowth::solveStepDeflation " << endl;
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

            COUT << "femGrowth"; printf("    %d \t %11.4e\n", (iter+1), rNorm);

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



int femGrowth::solveStepDeflation()
{
    cout << "femGrowth::solveStepDeflation " << endl;
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

            printf(" femGrowth   %d \t %11.4e\n", (iter+1), rNorm);

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
int femGrowth::solveWithNewtonRaphson()
{
    cout << "femGrowth::solveWithNewtonRaphson " << endl;

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

        printf(" femGrowth ...  %d \t %11.4e \n", (iter+1), rNorm);
        //printf(" femGrowth ...  %d \t %11.4e \t %11.4e \n", (iter+1), rNormPrev, rNorm);

        // check for convergence and divergence of the iterations
        if( converged() )
        {
          printf(" femGrowth ...  Iterations CONVERGED \n");
          break;
        }
        else if( (rNorm/rNormPrev) > 10000.0 )
        {
          printf(" femGrowth ...  Iterations are diverging. NR loop is terminated. \n");
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



int femGrowth::solveStepArcLengthGrowthModel()
{
    cout << "femGrowth::solveStepArcLengthGrowthModel " << endl;
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

        printf("  femGrowth  %d \t %11.4e\n", (iter+1), rNorm);

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





void  femGrowth::computerInternalForceDerivativeOfGrowth()
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




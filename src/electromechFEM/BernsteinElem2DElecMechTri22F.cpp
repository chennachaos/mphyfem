
#include "BernsteinElem2DElecMechTri22F.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "utilitiesmaterial.h"
#include "BasisFunctionsBernstein.h"

using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



BernsteinElem2DElecMechTri22F::BernsteinElem2DElecMechTri22F()
{
  ndim   = 2;
  degree = 2;
  npElem = 6;
  nlbfU  = 6;
  nlbfF  = 6;
  ndof   = 2;
  nsize  = npElem*ndof;

  Nc.resize(nlbfU);
  dNc_dx.resize(nlbfU);
  dNc_dy.resize(nlbfU);

  pres = presPrev = 0.0;
  AlgoType = 1;
}


BernsteinElem2DElecMechTri22F::~BernsteinElem2DElecMechTri22F()
{
}


void BernsteinElem2DElecMechTri22F::prepareElemData()
{
  ElementBase::prepareElemData();

  //AlgoType = (int) elmDat[9];
  AlgoType = 1;

  return;
}



double BernsteinElem2DElecMechTri22F::computeVolume(bool init)
{
  double  dvol, Jac, param[2];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

  int   ii, gp;

  double xNode[6], yNode[6], xx, yy;

  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
  }

  vector<double>  gausspoints1, gausspoints2, gaussweights;

  getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    elemVol=0.0;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          if(init)
            GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          else
            GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);

          xx = yy= 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
          }

          if(axsy)
            dvol *= 2.0*PI*xx;

          elemVol += dvol;
  }//gp

  //printMatrix(Kuu);  printf("\n\n\n");  printVector(FlocalU);

  return elemVol;
}



int BernsteinElem2DElecMechTri22F::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    int  ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  fact, dvol0, Jac, bb1, cc1, param[2];
    double rho0 = elmDat[5] ;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Mlocal.rows() != nsize)
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int  nGPt=7;
    getGaussPointsTriangle(nGPt, gausspoints1, gausspoints2, gaussweights);

    elemVol=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(thick*Jac);

        xx = yy= 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        if(axsy)
          dvol0 *= 2.0*PI*yy;

        elemVol += dvol0;

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = (dvol0*rho0)*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          for(jj=0; jj<nlbfU; jj++)
          {
              TJ   = 2*jj;
              TJp1 = TJ+1;

              fact  = bb1*N[jj];

              Mlocal(TI,   TJ)    += fact ;
              Mlocal(TIp1, TJp1)  += fact ;
          }
        }
    } //gp
    //cout << " elemVol = " << elemVol << endl;
    //printMatrix(Kuu);  printf("\n\n\n");

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      for(ii=0; ii<nsize; ii++)
      {
        fact=0.0;
        for(jj=0; jj<nsize; jj++)
        {
          fact += Mlocal(ii,jj);
          Mlocal(ii,jj) = 0.0;
        }
        Mlocal(ii,ii) = fact;
      }
    }

    return 0;
}



double  BernsteinElem2DElecMechTri22F::calcCriticalTimeStep(bool flag)
{
  double  dtCric;

  return  dtCric;
}




int BernsteinElem2DElecMechTri22F::calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter)
{
  if(finite)
  {
    if(AlgoType == 1)
      calcStiffnessAndResidualFS1(Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, firstIter);
    else
      calcStiffnessAndResidualFS2(Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, firstIter);
  }
  else
    calcStiffnessAndResidualSS(Kuu, Kuf, Kfu, Kff, FlocalU, FlocalF, firstIter);

  return 0;
}



int BernsteinElem2DElecMechTri22F::calcResidual(VectorXd& FlocalU)
{
  if(finite)
  {
    if(AlgoType == 1)
      calcResidualFS1(FlocalU);
    else
      calcResidualFS2(FlocalU);
  }
  else
    calcResidualSS(FlocalU);

  return 0;
}


int BernsteinElem2DElecMechTri22F::calcStiffnessAndResidualSS(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, F33, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    double rho  = elmDat[5] ;
    double rho0 = rho ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

    VectorXd  Nc(nlbfU), dNc_dx(nlbfU), dNc_dy(nlbfU);
    double  Fc[4], trFc, detFc, Bbar[4][2];
    double  fact1, fact2;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    vector<double>  gausspoints1, gausspoints2, gaussweights;

    if(AlgoType == 1)
    {
      //quadrature point for the element centroid
      getGaussPointsTriangle(1, gausspoints1, gausspoints2, gaussweights);

      param[0] = gausspoints1[0];
      param[1] = gausspoints2[0];

      GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nc(0), &dNc_dx(0), &dNc_dy(0), Jac);

      //quadrature points for the element integration
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);
    }
    else
    {
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

      elemVol=0.0;
      dNc_dx.setZero();
      dNc_dy.setZero();
      for(gp=0; gp<nGP; gp++)
      {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        elemVol  += dvol0;

        dNc_dx += dN_dx*dvol0;
        dNc_dy += dN_dy*dvol0;
      }
      dNc_dx /= elemVol;
      dNc_dy /= elemVol;
    }

    Fc[0] = computeValueCur(0, dNc_dx) + 1.0;
    Fc[2] = computeValueCur(0, dNc_dy);
    Fc[1] = computeValueCur(1, dNc_dx);
    Fc[3] = computeValueCur(1, dNc_dy) + 1.0;

    detFc = Fc[0]*Fc[3] - Fc[1]*Fc[2];
    trFc  = Fc[0] + Fc[3];


    // resize local matrices and initialise them to zero
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)     Kuu.resize(nsize, nsize);
    Kuu.setZero();
    if(Kuf.rows() < nsize || Kuf.cols() < nlbfF)     Kuf.resize(nsize, nlbfF);
    Kuf.setZero();
    if(Kfu.rows() < nlbfF || Kfu.cols() < nsize)     Kfu.resize(nlbfF, nsize);
    Kfu.setZero();
    if(Kff.rows() < nlbfF || Kff.cols() < nlbfF)     Kff.resize(nlbfF, nlbfF);
    Kff.setZero();

    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalF.rows() < nlbfF)   FlocalF.resize(nlbfF);
    FlocalF.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        veloCur[0] = computeValueDotCur(0, N);
        veloCur[1] = computeValueDotCur(1, N);

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;

        fact = (Fc[0] + Fc[3] - F[0] - F[3])/3.0;

        F[0] += fact ;
        F[3] += fact ;
        F33 = 1.0 + fact ;
        //F33 = 1.0;


        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol;
        }

        // Calculate Stiffness and Residual
        //==============================================

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dN_dx[ii];
          bb2 = dN_dy[ii];
          bb3 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          fact1 = (dNc_dx(ii) - bb1)/3.0;
          fact2 = (dNc_dy(ii) - bb2)/3.0;

          Bbar[0][0] = bb1+fact1;
          Bbar[1][0] = fact1;
          Bbar[2][0] = fact1;
          //Bbar[2][0] = 0.0;
          Bbar[3][0] = bb2;

          Bbar[0][1] = fact2;
          Bbar[1][1] = bb2+fact2;
          Bbar[2][1] = fact2;
          //Bbar[2][1] = 0.0;
          Bbar[3][1] = bb1;

          bc[0][0] = (Bbar[0][0] * cc[0][0] + Bbar[1][0] * cc[1][0] + Bbar[2][0] * cc[2][0] + Bbar[3][0] * cc[3][0]);
          bc[0][1] = (Bbar[0][0] * cc[0][1] + Bbar[1][0] * cc[1][1] + Bbar[2][0] * cc[2][1] + Bbar[3][0] * cc[3][1]);
          bc[0][2] = (Bbar[0][0] * cc[0][2] + Bbar[1][0] * cc[1][2] + Bbar[2][0] * cc[2][2] + Bbar[3][0] * cc[3][2]);
          bc[0][3] = (Bbar[0][0] * cc[0][3] + Bbar[1][0] * cc[1][3] + Bbar[2][0] * cc[2][3] + Bbar[3][0] * cc[3][3]);

          bc[1][0] = (Bbar[0][1] * cc[0][0] + Bbar[1][1] * cc[1][0] + Bbar[2][1] * cc[2][0] + Bbar[3][1] * cc[3][0]);
          bc[1][1] = (Bbar[0][1] * cc[0][1] + Bbar[1][1] * cc[1][1] + Bbar[2][1] * cc[2][1] + Bbar[3][1] * cc[3][1]);
          bc[1][2] = (Bbar[0][1] * cc[0][2] + Bbar[1][1] * cc[1][2] + Bbar[2][1] * cc[2][2] + Bbar[3][1] * cc[3][2]);
          bc[1][3] = (Bbar[0][1] * cc[0][3] + Bbar[1][1] * cc[1][3] + Bbar[2][1] * cc[2][3] + Bbar[3][1] * cc[3][3]);

          TI   = 2*ii;
          TIp1 = TI+1;

          FlocalU(TI)   += bb5*force[0];
          FlocalU(TIp1) += bb5*force[1];

          FlocalU(TI)   -= (Bbar[0][0]*stre[0] + Bbar[1][0]*stre[1] + Bbar[2][0]*stre[2] + Bbar[3][0]*stre[3]) ;
          FlocalU(TIp1) -= (Bbar[0][1]*stre[0] + Bbar[1][1]*stre[1] + Bbar[2][1]*stre[2] + Bbar[3][1]*stre[3]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = N[jj];

            fact1 = (dNc_dx(jj) - cc1)/3.0;
            fact2 = (dNc_dy(jj) - cc2)/3.0;

            Bbar[0][0] = cc1+fact1;
            Bbar[1][0] = fact1;
            Bbar[2][0] = fact1;
            //Bbar[2][0] = 0.0;
            Bbar[3][0] = cc2;

            Bbar[0][1] = fact2;
            Bbar[1][1] = cc2+fact2;
            Bbar[2][1] = fact2;
            //Bbar[2][1] = 0.0;
            Bbar[3][1] = cc1;

            TJ   = 2*jj;
            TJp1 = TJ+1;

            // acceleration term
            acceFact2 = acceFact1*cc3*rho0;

            fact  = bb5*acceFact2;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;

            Kuu(TI,TJ)     += af*(bc[0][0] * Bbar[0][0] + bc[0][1] * Bbar[1][0] + bc[0][2] * Bbar[2][0] + bc[0][3] * Bbar[3][0]) ;
            Kuu(TI,TJp1)   += af*(bc[0][0] * Bbar[0][1] + bc[0][1] * Bbar[1][1] + bc[0][2] * Bbar[2][1] + bc[0][3] * Bbar[3][1]) ;
            Kuu(TIp1,TJ)   += af*(bc[1][0] * Bbar[0][0] + bc[1][1] * Bbar[1][0] + bc[1][2] * Bbar[2][0] + bc[1][3] * Bbar[3][0]) ;
            Kuu(TIp1,TJp1) += af*(bc[1][0] * Bbar[0][1] + bc[1][1] * Bbar[1][1] + bc[1][2] * Bbar[2][1] + bc[1][3] * Bbar[3][1]) ;
          }
        }
    }//gp

    return 0;
}




int BernsteinElem2DElecMechTri22F::calcResidualSS(VectorXd& FlocalU)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, F33, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    double rho  = elmDat[5] ;
    double rho0 = rho ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;

    VectorXd  Nc(nlbfU), dNc_dx(nlbfU), dNc_dy(nlbfU);
    double  Fc[4], detFc, trFc, Bbar[4][2], fact1, fact2;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    vector<double>  gausspoints1, gausspoints2, gaussweights;

    if(AlgoType == 1)
    {
      getGaussPointsTriangle(1, gausspoints1, gausspoints2, gaussweights);

      param[0] = gausspoints1[0];
      param[1] = gausspoints2[0];

      GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nc(0), &dNc_dx(0), &dNc_dy(0), Jac);
    }
    else
    {
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

      elemVol=0.0;
      dNc_dx.setZero();
      dNc_dy.setZero();
      for(gp=0; gp<nGP; gp++)
      {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        elemVol  += dvol0;

        dNc_dx += dN_dx*dvol0;
        dNc_dy += dN_dy*dvol0;
      }
      dNc_dx /= elemVol;
      dNc_dy /= elemVol;
    }

    //Fc[0] = computeValueCur(0, dNc_dx) + 1.0;
    //Fc[2] = computeValueCur(0, dNc_dy);
    //Fc[1] = computeValueCur(1, dNc_dx);
    //Fc[3] = computeValueCur(1, dNc_dy) + 1.0;

    Fc[0] = computeValue(0, dNc_dx) + 1.0;
    Fc[2] = computeValue(0, dNc_dy);
    Fc[1] = computeValue(1, dNc_dx);
    Fc[3] = computeValue(1, dNc_dy) + 1.0;

    detFc = Fc[0]*Fc[3] - Fc[1]*Fc[2];
    trFc  = Fc[0] + Fc[3];

    if(AlgoType == 1)
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);


    if(FlocalU.rows() != nsize)
      FlocalU.resize(nsize);
    FlocalU.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
          F33 = 1.0;

        fact = (Fc[0] + Fc[3] - F[0] - F[3])/3.0;

        F[0] += fact ;
        F[3] += fact ;
        F33 = 1.0 + fact ;


        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        // Calculate Residual
        //==============================================

        // body force terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        for(ii=0;ii<4;ii++)
          stre[ii] *= dvol;


        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dN_dx[ii];
          bb2 = dN_dy[ii];
          bb3 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          fact1 = (dNc_dx(ii) - bb1)/3.0;
          fact2 = (dNc_dy(ii) - bb2)/3.0;

          Bbar[0][0] = bb1+fact1;
          Bbar[1][0] = fact1;
          Bbar[2][0] = fact1;
          Bbar[3][0] = bb2;

          Bbar[0][1] = fact2;
          Bbar[1][1] = bb2+fact2;
          Bbar[2][1] = fact2;
          Bbar[3][1] = bb1;

          TI   = 2*ii;
          TIp1 = TI+1;

          FlocalU(TI)   += bb5*force[0];
          FlocalU(TIp1) += bb5*force[1];

          FlocalU(TI)   -= (Bbar[0][0]*stre[0] + Bbar[1][0]*stre[1] + Bbar[2][0]*stre[2] + Bbar[3][0]*stre[3]) ;
          FlocalU(TIp1) -= (Bbar[0][1]*stre[0] + Bbar[1][1]*stre[1] + Bbar[2][1]*stre[2] + Bbar[3][1]*stre[3]) ;
        }
    }//gp

    return 0;
}



int BernsteinElem2DElecMechTri22F::calcStiffnessAndResidualFS1(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter)
{
    // F-bar formulation by Eduardo deSauza Neto and Peric and co.
    //cout << "  F-bar formulation by Eduardo deSauza Neto and Peric and co. " << endl;

    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), dNu_dz(nlbfU);
    VectorXd  Nf(nlbfU), dNf_dx(nlbfU), dNf_dy(nlbfU), dNf_dz(nlbfU);
    VectorXd  Nuc(nlbfU), dNuc_dx(nlbfU), dNuc_dy(nlbfU), dNuc_dz(nlbfU);

    double  detF, detFc, trF, trFc, fact, fact1, dvol, dvol0, Jac, volstrain, phi, r2d3 = 2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[3], bforce[3], force[3], tarr[6],  tarr2[6];
    double  veloCur[3], acceCur[3], sig[3];
    MatrixXd  F(3,3), Fn(3,3), Fc(3,3), Cmat(9,9), Bmat(3,9), Amat(3,3),  Gc(3,9);
    Amat.setZero();
    Bmat.setZero();
    Cmat.setZero();
    Gc.setZero();
    F.setZero();
    Fn.setZero();
    Fc.setZero();

    VectorXd  stre(9),  elecDisp(3),  elecField(3), qVec(9);
    stre.setZero();

    double BULK = matDat[0];
    double mu   = matDat[1];
    double LAMBDA = BULK-r2d3*mu;
    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    //  compute deformation gradient at the element center
    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(1, gausspoints1, gausspoints2, gaussweights);

    param[0] = gausspoints1[gp];
    param[1] = gausspoints2[gp];

    // basis functions wrt to the original configuration
    GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nuc(0), &dNuc_dx(0), &dNuc_dy(0), Jac);

    Fc(0,0) = computeValueCur(0, dNuc_dx) + 1.0;            // eps_xx
    Fc(0,1) = computeValueCur(1, dNuc_dx);                  // eps_yx
    Fc(1,0) = computeValueCur(0, dNuc_dy);                  // eps_xy
    Fc(1,1) = computeValueCur(1, dNuc_dy) + 1.0;            // eps_yy
    Fc(2,2) = 1.0;                                          // eps_zz

    trFc  = Fc(0,0) + Fc(1,1);
    detFc = Fc(0,0)*Fc(1,1) - Fc(0,1)*Fc(1,0);

    // basis functions wrt to the current configuration
    GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nuc(0), &dNuc_dx(0), &dNuc_dy(0), Jac);


    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    elemVolOrig=0.0;
    elemVolCur=0.0;

    // resize local matrices and initialise them to zero
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)     Kuu.resize(nsize, nsize);
    Kuu.setZero();
    if(Kuf.rows() < nsize || Kuf.cols() < nlbfF)     Kuf.resize(nsize, nlbfF);
    Kuf.setZero();
    if(Kfu.rows() < nlbfF || Kfu.cols() < nsize)     Kfu.resize(nlbfF, nsize);
    Kfu.setZero();
    if(Kff.rows() < nlbfF || Kff.cols() < nlbfF)     Kff.resize(nlbfF, nlbfF);
    Kff.setZero();

    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalF.rows() < nlbfF)   FlocalF.resize(nlbfF);
    FlocalF.setZero();


    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

        Nf     = Nu;
        dNf_dx = dNu_dx;
        dNf_dy = dNu_dy;

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F(0,0) = computeValueCur(0, dNu_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(1, dNu_dx);                  // eps_yx
        F(1,0) = computeValueCur(0, dNu_dy);                  // eps_xy
        F(1,1) = computeValueCur(1, dNu_dy) + 1.0;            // eps_yy

        volstrain = F(0) + F(4) - 2.0;
        detF = F(0)*F(4) - F(1)*F(3);

        if(detF < 0.0)   return 1;

        acceCur[0] = computeValueDotDotCur(0, Nu);
        acceCur[1] = computeValueDotDotCur(1, Nu);

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);

          dNf_dx = dNu_dx;
          dNf_dy = dNu_dy;
        }
        elemVolCur  += dvol;


        // evaluate phi at the quadrature points
        elecField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          elecField[0] -= (dNf_dx[ii]*SolnData->var3[nodeNums[ii]]);
          elecField[1] -= (dNf_dy[ii]*SolnData->var3[nodeNums[ii]]);
        }

        //  modify F
        if(sss == 1)  // plane stress
        {
          if(finite)
            F(2,2) = 1.0/sqrt(detF);
          else
            F(2,2) = 3.0 - F(0,0) - F(1,1);
        }
        else if(sss == 2)    // plane strain
        {
          fact = sqrt(detFc/detF);

          F *= fact ;

          F(2,2) = 1.0;
        }

        pres = 0.0;
        MatlData->computeStressAndTangent(true, sss, Fn, F, elecField, pres, stre, elecDisp, Cmat, Bmat, Amat, ivar, gp, dt);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        // q tensor for the F-bar formulation
        qVec = -0.5*stre;
        for(ii=0; ii<5; ii++)
        {
            qVec(ii) += ( Cmat(ii,0)+Cmat(ii,4)+stre(ii) )/2.0;
        }

        // Calculate Stiffness and Residual
        //==============================================

        xx = yy = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += Nu[ii]*xNode[ii];
          yy += Nu[ii]*yNode[ii];
        }

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        for(ii=0; ii<nlbfU; ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb4 = dvol*Nu[ii];
          bb5 = dvol0*Nu[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          //  Kuu matrix
          for(jj=0; jj<5; jj++)
          {
            Gc(0,jj) = bb1*Cmat(0,jj) + bb2*Cmat(3,jj) ;
            Gc(1,jj) = bb1*Cmat(1,jj) + bb2*Cmat(4,jj) ;
          }

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[1] + bb2*stre[4];

          FlocalU(TI)   += (bb5*force[0] - sig[0]) ;
          FlocalU(TIp1) += (bb5*force[1] - sig[1]) ;

          tarr[0] = bb1*qVec[0] + bb2*qVec[3] ;
          tarr[1] = bb1*qVec[1] + bb2*qVec[4] ;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dNu_dx[jj];
            cc2 = dNu_dy[jj];
            cc4 = Nu[jj];

            TJ   = 2*jj;
            TJp1 = TJ+1;

            // acceleration term
            acceFact2 = acceFact1*cc4*rho0;

            fact  = bb5*acceFact2;

            // material Stiffness contribution
            fact  += af*(sig[0]*cc1+sig[1]*cc2)*FiniteFact;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;

            Kuu(TI,   TJ)    +=  af*(Gc(0,0)*cc1 + Gc(0,3)*cc2) ;
            Kuu(TI,   TJp1)  +=  af*(Gc(0,1)*cc1 + Gc(0,4)*cc2) ;

            Kuu(TIp1, TJ)    +=  af*(Gc(1,0)*cc1 + Gc(1,3)*cc2) ;
            Kuu(TIp1, TJp1)  +=  af*(Gc(1,1)*cc1 + Gc(1,4)*cc2) ;

            //  additional terms due to F-bar formulation
            cc1 = af*(dNuc_dx(jj) - cc1);
            cc2 = af*(dNuc_dy(jj) - cc2);

            Kuu(TI,   TJ)    +=  tarr[0]*cc1;
            Kuu(TI,   TJp1)  +=  tarr[0]*cc2;

            Kuu(TIp1, TJ)    +=  tarr[1]*cc1;
            Kuu(TIp1, TJp1)  +=  tarr[1]*cc2;
          }
          //printMatrix(Kuu);

          //  Kuf & Kfu matrices
          Gc(0,0) = bb1*Bmat(0,0) + bb2*Bmat(0,3) ;
          Gc(0,1) = bb1*Bmat(1,0) + bb2*Bmat(1,3) ;

          Gc(1,0) = bb1*Bmat(0,1) + bb2*Bmat(0,4) ;
          Gc(1,1) = bb1*Bmat(1,1) + bb2*Bmat(1,4) ;

          for(jj=0; jj<nlbfF; jj++)
          {
            tarr2[0] = Gc(0,0)*dNf_dx[jj] + Gc(0,1)*dNf_dy[jj] ;
            tarr2[1] = Gc(1,0)*dNf_dx[jj] + Gc(1,1)*dNf_dy[jj] ;

            Kuf(TI,   jj)  += tarr2[0];
            Kuf(TIp1, jj)  += tarr2[1];

            Kfu(jj, TI)    += tarr2[0];
            Kfu(jj, TIp1)  += tarr2[1];
          }
        }

        //  Kff matrix
        for(ii=0; ii<nlbfF; ii++)
        {
          bb1 = dvol*dNf_dx[ii];
          bb2 = dvol*dNf_dy[ii];
          bb4 = dvol*Nf[ii];
          bb5 = dvol0*Nf[ii];

          FlocalF(ii) -= (bb1*elecDisp(0)+bb2*elecDisp(1));

          tarr[0] = bb1*Amat(0,0) + bb2*Amat(1,0) ;
          tarr[1] = bb1*Amat(0,1) + bb2*Amat(1,1) ;

          for(jj=0; jj<nlbfF; jj++)
          {
            Kff(ii, jj) -= (dNf_dx[jj]*tarr[0] + dNf_dy[jj]*tarr[1]);
          }
        }
    }//gp

    //if(elenum < 10)
    //{
      //printMatrix(Kuu);
      //printMatrix(Kuf);
      //printMatrix(Kfu);
      //printMatrix(Kff);
    //}

    return 0;
}




int BernsteinElem2DElecMechTri22F::calcResidualFS1(VectorXd& FlocalU)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, F33, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], cc[4][4], bc[2][4], param[2], bforce[2], force[2];

    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;

    double  aa[5][5], cch[5], Fc[4], fact1, fact2;
    double  detFc, trFc;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    //
    vector<double>  gausspoints1, gausspoints2, gaussweights;


    if(AlgoType == 1)
    {
      getGaussPointsTriangle(1, gausspoints1, gausspoints2, gaussweights);

      param[0] = gausspoints1[0];
      param[1] = gausspoints2[0];

      GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nc(0), &dNc_dx(0), &dNc_dy(0), Jac);

      Fc[0] = computeValue(0, dNc_dx) + 1.0;
      Fc[2] = computeValue(0, dNc_dy);
      Fc[1] = computeValue(1, dNc_dx);
      Fc[3] = computeValue(1, dNc_dy) + 1.0;

      detFc = Fc[0]*Fc[3] - Fc[1]*Fc[2];
      trFc  = Fc[0] + Fc[3];

      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);
    }
    else
    {
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

      elemVol=0.0;
      dNc_dx.setZero();
      dNc_dy.setZero();
      detFc = 0.0;
      for(gp=0; gp<nGP; gp++)
      {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        // volume in the original configuration
        dvol0 = gaussweights[gp]*(Jac*thick);

        Fc[0] = computeValue(0, dN_dx) + 1.0;
        Fc[2] = computeValue(0, dN_dy);
        Fc[1] = computeValue(1, dN_dx);
        Fc[3] = computeValue(1, dN_dy) + 1.0;

        // determinant of the deformation gradient
        detF = Fc[0]*Fc[3] - Fc[1]*Fc[2];

        // volume in the deformed configuration
        dvol = dvol0*detF;

        elemVol  += dvol;

        detFc += detF*dvol;
      }
      detFc /= elemVol;
    }


    // resize local matrices and initialise them to zero
    if(FlocalU.rows() != nsize)
      FlocalU.resize(nsize);
    FlocalU.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)          return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);
        }

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          if(finite)
          {
            if(AlgoType == 1)
              fact = sqrt(detFc/detF);
            else
              fact = pow(detFc/detF, 1.0/3.0);

            for(ii=0; ii<4; ii++)
              F[ii] = F[ii] * fact ;

            F33 = 1.0;
          }
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)          return 1;

        if(finite)
          dvol *= F33;


        // Calculate Residual
        //==============================================

        // body force terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          FlocalU(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3]) ;
          FlocalU(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1]) ;
        }
    } //gp

  return 0;
}



/*
int BernsteinElem2DElecMechTri22F::calcStiffnessAndResidualFS2(MatrixXd& Kuu, VectorXd& FlocalU)
{
    //F-bar formulation by Ejgudej

    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;

    double  F33, fact, dvol, dvol0, Jac, r1d3 = 1.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, fact1, fact2;
    double  fact3, fact4, fact5, fact6, fact7, fact8, fact9;
    double  stre[4], cc[4][4], bc[2][6], param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    double  BULK = matDat[0];
    double  eps  = 1.0/BULK;
    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;


    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);


    double  F[nGP][4], detF[nGP], volume[nGP], volume0[nGP];
    double  Fc[4], detFc, trFc, dummy, alpha;
    double  MAB, invMAB;

    double  Bbar[nlbfU][4][2], g[2][nlbfU], gmN[2][nlbfU];
    double  dN_dx[nGP][nlbfU], dN_dy[nGP][nlbfU];

    MatrixXd  N(nGP,nlbfU);//, dN_dx(nGP,nlbfU), dN_dy(nGP,nlbfU);
    VectorXd  NN(nlbfU), dNN_dx(nlbfU), dNN_dy(nlbfU);
    VectorXd  Cmat1(nlbfU), Cmat2(nlbfU), Wmat1(nlbfU), Wmat2(nlbfU);
    double  detFbar, detJbar1, detJbar2, g2[nsize][nsize], g3[nsize][nsize];
    double  Cmat3[nsize][nsize], Wmat3[nsize][nsize], Cmat4[nsize][nsize], Wmat4[nsize][nsize]  ;


    MAB=0.0;
    invMAB=0.0;
    detJbar1 = 0.0;
    detJbar2 = 0.0;
    Cmat1.setZero();
    Cmat2.setZero();
    Wmat1.setZero();
    Wmat2.setZero();

    elemVol=0.0;
    dNc_dx.setZero();
    dNc_dy.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &NN(0), &dNN_dx(0), &dNN_dy(0), Jac);

        // volume in the original configuration
        dvol0 = gaussweights[gp]*(Jac*thick);
        volume0[gp] = dvol0;

        // compute Mass matrix
        MAB += dvol0;

        F[gp][0] = computeValueCur(0, dNN_dx) + 1.0;
        F[gp][2] = computeValueCur(0, dNN_dy);
        F[gp][1] = computeValueCur(1, dNN_dx);
        F[gp][3] = computeValueCur(1, dNN_dy) + 1.0;

        // determinant of the deformation gradient
        detF[gp] = F[gp][0]*F[gp][3] - F[gp][1]*F[gp][2];

        GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &NN(0), &dNN_dx(0), &dNN_dy(0), Jac);

        for(ii=0; ii<nlbfU; ii++)
        {
          N(gp,ii) = NN(ii);
          dN_dx[gp][ii] = dNN_dx(ii);
          dN_dy[gp][ii] = dNN_dy(ii);
        }

        // volume in the deformed configuration
        dvol = gaussweights[gp]*(Jac*thick);
        volume[gp] = dvol;

        // compute (detF)^(1/3)
        dummy = pow(detF[gp],r1d3);

        // compute C matrix
        fact = dvol*dummy;

        detJbar1 += fact;
        for(jj=0;jj<nlbfU;jj++)
        {
          Cmat1(jj) += dN_dx[gp][jj] * fact;
          Cmat2(jj) += dN_dy[gp][jj] * fact;
        }
        //dNc_dx += dN_dx*dvol;
        //dNc_dy += dN_dy*dvol;

        for(ii=0;ii<nlbfU;ii++)
        {
          TI   = 2*ii;
          TIp1 = TI+1;

          fact1 = fact * dN_dx[gp][ii];
          fact2 = fact * dN_dy[gp][ii];

          for(jj=0;jj<nlbfU;jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Cmat3[TI][TJ]     +=  ( fact1 * dN_dx[gp][jj] );
            Cmat3[TI][TJp1]   +=  ( fact1 * dN_dy[gp][jj] );
            Cmat3[TIp1][TJ]   +=  ( fact2 * dN_dx[gp][jj] );
            Cmat3[TIp1][TJp1] +=  ( fact2 * dN_dy[gp][jj] );

            Cmat4[TI][TJ]     +=  ( fact1 * dN_dx[gp][jj] );
            Cmat4[TI][TJp1]   +=  ( fact2 * dN_dx[gp][jj] );
            Cmat4[TIp1][TJ]   +=  ( fact1 * dN_dy[gp][jj] );
            Cmat4[TIp1][TJp1] +=  ( fact2 * dN_dy[gp][jj] );
          }
        }
    }

    //cout << " MAB = " << MAB << '\t' << detJbar1 << endl;

    // Compute the inverse of MASS Matrix of the projected space
    invMAB = 1.0/MAB;

    // compute W matrix W = invMAB * C'

    Wmat1 = invMAB * Cmat1;
    Wmat2 = invMAB * Cmat2;

    for(ii=0;ii<nsize;ii++)
    {
      for(jj=0;jj<nsize;jj++)
      {
        Wmat3[ii][jj] = invMAB * Cmat3[ii][jj];
        Wmat4[ii][jj] = invMAB * Cmat4[ii][jj];
      }
    }

    detJbar2 = invMAB*detJbar1;
    detFbar  = detJbar2;

    for(ii=0;ii<nlbfU;ii++)
    {
      g[0][ii] = Wmat1(ii) / detFbar;
      g[1][ii] = Wmat2(ii) / detFbar ;
    }

    for(ii=0;ii<nsize;ii++)
    {
      for(jj=0;jj<nsize;jj++)
      {
        g2[ii][jj] =  Wmat3[ii][jj] / detFbar ;
        g3[ii][jj] =  Wmat4[ii][jj] / detFbar ;
      }
    }


    // resize local matrices and initialise them to zero
    if(Kuu.rows() != nsize)
    {
      Kuu.resize(nsize, nsize);
      FlocalU.resize(nsize);
    }
    Kuu.setZero();
    FlocalU.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        dvol0 = volume0[gp];
        dvol  = volume[gp];

        // calculate Fbar
        alpha = detFbar/pow(detF[gp], r1d3);

        for(ii=0;ii<4;ii++)
          F[gp][ii] = alpha * F[gp][ii];

        F33 = alpha ;

        //matlib2d_(matDat, F[gp], &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)          return 1;

        dvol *= F33;

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
          for(jj=0;jj<4;jj++)
            cc[ii][jj] *= dvol;
        }

        //==============================================
        // CALCULATE Bbar Matrices
        //==============================================

        for(ii=0;ii<nlbfU;ii++)
        {
           gmN[0][ii] = g[0][ii] - dN_dx[gp][ii];
           gmN[1][ii] = g[1][ii] - dN_dy[gp][ii];

           bb1 = gmN[0][ii]/3.0;
           bb2 = gmN[1][ii]/3.0;

           Bbar[ii][0][0] = dN_dx[gp][ii] + bb1;
           Bbar[ii][1][0] = bb1;
           Bbar[ii][2][0] = bb1;
           Bbar[ii][3][0] = dN_dy[gp][ii];

           Bbar[ii][0][1] = bb2;
           Bbar[ii][1][1] = dN_dy[gp][ii] + bb2;
           Bbar[ii][2][1] = bb2;
           Bbar[ii][3][1] = dN_dx[gp][ii];
       }

        // Calculate Stiffness and Residual
        //==============================================

        NN = N.col(gp);

        acceCur[0] = computeValueDotDotCur(0, NN);
        acceCur[1] = computeValueDotDotCur(1, NN);

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        // 1.) Material Contribution
        //==============================================

        for(ii=0;ii<nlbfU;ii++)
        {
           bb5 = dvol0*N(gp,ii);

           bc[0][0] = (Bbar[ii][0][0] * cc[0][0] + Bbar[ii][1][0] * cc[1][0] + Bbar[ii][2][0] * cc[2][0] + Bbar[ii][3][0] * cc[3][0]);
           bc[0][1] = (Bbar[ii][0][0] * cc[0][1] + Bbar[ii][1][0] * cc[1][1] + Bbar[ii][2][0] * cc[2][1] + Bbar[ii][3][0] * cc[3][1]);
           bc[0][2] = (Bbar[ii][0][0] * cc[0][2] + Bbar[ii][1][0] * cc[1][2] + Bbar[ii][2][0] * cc[2][2] + Bbar[ii][3][0] * cc[3][2]);
           bc[0][3] = (Bbar[ii][0][0] * cc[0][3] + Bbar[ii][1][0] * cc[1][3] + Bbar[ii][2][0] * cc[2][3] + Bbar[ii][3][0] * cc[3][3]);

           bc[1][0] = (Bbar[ii][0][1] * cc[0][0] + Bbar[ii][1][1] * cc[1][0] + Bbar[ii][2][1] * cc[2][0] + Bbar[ii][3][1] * cc[3][0]);
           bc[1][1] = (Bbar[ii][0][1] * cc[0][1] + Bbar[ii][1][1] * cc[1][1] + Bbar[ii][2][1] * cc[2][1] + Bbar[ii][3][1] * cc[3][1]);
           bc[1][2] = (Bbar[ii][0][1] * cc[0][2] + Bbar[ii][1][1] * cc[1][2] + Bbar[ii][2][1] * cc[2][2] + Bbar[ii][3][1] * cc[3][2]);
           bc[1][3] = (Bbar[ii][0][1] * cc[0][3] + Bbar[ii][1][1] * cc[1][3] + Bbar[ii][2][1] * cc[2][3] + Bbar[ii][3][1] * cc[3][3]);

           TI   = 2*ii;
           TIp1 = TI+1;

           FlocalU(TI)   += bb5*force[0] ;
           FlocalU(TIp1) += bb5*force[1] ;

           FlocalU(TI)   -= (Bbar[ii][0][0]*stre[0] + Bbar[ii][1][0]*stre[1] + Bbar[ii][2][0]*stre[2] + Bbar[ii][3][0]*stre[3]) ;
           FlocalU(TIp1) -= (Bbar[ii][0][1]*stre[0] + Bbar[ii][1][1]*stre[1] + Bbar[ii][2][1]*stre[2] + Bbar[ii][3][1]*stre[3]) ;

           for(jj=0;jj<nlbfU;jj++)
           {
              TJ   = 2*jj;
              TJp1 = TJ+1;

              cc3 = N(gp,jj);

              // acceleration term
              acceFact2 = acceFact1*cc3*rho0;

              fact  = bb5*acceFact2;

              Kuu(TI,   TJ)    += fact ;
              Kuu(TIp1, TJp1)  += fact ;

              Kuu(TI, TJ)     +=  (bc[0][0]*Bbar[jj][0][0] + bc[0][1]*Bbar[jj][1][0] + bc[0][2]*Bbar[jj][2][0] + bc[0][3]*Bbar[jj][3][0]) ;
              Kuu(TI, TJp1)   +=  (bc[0][0]*Bbar[jj][0][1] + bc[0][1]*Bbar[jj][1][1] + bc[0][2]*Bbar[jj][2][1] + bc[0][3]*Bbar[jj][3][1]) ;
              Kuu(TIp1, TJ)   +=  (bc[1][0]*Bbar[jj][0][0] + bc[1][1]*Bbar[jj][1][0] + bc[1][2]*Bbar[jj][2][0] + bc[1][3]*Bbar[jj][3][0]) ;
              Kuu(TIp1, TJp1) +=  (bc[1][0]*Bbar[jj][0][1] + bc[1][1]*Bbar[jj][1][1] + bc[1][2]*Bbar[jj][2][1] + bc[1][3]*Bbar[jj][3][1]) ;
           }
        }

        // 2.) Geometric Contribution
        //==============================================

        for(ii=0;ii<nlbfU;ii++)
        {
           bc[0][0] =  Bbar[ii][0][0] * stre[0] + Bbar[ii][3][0] * stre[3] ;
           bc[0][1] =  Bbar[ii][0][0] * stre[3] + Bbar[ii][3][0] * stre[1] ;
           bc[0][2] =  Bbar[ii][1][0] * stre[3] ;
           bc[0][3] =  Bbar[ii][1][0] * stre[1] ;
           bc[0][4] =  Bbar[ii][2][0] * stre[2] ;

           bc[1][0] =  Bbar[ii][0][1] * stre[0] ;
           bc[1][1] =  Bbar[ii][0][1] * stre[3] ;
           bc[1][2] =  Bbar[ii][3][1] * stre[0] + Bbar[ii][1][1] * stre[3] ;
           bc[1][3] =  Bbar[ii][3][1] * stre[3] + Bbar[ii][1][1] * stre[1] ;
           bc[1][4] =  Bbar[ii][2][1] * stre[2] ;

           TI   = 2*ii;
           TIp1 = TI+1;

           for(jj=0;jj<nlbfU;jj++)
           {
              TJ   = 2*jj;
              TJp1 = TJ+1;

              Kuu(TI, TJ)     +=  (bc[0][0] * Bbar[jj][0][0] + bc[0][1] * Bbar[jj][3][0] + bc[0][3] * Bbar[jj][1][0] + bc[0][4] * Bbar[jj][2][0]) ;
              Kuu(TI, TJp1)   +=  (bc[0][0] * Bbar[jj][0][1] + bc[0][2] * Bbar[jj][3][1] + bc[0][3] * Bbar[jj][1][1] + bc[0][4] * Bbar[jj][2][1]) ;
              Kuu(TIp1, TJ)   +=  (bc[1][0] * Bbar[jj][0][0] + bc[1][1] * Bbar[jj][3][0] + bc[1][3] * Bbar[jj][1][0] + bc[1][4] * Bbar[jj][2][0]) ;
              Kuu(TIp1, TJp1) +=  (bc[1][0] * Bbar[jj][0][1] + bc[1][2] * Bbar[jj][3][1] + bc[1][3] * Bbar[jj][1][1] + bc[1][4] * Bbar[jj][2][1]) ;
           }
        }


        // 3.) COMPUTE THE EXTRA TERMS
        //==============================================

        //fact = ( stre[0] + stre[1] + stre[2] )/9.0 ;
        fact = ( stre[0] + stre[1] )/9.0 ;
        fact9 = 3.0 * fact ;

        for(ii=0;ii<nlbfU;ii++)
        {
          fact1 = ( stre[0] * dN_dx[gp][ii] + stre[3] * dN_dy[gp][ii] )/3.0 + fact * gmN[0][ii] ;
          fact2 = ( stre[3] * dN_dx[gp][ii] + stre[1] * dN_dy[gp][ii] )/3.0 + fact * gmN[1][ii] ;

          fact5 = fact * g[0][ii] ;
          fact6 = fact * g[1][ii] ;

          fact7 = fact9 * dN_dx[gp][ii] ;
          fact8 = fact9 * dN_dy[gp][ii] ;

          TI   = 2*ii;
          TIp1 = TI+1;

          for(jj=0;jj<nlbfU;jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            fact3 = ( stre[0] * dN_dx[gp][jj] + stre[3] * dN_dy[gp][jj] )/3.0;
            fact4 = ( stre[3] * dN_dx[gp][jj] + stre[1] * dN_dy[gp][jj] )/3.0;

            Kuu(TI,  TJ)   += (gmN[0][ii] * fact3 + fact1 * gmN[0][jj] - fact5 * g[0][jj] + fact7 * dN_dx[gp][jj]) ;
            Kuu(TI,  TJp1) += (gmN[0][ii] * fact4 + fact1 * gmN[1][jj] - fact5 * g[1][jj] + fact8 * dN_dx[gp][jj]) ;
            Kuu(TIp1,TJ)   += (gmN[1][ii] * fact3 + fact2 * gmN[0][jj] - fact6 * g[0][jj] + fact7 * dN_dy[gp][jj]) ;
            Kuu(TIp1,TJp1) += (gmN[1][ii] * fact4 + fact2 * gmN[1][jj] - fact6 * g[1][jj] + fact8 * dN_dy[gp][jj]) ;
          }
        }

        //
        for(ii=0;ii<nsize;ii++)
        {
           for(jj=0;jj<nsize;jj++)
           {
             Kuu(ii,jj)  +=  fact * ( g2[ii][jj] - 3.0 * g3[ii][jj] );
           }
        }
        //
    } //gp

    //printMatrix(Kuu);  printVector(FlocalU);

    return 0;
}
*/




int BernsteinElem2DElecMechTri22F::calcStiffnessAndResidualFS2(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter)
{
    //F-bar formulation by Ejgudej

    int   err,  isw,  count,  count1, index, ll, ii, jj, kk, gp;
    int   TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  fact, dvol, dvol0, Jac, r1d3 = 1.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, fact1, fact2;
    double  fact3, fact4, fact5, fact6, fact7, fact8, fact9;
    double  stre[4], cc[4][4], bc[3][6], param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];
    double  fact11, fact12, fact13, fact21, fact22, fact23;
    double  fact31, fact32, fact33, fact41, fact42, fact43;


    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;


    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);


    double  F[nGP][9], FF[4], detF[nGP], Fc[4], detFc, F33, volume[nGP], volume0[nGP];
    double  dN_dx[nGP][nlbfU], dN_dy[nGP][nlbfU], dN_dz[nGP][nlbfU];

    MatrixXd  N(nGP,nlbfU);//, dN_dx(nGP,nlbfU), dN_dy(nGP,nlbfU);
    VectorXd  NN(nlbfU), dNN_dx(nlbfU), dNN_dy(nlbfU), dNN_dz(nlbfU);
    double  Jp1d3, Jp1d3bar, dummy, alpha, sii, siid3, siid9;
    double  g2[nsize][nsize], g3[nsize][nsize];
    MatrixXd  Cmat5(nsize,nsize), Wmat5(nsize,nsize), Cmat6(nsize,nsize), Wmat6(nsize,nsize) ;
    MatrixXd  Bbar1(2,4), Dmat(4,4), Bbar2(4,2), BDmat(2,4), BDBmat(2,2);
    MatrixXd  Wbar1(2,5), Wbar2(5,2);


    Jp1d3bar = 0.0;
    elemVolOrig=0.0;
    elemVolCur=0.0;
    dNc_dx.setZero();
    dNc_dy.setZero();

    Cmat5.setZero();
    Wmat5.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &NN(0), &dNN_dx(0), &dNN_dy(0), Jac);

        // volume in the original configuration
        dvol0 = gaussweights[gp]*Jac;
        elemVolOrig += dvol0;
        volume0[gp] = dvol0;

        F[gp][0] = computeValueCur(0, dNN_dx) + 1.0;
        F[gp][2] = computeValueCur(0, dNN_dy);

        F[gp][1] = computeValueCur(1, dNN_dx);
        F[gp][3] = computeValueCur(1, dNN_dy) + 1.0;

        // determinant of the deformation gradient
        detF[gp] = F[gp][0]*F[gp][3] - F[gp][1]*F[gp][2] ;

        GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &NN(0), &dNN_dx(0), &dNN_dy(0), Jac);

        // volume in the deformed configuration
        dvol = gaussweights[gp]*Jac;
        elemVolCur += dvol;
        volume[gp] = dvol;

        // compute (J)^(1/3)
        Jp1d3 = pow(detF[gp],r1d3);

        Jp1d3bar += Jp1d3*dvol0;

        // compute C matrix
        fact = Jp1d3*dvol0;

        for(ii=0; ii<nlbfU; ii++)
        {
          N(gp,ii) = NN(ii);

          dN_dx[gp][ii] = dNN_dx(ii);
          dN_dy[gp][ii] = dNN_dy(ii);

          fact1 = fact * dNN_dx(ii);
          fact2 = fact * dNN_dy(ii);

          dNc_dx(ii) += fact1;
          dNc_dy(ii) += fact2;

          TI   = 2*ii;
          TIp1 = TI+1;

          for(jj=0;jj<nlbfU;jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Cmat5(TI,   TJ)    +=  ( fact1 * dNN_dx(jj) );
            Cmat5(TI,   TJp1)  +=  ( fact1 * dNN_dy(jj) );

            Cmat5(TIp1, TJ)    +=  ( fact2 * dNN_dx(jj) );
            Cmat5(TIp1, TJp1)  +=  ( fact2 * dNN_dy(jj) );
          }
        }
    }

    //cout << " elemVol = " << elemVolOrig << '\t' << elemVolCur << '\t' << Jp1d3bar << endl;

    // compute pis
    Jp1d3bar = Jp1d3bar / elemVolOrig;
    dNc_dx   = dNc_dx/elemVolOrig;
    dNc_dy   = dNc_dy/elemVolOrig;
    Cmat5    = Cmat5/elemVolOrig;

    // compute barred quantities
    dNc_dx = dNc_dx/Jp1d3bar;
    dNc_dy = dNc_dy/Jp1d3bar;
    Cmat5  = Cmat5/Jp1d3bar;

    //printMatrix(Cmat5);

    // resize local matrices and initialise them to zero
    // resize local matrices and initialise them to zero
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)     Kuu.resize(nsize, nsize);
    Kuu.setZero();
    if(Kuf.rows() < nsize || Kuf.cols() < nlbfF)     Kuf.resize(nsize, nlbfF);
    Kuf.setZero();
    if(Kfu.rows() < nlbfF || Kfu.cols() < nsize)     Kfu.resize(nlbfF, nsize);
    Kfu.setZero();
    if(Kff.rows() < nlbfF || Kff.cols() < nlbfF)     Kff.resize(nlbfF, nlbfF);
    Kff.setZero();

    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalF.rows() < nlbfF)   FlocalF.resize(nlbfF);
    FlocalF.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        dvol0 = volume0[gp];
        dvol  = volume[gp];

        // calculate Fbar
        alpha = Jp1d3bar/pow(detF[gp], r1d3);

        for(ii=0;ii<4;ii++)
          FF[ii] = alpha * F[gp][ii];

        F33 = alpha;

        //cout << alpha << '\t' << dvol << endl;

        //matlib2d_(matDat, FF, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        dvol *= F33;

        //cout << stre[0] << '\t' << stre[1] << '\t' << stre[2] << '\t' << stre[3] << '\t' << stre[4] << '\t' << stre[5] << endl;

        if(err !=0)          return 1;

        for(ii=0; ii<4; ii++)
        {
          stre[ii] *= dvol;
          for(jj=0; jj<4; jj++)
          {
            Dmat(ii,jj) = cc[ii][jj]*dvol;
          }
        }

        //cout << " AlgoType = " << AlgoType << endl;
        //cout << stre[0] << '\t' << stre[1] << '\t' << stre[2] << '\t' << stre[3] << endl;

        // Calculate Stiffness and Residual
        //==============================================

        NN = N.col(gp);

        acceCur[0] = computeValueDotDotCur(0, NN);
        acceCur[1] = computeValueDotDotCur(1, NN);

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        // 1.) Material Contribution
        //==============================================

        for(ii=0;ii<nlbfU;ii++)
        {
           bb5 = dvol0*N(gp,ii);

           TI   = 2*ii;
           TIp1 = TI+1;

           bb1 = (dNc_dx[ii] - dN_dx[gp][ii])/3.0;
           bb2 = (dNc_dy[ii] - dN_dy[gp][ii])/3.0;

           Bbar1(0,0) = bb1+dN_dx[gp][ii];
           Bbar1(0,1) = bb1;
           Bbar1(0,2) = bb1;
           Bbar1(0,3) = dN_dy[gp][ii];

           Bbar1(1,0) = bb2;
           Bbar1(1,1) = bb2+dN_dy[gp][ii];
           Bbar1(1,2) = bb2;
           Bbar1(1,3) = dN_dx[gp][ii];

           BDmat = Bbar1*Dmat;

           FlocalU(TI)   += bb5*force[0] ;
           FlocalU(TIp1) += bb5*force[1] ;

           FlocalU(TI)   -= (Bbar1(0,0)*stre[0] + Bbar1(0,1)*stre[1] + Bbar1(0,2)*stre[2] + Bbar1(0,3)*stre[3] );
           FlocalU(TIp1) -= (Bbar1(1,0)*stre[0] + Bbar1(1,1)*stre[1] + Bbar1(1,2)*stre[2] + Bbar1(1,3)*stre[3] );

           Wbar1(0,0) = Bbar1(0,0)*stre[0] + Bbar1(0,3)*stre[3];
           Wbar1(0,1) = Bbar1(0,1)*stre[3];
           Wbar1(0,2) = Bbar1(0,0)*stre[3] + Bbar1(0,3)*stre[1];
           Wbar1(0,3) = Bbar1(0,1)*stre[1];
           Wbar1(0,4) = Bbar1(0,2)*stre[2]; 

           Wbar1(1,0) = Bbar1(1,0)*stre[0]; 
           Wbar1(1,1) = Bbar1(1,3)*stre[0] + Bbar1(1,1)*stre[3];
           Wbar1(1,2) = Bbar1(1,0)*stre[3];
           Wbar1(1,3) = Bbar1(1,3)*stre[3] + Bbar1(1,1)*stre[1];
           Wbar1(1,4) = Bbar1(1,2)*stre[2];

           for(jj=0;jj<nlbfU;jj++)
           {
              cc1 = (dNc_dx[jj] - dN_dx[gp][jj])/3.0;
              cc2 = (dNc_dy[jj] - dN_dy[gp][jj])/3.0;

              Bbar2(0,0) = cc1+dN_dx[gp][jj];
              Bbar2(1,0) = cc1;
              Bbar2(2,0) = cc1;
              Bbar2(3,0) = dN_dy[gp][jj];

              Bbar2(0,1) = cc2;
              Bbar2(1,1) = cc2+dN_dy[gp][jj];
              Bbar2(2,1) = cc2;
              Bbar2(3,1) = dN_dx[gp][jj];

              TJ   = 2*jj;
              TJp1 = TJ+1;

              cc5 = N(gp,jj);

              // acceleration term
              acceFact2 = acceFact1*cc5*rho0;

              fact  = bb5*acceFact2;

              Kuu(TI,   TJ)    +=  fact;
              Kuu(TIp1, TJp1)  +=  fact;

              // material contribution
              BDBmat = (af*BDmat)*Bbar2;

              Kuu(TI, TJ)     +=  BDBmat(0,0) ;
              Kuu(TI, TJp1)   +=  BDBmat(0,1) ;

              Kuu(TIp1, TJ)   +=  BDBmat(1,0) ;
              Kuu(TIp1, TJp1) +=  BDBmat(1,1) ;

              // Geometric (stress) contribution
              Wbar2(0,0) = Bbar2(0,0);
              Wbar2(1,0) = 0.0;
              Wbar2(2,0) = Bbar2(3,0);
              Wbar2(3,0) = Bbar2(1,0);
              Wbar2(4,0) = Bbar2(2,0); 

              Wbar2(0,1) = Bbar2(0,1); 
              Wbar2(1,1) = Bbar2(3,1);
              Wbar2(2,1) = 0.0;
              Wbar2(3,1) = Bbar2(1,1);
              Wbar2(4,1) = Bbar2(2,1);

              BDBmat = (af*Wbar1)*Wbar2;

              Kuu(TI, TJ)     +=  BDBmat(0,0) ;
              Kuu(TI, TJp1)   +=  BDBmat(0,1) ;

              Kuu(TIp1, TJ)   +=  BDBmat(1,0) ;
              Kuu(TIp1, TJp1) +=  BDBmat(1,1) ;
           }
        }


        // COMPUTE THE EXTRA TERMS
        //==============================================
        sii = stre[0] + stre[1] + stre[2] ;
        siid3 = sii/3.0;
        siid9 = sii/9.0;

        for(ii=0;ii<nlbfU;ii++)
        {
           bb5 = dvol0*N(gp,ii);

           TI   = 2*ii;
           TIp1 = TI+1;

           bb1 = (dNc_dx[ii] - dN_dx[gp][ii])/3.0;
           bb2 = (dNc_dy[ii] - dN_dy[gp][ii])/3.0;

           fact11 = dN_dx[gp][ii]*stre[0] + dN_dy[gp][ii]*stre[3] ;
           fact12 = dN_dx[gp][ii]*stre[3] + dN_dy[gp][ii]*stre[1] ;

           fact21 = bb1*stre[0] + bb2*stre[3] ;
           fact22 = bb1*stre[3] + bb2*stre[1] ;

           for(jj=0;jj<nlbfU;jj++)
           {
              TJ   = 2*jj;
              TJp1 = TJ+1;

              cc1 = (dNc_dx[jj] - dN_dx[gp][jj])/3.0;
              cc2 = (dNc_dy[jj] - dN_dy[gp][jj])/3.0;

              fact31 = af*(dN_dx[gp][jj]*stre[0] + dN_dy[gp][jj]*stre[3]) ;
              fact32 = af*(dN_dx[gp][jj]*stre[3] + dN_dy[gp][jj]*stre[1]) ;

              fact41 = af*(cc1*stre[0] + cc2*stre[3] ) ;
              fact42 = af*(cc1*stre[3] + cc2*stre[1] ) ;

              // term 3 and term 4
              Kuu(TI, TJ)     +=  fact11*cc1 + bb1*fact31 ;
              Kuu(TI, TJp1)   +=  fact11*cc2 + bb1*fact32 ;

              Kuu(TIp1, TJ)   +=  fact12*cc1 + bb2*fact31 ;
              Kuu(TIp1, TJp1) +=  fact12*cc2 + bb2*fact32 ;

              //
              // term 5
              cc1 *= (af*sii);
              cc2 *= (af*sii);

              Kuu(TI, TJ)     +=  bb1*cc1 ;
              Kuu(TI, TJp1)   +=  bb1*cc2 ;

              Kuu(TIp1, TJ)   +=  bb2*cc1 ;
              Kuu(TIp1, TJp1) +=  bb2*cc2 ;

              // term 6
              fact = af*siid9;

              Kuu(TI,   TJ)   +=  fact*Cmat5(TI,   TJ) ;
              Kuu(TI,   TJp1) +=  fact*Cmat5(TI,   TJp1) ;

              Kuu(TIp1, TJ)   +=  fact*Cmat5(TIp1, TJ) ;
              Kuu(TIp1, TJp1) +=  fact*Cmat5(TIp1, TJp1) ;
              //

              // term 7
              cc1 = af*siid9*dNc_dx(jj);
              cc2 = af*siid9*dNc_dy(jj);

              Kuu(TI, TJ)     -=  dNc_dx(ii)*cc1 ;
              Kuu(TI, TJp1)   -=  dNc_dx(ii)*cc2 ;

              Kuu(TIp1, TJ)   -=  dNc_dy(ii)*cc1 ;
              Kuu(TIp1, TJp1) -=  dNc_dy(ii)*cc2 ;

              //
              // term 8
              fact = af*siid3;

              Kuu(TI,   TJ)   +=  fact*dN_dx[gp][ii]*dN_dx[gp][jj] ;
              Kuu(TI,   TJp1) +=  fact*dN_dy[gp][ii]*dN_dx[gp][jj] ;

              Kuu(TIp1, TJ)   +=  fact*dN_dx[gp][ii]*dN_dy[gp][jj] ;
              Kuu(TIp1, TJp1) +=  fact*dN_dy[gp][ii]*dN_dy[gp][jj] ;

              // term 9
              fact = af*siid3;

              Kuu(TI,   TJ)   -=  fact*Cmat5(TI,   TJ) ;
              Kuu(TI,   TJp1) -=  fact*Cmat5(TIp1, TJ) ;

              Kuu(TIp1, TJ)   -=  fact*Cmat5(TI,   TJp1) ;
              Kuu(TIp1, TJp1) -=  fact*Cmat5(TIp1, TJp1) ;
              //
           }
        }
    } //gp

    //printMatrix(Kuu);
    //printVector(FlocalU);
    return 0;
}



int BernsteinElem2DElecMechTri22F::calcResidualFS2(VectorXd& FlocalU)
{
    //F-bar formulation by Ejgudej

    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;

    double  F33, fact, dvol, dvol0, Jac, r1d3 = 1.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, fact1, fact2;
    double  fact3, fact4, fact5, fact6, fact7, fact8, fact9;
    double  stre[4], cc[4][4], param[2], bforce[2], force[2];

    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);


    double  F[nGP][4], detF[nGP], volume[nGP], volume0[nGP];
    double  Bbar[2][4], dN_dx[nGP][nlbfU], dN_dy[nGP][nlbfU];
    double  Fc[4], detFc, Jp1d3, Jp1d3bar, dummy, alpha;

    MatrixXd  N(nGP,nlbfU);//, dN_dx(nGP,nlbfU), dN_dy(nGP,nlbfU);
    VectorXd  NN(nlbfU), dNN_dx(nlbfU), dNN_dy(nlbfU);

    Jp1d3bar = 0.0;
    elemVolOrig=0.0;
    elemVolCur=0.0;

    dNc_dx.setZero();
    dNc_dy.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &NN(0), &dNN_dx(0), &dNN_dy(0), Jac);

        // volume in the original configuration
        dvol0 = gaussweights[gp]*(Jac*thick);
        elemVolOrig += dvol0;
        volume0[gp] = dvol0;

        F[gp][0] = computeValue(0, dNN_dx) + 1.0;
        F[gp][2] = computeValue(0, dNN_dy);
        F[gp][1] = computeValue(1, dNN_dx);
        F[gp][3] = computeValue(1, dNN_dy) + 1.0;

        // determinant of the deformation gradient
        detF[gp] = F[gp][0]*F[gp][3] - F[gp][1]*F[gp][2];

        GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &NN(0), &dNN_dx(0), &dNN_dy(0), Jac);
        // volume in the deformed configuration
        dvol = gaussweights[gp]*Jac;
        elemVolCur += dvol;
        volume[gp] = dvol;

        // compute (J)^(1/3)
        Jp1d3 = pow(detF[gp],r1d3);

        Jp1d3bar += Jp1d3*dvol0;

        // compute C matrix
        fact = Jp1d3*dvol0;

        for(ii=0; ii<nlbfU; ii++)
        {
          N(gp,ii) = NN(ii);

          dN_dx[gp][ii] = dNN_dx(ii);
          dN_dy[gp][ii] = dNN_dy(ii);

          dNc_dx(ii) += dNN_dx(ii) * fact;
          dNc_dy(ii) += dNN_dy(ii) * fact;
        }
    }

    //cout << " MAB = " << MAB << '\t' << detJbar1 << endl;

    // computing pis
    Jp1d3bar = Jp1d3bar/elemVolOrig;
    dNc_dx   = dNc_dx/elemVolOrig;
    dNc_dy   = dNc_dy/elemVolOrig;

    // compute barred quantities
    dNc_dx = dNc_dx/Jp1d3bar;
    dNc_dy = dNc_dy/Jp1d3bar;

    // resize local matrices and initialise them to zero
    if(FlocalU.rows() != nsize)
      FlocalU.resize(nsize);
    FlocalU.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        dvol0 = volume0[gp];
        dvol  = volume[gp];

        // calculate Fbar
        alpha = Jp1d3bar/pow(detF[gp], r1d3);

        for(ii=0;ii<4;ii++)
          F[gp][ii] = alpha * F[gp][ii];

        F33 = alpha ;

        //matlib2d_(matDat, F[gp], &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)          return 1;

        dvol *= F33;

        for(ii=0;ii<4;ii++)
        {
          stre[ii] *= dvol;
        }

        //==============================================
        // CALCULATE Bbar Matrices
        //==============================================

        // Calculate Residual
        //====================

        NN = N.col(gp);

        // body force terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        for(ii=0;ii<nlbfU;ii++)
        {
           bb5 = dvol0*N(gp,ii);

           TI   = 2*ii;
           TIp1 = TI+1;

           fact1 = (dNc_dx[ii] - dN_dx[gp][ii])/3.0;
           fact2 = (dNc_dy[ii] - dN_dy[gp][ii])/3.0;

           Bbar[0][0] = fact1+dN_dx[gp][ii];
           Bbar[0][1] = fact1;
           Bbar[0][2] = fact1;
           Bbar[0][3] = dN_dy[gp][ii];

           Bbar[1][0] = fact2;
           Bbar[1][1] = fact2+dN_dy[gp][ii];
           Bbar[1][2] = fact2;
           Bbar[1][3] = dN_dx[gp][ii];

           FlocalU(TI)   += bb5*force[0] ;
           FlocalU(TIp1) += bb5*force[1] ;

           FlocalU(TI)   -= (Bbar[0][0]*stre[0] + Bbar[0][1]*stre[1] + Bbar[0][2]*stre[2] + Bbar[0][3]*stre[3] ) ;
           FlocalU(TIp1) -= (Bbar[1][0]*stre[0] + Bbar[1][1]*stre[1] + Bbar[1][2]*stre[2] + Bbar[1][3]*stre[3] ) ;
        }
    } //gp

    //printMatrix(Kuu);  printVector(FlocalU);

    return 0;
}





int BernsteinElem2DElecMechTri22F::calcLoadVector(VectorXd& FlocalU)
{
  if( NeumannBCs.size() == 0)
    return 0;


  int side, dir, ii, jj, nn, TI, TIp1, TJ, TJp1, gp, nGP=3;
  int face_node_map[3][3]={{0,3,1},{1,4,2},{2,5,0}};
  double  specVal, xNode[3], yNode[3], xx, yy, param[2], trac[2], normal[2], tang[2], N[3];
  double  Jac, dvol;

  vector<double>  gausspoints1, gaussweights;
  getGaussPoints1D(nGP, gausspoints1, gaussweights);

  //printVector(FlocalU);  printf("\n\n\n");

  // side 0 - 1-4-2
  // side 1 - 2-5-3
  // side 2 - 3-6-1

  // loop over all the sides and compute force vector due to applied traction
  for(nn=0; nn<NeumannBCs.size(); nn++)
  {
    //printVector(NeumannBCs[nn]);

    side    = NeumannBCs[nn][0]; // edge/face number
    dir     = NeumannBCs[nn][1]; // element number
    specVal = NeumannBCs[nn][2]; // specified value
    specVal *= timeFunction[0].prop; // load factor

    for(ii=0; ii<3; ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[face_node_map[side][ii]]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[face_node_map[side][ii]]][1];
    }

    //elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
      param[0] = gausspoints1[gp];

      //BernsteinBasisFunsEdge2D(degree, param, xNode, yNode, N, normal, Jac);

      dvol = gaussweights[gp]*(Jac*thick);
      elemVol += dvol;

      //cout << " normal = " << normal[0] << '\t' << normal[1] << '\t' << Jac << '\t' << dvol << endl;
      //cout << " Funcs = " << N[0] << '\t' << N[1] << '\t' << N[2] << endl;

        xx = yy= 0.0;
        for(ii=0;ii<3;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        if(axsy)
          dvol *= 2.0*PI*xx;

        tang[0] = -normal[1];
        tang[1] =  normal[0];

        if(dir == 0)
        {
          trac[0] = specVal*normal[0]*dvol;
          trac[1] = specVal*normal[1]*dvol;
        }
        if(dir == 1)
        {
          trac[0] = specVal*tang[0]*dvol;
          trac[1] = specVal*tang[1]*dvol;
        }

        //cout << " trac = " << trac[0] << '\t' << trac[1] << endl;

        for(ii=0;ii<3;ii++)
        {
          TI   = face_node_map[side][ii]*ndof;
          TIp1 = TI+1;

          FlocalU(TI)   += N[ii]*trac[0];
          FlocalU(TIp1) += N[ii]*trac[1];
        }
    } // for(gp=0; gp<nGP; gp++)
    //cout << " Volume = " << elemVol << endl;
  }

  //printVector(FlocalU);  printf("\n\n\n");

  return 0;
}



int BernsteinElem2DElecMechTri22F::calcInternalForces()
{

  return 0;
}



void BernsteinElem2DElecMechTri22F::elementContourplot(int vartype, int varindex, int index)
{
   double outval[50];

   switch(vartype)
   {
       case 1:  // plot total strain
       case 2:  // plot elastic strain
       case 3:  // plot plastic strain

                projectStrain(true, vartype, varindex, index, outval);

              break;

       case 4:  // plot stress

                //intVar2 = intVar1;
                projectStress(true, vartype, varindex, index, outval);

              break;

       case 5:  // plot element internal variables

                projectInternalVariable(true, vartype, varindex, index, outval);

              break;

       default:

              cout  << " Invalid Variable Type to project in 'BernsteinElem2DSolidTria6Node::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
      vals2project[0] += outval[ii];

    vals2project[0] /= nGP;

    return;
}


void BernsteinElem2DElecMechTri22F::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
{
   double outval[50];

   switch(vartype)
   {
       case 1:  // plot total strain
       case 2:  // plot elastic strain
       case 3:  // plot plastic strain

                projectStrain(extrapolateFlag, vartype, varindex, index, outval);

              break;

       case 4:  // plot stress

                //intVar2 = intVar1;
                projectStress(extrapolateFlag, vartype, varindex, index, outval);

              break;

       case 5:  // plot element internal variables

                projectInternalVariable(extrapolateFlag, vartype, varindex, index, outval);

              break;

       default:

              cout  << "           Invalid Variable Type to project " << endl;
              break;
    }


    assert(vals2project.size() == npElem);

    // project from quadrature points to nodes/control points
    stressrecovery_extrapolate_Triangle(degree, outval, &vals2project[0]);

    //printVector(vals2project);

    return;
}


void BernsteinElem2DElecMechTri22F::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    int   err,  isw,  count,  count1, ll = 0, ii, jj, gp;

    double  detF, detFc, trFc, dvol, Jac, fact, dvol0, param[3], r1d3 = 1.0/3.0;
    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), Bmat(3,9), Amat(3,3);
    VectorXd  stre(9),  elecDisp(3),  elecField(3);
    F.setZero();
    Fn.setZero();
    stre.setZero();

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU);
    VectorXd  Nuc(nlbfU), dNuc_dx(nlbfU), dNuc_dy(nlbfU);

    double dt   = myTime.dt;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    vector<double>  gausspoints1, gausspoints2, gaussweights;

    if(AlgoType == 1)
    {
      getGaussPointsTriangle(1, gausspoints1, gausspoints2, gaussweights);

      param[0] = gausspoints1[0];
      param[1] = gausspoints2[0];

      GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nuc(0), &dNuc_dx(0), &dNuc_dy(0), Jac);

      F(0,0) = computeValue(0, dNuc_dx) + 1.0;            // eps_xx
      F(0,1) = computeValue(1, dNuc_dx);                  // eps_yx
      F(1,0) = computeValue(0, dNuc_dy);                  // eps_xy
      F(1,1) = computeValue(1, dNuc_dy) + 1.0;            // eps_yy

      detFc = F(0,0)*F(1,1) - F(0,1)*F(1,0);
      trFc  = F(0,0) + F(1,1);

      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);
    } //if(AlgoType == 1)
    else
    {
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

      elemVolOrig=0.0;
      detFc = 0.0;
      trFc  = 0.0;

      for(gp=0; gp<nGP; gp++)
      {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nuc(0), &dNuc_dx(0), &dNuc_dy(0), Jac);

          // volume in the original configuration
          dvol0 = gaussweights[gp]*(Jac*thick);
          elemVolOrig += dvol0;

          F(0,0) = computeValue(0, dNuc_dx) + 1.0;       // eps_xx
          F(0,1) = computeValue(1, dNuc_dx);             // eps_yx
          F(1,0) = computeValue(0, dNuc_dy);             // eps_xy
          F(1,1) = computeValue(1, dNuc_dy) + 1.0;       // eps_yy

          // determinant of the deformation gradient (for Fbar)
          detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

          // trace of the small-strain tensor (for Bbar)
          trFc  += (F(0,0) + F(1,1))*dvol0;
          detFc += pow(detF, r1d3)*dvol0;
        }
        trFc  /= elemVolOrig;
        detFc /= elemVolOrig;
    }//if(AlgoType == 1)


    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

        F(0,0) = computeValue(0, dNu_dx) + 1.0;       // eps_xx
        F(0,1) = computeValue(1, dNu_dx);             // eps_yx
        F(1,0) = computeValue(0, dNu_dy);             // eps_xy
        F(1,1) = computeValue(1, dNu_dy) + 1.0;       // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        // modify detF (tr(F)) with detFc (tr(Fc))

        if(sss == 1)  // plane stress
        {
          if(finite)
            F(2,2) = 1.0/sqrt(detF);
          else
            F(2,2) = 3.0 - F(0,0) - F(1,1);
        }
        else if(sss == 2)    // plane strain
        {
          if(finite)
          {
            if(AlgoType == 1)
            {
              fact = sqrt(detFc/detF);
              F *= fact ;

              F(2,2) = 1.0;
            }
            else
            {
              fact = detFc/pow(detF, 1.0/3.0);
              F *= fact ;

              F(2,2) = fact;
            }
          }
          else
          {
            //fact = 0.5 * (trFc - F[0] - F[3]) ;
            //F[0] = F[0] + fact ;
            //F[3] = F[3] + fact ;

            fact = (trFc - F(0,0) - F(1,1))/3.0;

            F(0,0) += fact ;
            F(1,1) += fact ;
            F(2,2)  = 1.0 + fact;
          }
        }

        double  dummy = MatlData->matData[2];
        MatlData->matData[2] = 0.0;

        pres = 0.0;
        MatlData->computeStressAndTangent(false, sss, Fn, F, elecField, pres, stre, elecDisp, Cmat, Bmat, Amat, ivar, gp, dt);
        count++;
        ll += nivGP;

        MatlData->matData[2] = dummy;

        if(varindex < 9)
          outval[gp] = stre[varindex];
        else if(varindex == 9)
          outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*stre[1]*stre[1])/2.0);
        else if(varindex == 10)
          outval[gp] = (stre[0]+stre[4]+stre[8])/3.0;
        else
          outval[gp] = 0.0;
    }//gp

    return;
}



void BernsteinElem2DElecMechTri22F::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    return;
}



void BernsteinElem2DElecMechTri22F::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp
    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}

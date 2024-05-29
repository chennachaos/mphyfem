
#include "BernsteinElem2DElecMechTri210.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "BasisFunctionsBernstein.h"


using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



BernsteinElem2DElecMechTri210::BernsteinElem2DElecMechTri210()
{
  ndim   = 2;
  degree = 2;
  npElem = 6;
  nlbfU   = 6;
  ndof   = 2;
  nsize  = npElem*ndof;

  dNc_dx.resize(nlbfU);
  dNc_dy.resize(nlbfU);

  pres = presPrev = 0.0;
  AlgoType = 2;
}


BernsteinElem2DElecMechTri210::~BernsteinElem2DElecMechTri210()
{
}


void BernsteinElem2DElecMechTri210::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double BernsteinElem2DElecMechTri210::computeVolume(bool configflag)
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

          if(configflag)
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

  return elemVol;
}



int BernsteinElem2DElecMechTri210::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != nsize) || (Mlocal.cols() != nsize) )
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = elmDat[5]*computeVolume(true)/6.0;

      for(int ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
  else
  {
    int  ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double fact, dvol0, Jac, bb1, cc1, param[2];
    double rho0 = elmDat[5] ;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int nGPt=7;
    getGaussPointsTriangle(nGPt, gausspoints1, gausspoints2, gaussweights);

    elemVolOrig=0.0;
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

        elemVolOrig += dvol0;

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
    cout << " elemVolOrig = " << elemVolOrig << endl;
    //printMatrix(Klocal);  printf("\n\n\n");
  }

  return 0;
}


int BernsteinElem2DElecMechTri210::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != 1) || (Mlocal.cols() != 1) )
      Mlocal.resize(1, 1);

    Mlocal(0,0) = computeVolume(true);

    return 0;
}


double  BernsteinElem2DElecMechTri210::calcCriticalTimeStep(bool flag)
{
    int  ii;

    double xNode[6], yNode[6];
    for(ii=0;ii<npElem;ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

      xNode[ii] = GeomData->NodePosNew[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosNew[nodeNums[ii]][1];
    }

    double  xd, yd, charlen=1.0e10;
    int  edge_nodes[3][2] = { {0,1}, {1,2}, {2,0}};

    for(ii=0; ii<3; ii++)
    {
      xd = xNode[edge_nodes[ii][0]] - xNode[edge_nodes[ii][1]];
      yd = yNode[edge_nodes[ii][0]] - yNode[edge_nodes[ii][1]]; 

      charlen = min(charlen, xd*xd+yd*yd);
    }
    charlen = 0.5*sqrt(charlen);

    double K    = matDat[0] ;
    double mu   = matDat[1] ;
    double rho  = elmDat[5] ;

    double  wave_speed = sqrt((K+4.0*mu/3.0)/rho);
    //double  wave_speed = sqrt((4.0*mu/3.0)/rho);

    double  dtCric = charlen/wave_speed;

    return  dtCric;
}



int BernsteinElem2DElecMechTri210::calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, volstrain, pbar, r2d3 = 2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, dUdJ, d2UdJ2;
    double  stre[4], streDev[4], cc[4][4], bc[2][4], param[2], bforce[2], force[2];
    double  Idev[4][4], D11[4][4], cctmp[4][4];
    double  veloCur[2], acceCur[2], sig[2];
    Idev2D(Idev);

    double BULK = matDat[0];
    double eps  = 0.0; if(!SolnData->TRULY_INCOMPRESSIBLE) eps=1.0/BULK;
    double mu   = matDat[1] ;
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

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal1.rows() <= nsize)   Flocal1.resize(nsize);
    Flocal1.setZero();
    if(Flocal2.rows() <= 1)   Flocal2.resize(1);
    Flocal2.setZero();
    if(Kup.rows() <= nsize || Kup.cols() <= 1)   Kup.resize(nsize, 1);
    Kup.setZero();
    if(Kpu.rows() <= 1 || Kpu.cols() <= nsize)   Kpu.resize(1, nsize);
    Kpu.setZero();
    if(Kpp.rows() <= 1 || Kpp.cols() <= 1)   Kpp.resize(1, 1);
    Kpp.setZero();
    if(Kuu.rows() <= nsize || Kuu.cols() <= nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    //Stokes2DEx1  analy;
    //Elasticity2DMixedEx1  analy(matDat[0], matDat[1]);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        volstrain = F[0] + F[3] - 2.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);


        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)
          dvol *= F33;

        //printf("strains   ... %5d \t %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", elenum, F[0]-1.0, F[1], F[2], F[3]-1.0);
        //printf("stresses  ... %5d \t %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", elenum, stre[0], stre[1], stre[2], stre[3]);
        //for(ii=0; ii<4; ii++)
          //printf("%14.12f \t %14.12f \t %14.12f \t %14.12f \n", cc[ii][0], cc[ii][1], cc[ii][2], cc[ii][3]);
        //printf("\n\n\n");

        // evaluate pressure at the quadrature points
        // for this elements its a single value
        pres = SolnData->var3[elenum];

        if(SolnData->TRULY_INCOMPRESSIBLE)
        {
          stre[0] += pres;
          stre[1] += pres;
          stre[2] += pres;

          for(ii=0;ii<4;ii++)
          {
            for(jj=0;jj<4;jj++)
            {
              D11[ii][jj] = cc[ii][jj];
            }
          }
          if(finite)
          {
            fact = -2.0*pres;
            fact1 = -pres;
  
            for(ii=0;ii<3;ii++)
            {
              for(jj=0;jj<3;jj++)
                D11[ii][jj] -= fact1;
  
              D11[ii][ii] += fact;
            }
            D11[3][3] += 0.5*fact;
          }
        }
        else
        {
          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          streDev[0] = stre[0] - pbar;
          streDev[1] = stre[1] - pbar;
          streDev[2] = stre[2] - pbar;
          streDev[3] = stre[3];

          stre[0] = streDev[0] + pres;
          stre[1] = streDev[1] + pres;
          stre[2] = streDev[2] + pres;

          //
          for(ii=0;ii<4;ii++)
          {
            for(jj=0;jj<4;jj++)
            {
              cctmp[ii][jj] = 0.0;
              for(mm=0;mm<4;mm++)
                cctmp[ii][jj] += Idev[ii][mm] * cc[mm][jj];
            }
          }
  
          for(ii=0;ii<4;ii++)
          {
            for(jj=0;jj<4;jj++)
            {
              D11[ii][jj] = 0.0;
              for(mm=0;mm<4;mm++)
                D11[ii][jj] += cctmp[ii][mm] * Idev[mm][jj];
            }
          }

          if(finite)
          {
            fact = 2.0*(pbar - pres);
            fact1 = (r2d3*pbar - pres);
  
            for(ii=0;ii<3;ii++)
            {
              for(jj=0;jj<4;jj++)
              {
                D11[ii][jj] -= r2d3*streDev[jj] ;
                D11[jj][ii] -= r2d3*streDev[jj] ;
              }
              for(jj=0;jj<3;jj++)
                D11[ii][jj] -= fact1;
  
              D11[ii][ii] += fact;
            }
            D11[3][3] += 0.5*fact;
          }
          //
          /*
          for(ii=0;ii<4;ii++)
          {
            for(jj=0;jj<4;jj++)
            {
              D11[ii][jj] = cc[ii][jj];
            }
          }
          if(finite)
          {
            fact = -2.0*pres;
            fact1 = -pres;
  
            for(ii=0;ii<3;ii++)
            {
              for(jj=0;jj<3;jj++)
                D11[ii][jj] -= fact1;
  
              D11[ii][ii] += fact;
            }
            D11[3][3] += 0.5*fact;
          }
          */
        }

        //for(ii=0;ii<4;ii++)
          //cout << D11[ii][0] << '\t' << D11[ii][1] << '\t' << D11[ii][2] << '\t' << D11[ii][3] << endl;
        //printf("\n\n\n");

        // Calculate Stiffness and Residual
        //==============================================

        xx = yy = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        //force[0] = analy.computeForce(0, xx, yy, 0.0) ;
        //force[1] = analy.computeForce(1, xx, yy, 0.0) ;

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        if(finite)
        {
          if( matDat[2] == 6 )
          {
            dUdJ   = LAMBDA*log(detF)/detF;
            d2UdJ2 = LAMBDA*(1.0-log(detF))/detF/detF;
          }
          else
          {
            dUdJ   = 0.5*BULK*(detF-1.0/detF);
            d2UdJ2 = 0.5*BULK*(1.0+1.0/detF/detF);
          }
          d2UdJ2 /= BULK;
        }
        else
        {
          d2UdJ2 = 1.0;
        }

        for(ii=0; ii<nlbfU; ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          bc[0][0] = bb1 * D11[0][0] + bb2 * D11[3][0];
          bc[0][1] = bb1 * D11[0][1] + bb2 * D11[3][1];
          bc[0][2] = bb1 * D11[0][3] + bb2 * D11[3][3];

          bc[1][0] = bb2 * D11[1][0] + bb1 * D11[3][0];
          bc[1][1] = bb2 * D11[1][1] + bb1 * D11[3][1];
          bc[1][2] = bb2 * D11[1][3] + bb1 * D11[3][3];

          TI   = 2*ii;
          TIp1 = TI+1;

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[3] + bb2*stre[1];

          Flocal1(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal1(TIp1) += (bb5*force[1] - sig[1]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = N[jj];

            TJ   = 2*jj;
            TJp1 = TJ+1;

            // acceleration term
            acceFact2 = acceFact1*cc3*rho0;

            fact  = bb5*acceFact2;

            // material Stiffness contribution
            fact  += af*(sig[0]*cc1+sig[1]*cc2)*FiniteFact;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;

            Kuu(TI,   TJ)   +=  af*(bc[0][0] * cc1 + bc[0][2] * cc2) ;
            Kuu(TI,   TJp1) +=  af*(bc[0][1] * cc2 + bc[0][2] * cc1) ;
            Kuu(TIp1, TJ)   +=  af*(bc[1][0] * cc1 + bc[1][2] * cc2) ;
            Kuu(TIp1, TJp1) +=  af*(bc[1][1] * cc2 + bc[1][2] * cc1) ;
          }

          Kup(TI,   0)  += (af*bb1);
          Kup(TIp1, 0)  += (af*bb2);

          Kpu(0, TI)    += (d2UdJ2*af*bb1);
          Kpu(0, TIp1)  += (d2UdJ2*af*bb2);
        }

        if(finite)
        {
          Flocal2(0) -= ((dUdJ-pres)*eps)*dvol0;
        }
        else
        {
          Flocal2(0) -= (volstrain-pres*eps)*dvol0;
        }
    }//gp

    //printMatrix(Kuu);
    //printMatrix(Kpu);

    Kpp(0,0) = -elemVolOrig*eps*af;

    return 0;
}


int BernsteinElem2DElecMechTri210::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, pbar, volstrain, r2d3 = 2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], streDev[4], cc[4][4], bc[2][4], param[2], bforce[2], force[2];
    double  Idev[4][4], D11[4][4], cctmp[4][4];

    double BULK = matDat[0];
    double eps  = 0.0; if(!SolnData->TRULY_INCOMPRESSIBLE) eps=1.0/BULK;
    double mu   = matDat[1];
    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal1.rows() <= nsize)   Flocal1.resize(nsize);
    Flocal1.setZero();
    if(Flocal2.rows() <= 1)   Flocal2.resize(1);
    Flocal2.setZero();
    if(Kpu.rows() <= 1 || Kpu.cols() <= nsize)   Kpu.resize(1, nsize);
    Kpu.setZero();
    if(Kpp.rows() <= 1 || Kpp.cols() <= 1)   Kpp.resize(1, 1);
    Kpp.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F[0] = computeValuePrev(0, dN_dx) + 1.0;
        F[2] = computeValuePrev(0, dN_dy);
        F[1] = computeValuePrev(1, dN_dx);
        F[3] = computeValuePrev(1, dN_dy) + 1.0;

        //F[0] = computeValue(0, dN_dx) + 1.0;
        //F[2] = computeValue(0, dN_dy);
        //F[1] = computeValue(1, dN_dx);
        //F[3] = computeValue(1, dN_dy) + 1.0;

        volstrain = F[0] + F[3] - 2.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        //cout << " volstrain = " << gp << '\t' << volstrain << '\t' << detF << endl;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)          return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)
          dvol *= F33;


        // evaluate pressure at the quadrature points
        // for this elements its a single value
        pres = SolnData->var3[elenum];

        if(SolnData->TRULY_INCOMPRESSIBLE)
        {
          //stre[0] = 2.0*mu*(F[0]-1.0) ;
          //stre[1] = 2.0*mu*(F[3]-1.0) ;
          //stre[2] = 0.0 ;
          //stre[3] = mu*(F[1]+F[2]);

          stre[0] += pres;
          stre[1] += pres;
          stre[2] += pres;
        }
        else
        {
          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          streDev[0] = stre[0] - pbar;
          streDev[1] = stre[1] - pbar;
          streDev[2] = stre[2] - pbar;

          stre[0] = streDev[0] + pres;
          stre[1] = streDev[1] + pres;
          stre[2] = streDev[2] + pres;
        }

        //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", stre[0], stre[1], stre[2], stre[3]);

        // Calculate Stiffness and Residual
        //==============================================

        xx = yy = 0.0;
        for(ii=0; ii<nlbfU; ii++)
        {
          xx += xNode[ii]*N[ii];
          yy += yNode[ii]*N[ii];
        }

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        //force[0] = analy.computeForce(0, xx, yy, 0.0, tPrev);
        //force[1] = analy.computeForce(1, xx, yy, 0.0, tPrev);

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal1(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3]) ;
          Flocal1(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1]) ;

          Kpu(0, TI)    += bb1;
          Kpu(0, TIp1)  += bb2;
        }

        if(finite)
        {
          Flocal2(0) -= (detF-1.0-pres*eps)*dvol0;
        }
        else
        {
          Flocal2(0) -= (volstrain-pres*eps)*dvol0;
        }
    }//gp

    Kpp(0,0) = -elemVolOrig*eps;

    return 0;
}


int BernsteinElem2DElecMechTri210::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, pbar, volstrain, r2d3 = 2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[4], cc[4][4], param[2], bforce[2], force[2];

    double BULK = matDat[0] ;
    double eps  = 1.0/BULK;
    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;
    double tPrev = myTime.cur - dt;

    //LinearElasticTransient2D  analy(rho0, matDat[0], matDat[1]);

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        xx = yy = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
        }

        F[0] = computeValuePrev(0, dN_dx) + 1.0;
        F[2] = computeValuePrev(0, dN_dy);
        F[1] = computeValuePrev(1, dN_dx);
        F[3] = computeValuePrev(1, dN_dy) + 1.0;

        //F[0] = computeValue(0, dN_dx) + 1.0;
        //F[2] = computeValue(0, dN_dy);
        //F[1] = computeValue(1, dN_dx);
        //F[3] = computeValue(1, dN_dy) + 1.0;

        volstrain = F[0] + F[3] - 2.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        //cout << " volstrain = " << gp << '\t' << volstrain << '\t' << detF << endl;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)          return 1;


        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);
          elemVolCur  += dvol;
        }

        if(finite)
          dvol *= F33;


        // evaluate pressure at the quadrature points
        // for this elements its a single value
        pres = SolnData->var3Prev[elenum];

        fact = pres - (stre[0]+stre[1]+stre[2])/3.0;

        stre[0] += fact;
        stre[1] += fact;
        stre[2] += fact;

        // Calculate Residual
        //==============================================

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        //force[0] = analy.computeForce(0, xx, yy, 0.0, tPrev);
        //force[1] = analy.computeForce(1, xx, yy, 0.0, tPrev);

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3]) ;
          Flocal(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1]) ;
        }
    }//gp

    return 0;
}


int BernsteinElem2DElecMechTri210::calcResidualPressure(VectorXd& Flocal)
{
    int gp, ii;
    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, volstrain, param[2], dvol0, dvol, fact, Jac, r2d3 = 2.0/3.0, dUdJ;
    double  BULK = matDat[0];
    double mu   = matDat[1];
    double LAMBDA = BULK-r2d3*mu;

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() != 1)
      Flocal.resize(1);
    Flocal.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(finite)
        {
          if( matDat[2] == 6 )
          {
            dUdJ   = LAMBDA*log(detF)/detF;
          }
          else
          {
            dUdJ   = 0.5*BULK*(detF-1.0/detF);
          }
          fact = dUdJ;
          //fact = BULK*(detF-1.0);
        }
        else
        {
          volstrain = F[0] + F[3] - 2.0;

          fact = BULK*volstrain;
        }

        Flocal(0) += (fact*dvol0);
    }

    return 0;
}




int BernsteinElem2DElecMechTri210::solveForPressure()
{
    double  volstrain = 0.0, volstrainPrev=0.0;
    double  BULK = matDat[0];

    // explicit scheme
    if(SolnData->tis >= 400)
    {
      if(finite)
      {
        int gp, ii, ll=0, err, isw, count;
        VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

        double  F[4], detF, trF, F33, param[2], dvol0, dvol, Jac, pbar;
        double  cc[4][4], stre[4], dt, fact;

        vector<double>  gausspoints1, gausspoints2, gaussweights;
        getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

        Jbar = 0.0;
        pres = 0.0;
        elemVolOrig = 0.0;
        for(gp=0; gp<nGP; gp++)
        {
            param[0] = gausspoints1[gp];
            param[1] = gausspoints2[gp];

            GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

            dvol0 = gaussweights[gp]*(Jac*thick);
            elemVolOrig  += dvol0;

            F[0] = computeValue(0, dN_dx) + 1.0;
            F[2] = computeValue(0, dN_dy);
            F[1] = computeValue(1, dN_dx);
            F[3] = computeValue(1, dN_dy) + 1.0;

            detF = F[0]*F[3] - F[1]*F[2];

            /*
            // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
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
                fact = pow(Jbar/detF, 1.0/3.0);
                for(ii=0;ii<4;ii++)
                  F[ii] *= fact;
                F33 = fact;
              }
              F33 = 1.0;
            }

            //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP,     &finiteInt, &sss, &isw, &err, &count, NULL);
            count++;
            ll += nivGP;
            pbar = (stre[0]+stre[1]+stre[2])/3.0;

            //pres += pbar*dvol0;
            //Jbar += detF*dvol0;
            */
            pres += (detF-1.0)*dvol0;
        }
        pres = BULK*pres/elemVolOrig;
        //Jbar /= elemVolOrig;
      }
      else
      {
        volstrain  = computeValue(0, dNc_dx);
        volstrain += computeValue(1, dNc_dy);

        pres = BULK*volstrain;
      }
    }
    else // implicit scheme
    {
      volstrain  = computeValueIncr(0, dNc_dx);
      volstrain += computeValueIncr(1, dNc_dy);

      pres += BULK*(volstrain+Rp(0))/elemVolOrig;
    }

    return 0;
}





/*
int BernsteinElem2DElecMechTri210::calcResidualTIC(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], cc[4][4], param[2], bforce[2], force[2];

    double mu   = matDat[1] ;
    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;
    double tPrev = myTime.cur;

    LinearElasticTransientSI2D  analy(rho0, mu);

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();

    Rp = 0.0;
    Jbar = 0.0;
    elemVolOrig=0.0;
    elemVolCur=0.0;
    dNc_dx.setZero();
    dNc_dy.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F[0] = computeValuePrev(0, dN_dx) + 1.0;
        F[2] = computeValuePrev(0, dN_dy);
        F[1] = computeValuePrev(1, dN_dx);
        F[3] = computeValuePrev(1, dN_dy) + 1.0;

        //F[0] = computeValue(0, dN_dx) + 1.0;
        //F[2] = computeValue(0, dN_dy);
        //F[1] = computeValue(1, dN_dx);
        //F[3] = computeValue(1, dN_dy) + 1.0;

        volstrain = F[0] + F[3] - 2.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        //cout << " volstrain = " << gp << '\t' << volstrain << '\t' << detF << endl;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        ////matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        // compute stress and material elasticity tensor
        stre[0] = 2.0*mu*(F[0]-1.0) ;
        stre[1] = 2.0*mu*(F[3]-1.0) ;
        stre[2] = 0.0 ;
        stre[3] = mu*(F[1]+F[2]);

        stre[0] += pres;
        stre[1] += pres;
        stre[2] += pres;

        //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", stre[0], stre[1], stre[2], stre[3]);

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);
          elemVolCur  += dvol;
        }

        if(finite)
          dvol *= F33;

        // Calculate Stiffness and Residual
        //==============================================

        // entries for Kup matrix
        if(finite)
        {
          dNc_dx += dN_dx*dvol;
          dNc_dy += dN_dy*dvol;

          Rp -= (detF-1.0)*dvol0;
        }
        else
        {
          dNc_dx += dN_dx*dvol0;
          dNc_dy += dN_dy*dvol0;

          Rp -= (volstrain)*dvol0;
        }

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        xx = yy = 0.0;
        for(ii=0; ii<nlbfU; ii++)
        {
          xx += xNode[ii]*N[ii];
          yy += yNode[ii]*N[ii];
        }

        force[0] = analy.computeForce(0, xx, yy, 0.0, tPrev);
        force[1] = analy.computeForce(1, xx, yy, 0.0, tPrev);

        //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \n", xx, yy, force[0], force[1]);

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3]) ;
          Flocal(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1]) ;
        }
    }//gp

    return 0;
}
*/


int  BernsteinElem2DElecMechTri210::solveForPressureTIC(VectorXd& matK1, double beta, double dt, VectorXd& Flocal)
{
    double  volstrain = 0.0, volstrainPrev=0.0;

    // explicit scheme
    if(SolnData->tis >= 400)
    {
        //cout << "elemVolOrig=" << elemVolOrig << endl;
        //printVector(forAssyVec);
        int ii, TI, TIp1, ind1, ind2;

        /*
        // compute lhs
        ////////////////////////////////////

        double  Keff=0.0;
        for(ii=0; ii<npElem; ii++)
        {
          ind1 = 2*ii;
          ind2 = 2*nodeNums[ii];

          Keff   +=  dNc_dx[ii]*dNc_dx[ii];
          Keff   +=  dNc_dy[ii]*dNc_dy[ii];
          //cout << matK1[ind2] << '\t' << matK1[ind2+1] << endl;

          // account for the Dirichlet BCs
          //if(forAssyVec[ind1]   != -1)
            //Keff  +=  dNc_dx[ii]*(dNc_dx[ii]/matK1[ind2]);

          //if(forAssyVec[ind1+1] != -1)
            //Keff  +=  dNc_dy[ii]*(dNc_dy[ii]/matK1[ind2+1]);
        }

        double  fact = (elmDat[5]*elemVolOrig/6.0)/beta/dt/dt;
        //cout << " fact = " << fact << endl;
        Keff /= fact;

        // compute rhs
        ////////////////////////////////////

        volstrain  = computeValueIncr(0, dNc_dx);
        volstrain += computeValueIncr(1, dNc_dy);

        volstrainPrev  = computeValuePrev(0, dNc_dx);
        volstrainPrev += computeValuePrev(1, dNc_dy);

        //volstrain += volstrainPrev;

        // compute p_{n+1}
        ////////////////////////////////////

        double presIncr = volstrain/Keff;

        // compute contribution from p_{n+1}
        ////////////////////////////////////

        Flocal.setZero();
        for(ii=0; ii<nlbfU; ++ii)
        {
          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)    -=  dNc_dx[ii]*presIncr;
          Flocal(TIp1)  -=  dNc_dy[ii]*presIncr;
        }

        pres += presIncr;
        */

        //
        // Uzawa algorithm
        ////////////////////////////
        volstrain  = computeValue(0, dNc_dx);
        volstrain += computeValue(1, dNc_dy);

        volstrainPrev  = computeValuePrev(0, dNc_dx);
        volstrainPrev += computeValuePrev(1, dNc_dy);
        //cout << '\t' << matDat[1] << '\t' << volstrain << '\t' << volstrainPrev << endl;
        //volstrainPrev = 0.0;

        //pres += (matDat[1]/elemVolOrig)*(volstrain+volstrainPrev);

        pres += (100.0/elemVolOrig)*(volstrain+volstrainPrev);

        Flocal.setZero();
        for(ii=0; ii<nlbfU; ++ii)
        {
          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)    +=  dNc_dx[ii]*pres;
          Flocal(TIp1)  +=  dNc_dy[ii]*pres;
        }

        //cout << "pres= " << pres << endl;
        //
    }
    else // implicit scheme
    {
      volstrain  = computeValueIncr(0, dNc_dx);
      volstrain += computeValueIncr(1, dNc_dy);
    }

    return 0;
}





int BernsteinElem2DElecMechTri210::calcLoadVector(VectorXd& Flocal)
{
  if( NeumannBCs.size() == 0)
    return 0;


  int side, dir, ii, jj, nn, TI, TIp1, TJ, TJp1, gp, nGP=3;
  int face_node_map[3][3]={{0,3,1},{1,4,2},{2,5,0}};
  double  specVal, xNode[3], yNode[3], xx, yy, param[2], trac[2], normal[2], tang[2], N[3];
  double  Jac, dvol;

  vector<double>  gausspoints1, gaussweights;

  getGaussPoints1D(nGP, gausspoints1, gaussweights);


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

    //cout << side << '\t' << dir << '\t' << specVal << endl;

    for(ii=0; ii<3; ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

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

        for(ii=0;ii<3;ii++)
        {
          TI   = face_node_map[side][ii]*ndof;
          TIp1 = TI+1;

          Flocal(TI)   += N[ii]*trac[0];
          Flocal(TIp1) += N[ii]*trac[1];
        }
    } // for(gp=0; gp<nGP; gp++)
    //cout << " Volume = " << elemVol << endl;
  }

  //printVector(Flocal);

  return 0;
}



int BernsteinElem2DElecMechTri210::calcInternalForces()
{

  return 0;
}



void BernsteinElem2DElecMechTri210::elementContourplot(int vartype, int varindex, int index)
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


void BernsteinElem2DElecMechTri210::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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

      //vector<double>  valtemp;
      //valtemp = vals2project;

      //vals2project[0] = valtemp[0];
      //vals2project[1] = valtemp[1];
      //vals2project[2] = valtemp[2];
      //vals2project[3] = 2.0*( valtemp[3] - 0.25*(valtemp[0]+valtemp[1]) );
      //vals2project[4] = 2.0*( valtemp[4] - 0.25*(valtemp[1]+valtemp[2]) );
      //vals2project[5] = 2.0*( valtemp[5] - 0.25*(valtemp[2]+valtemp[0]) );

      //vals2project[3] = 0.5*valtemp[3] + 0.25*(valtemp[0]+valtemp[1]);
      //vals2project[4] = 0.5*valtemp[4] + 0.25*(valtemp[1]+valtemp[2]);
      //vals2project[5] = 0.5*valtemp[5] + 0.25*(valtemp[2]+valtemp[0]);

    //printVector(vals2project);

    return;
}


void BernsteinElem2DElecMechTri210::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    int   err,  isw,  count,  count1, ll = 0, ii, jj, gp;

    double F[4], Fc[4], detFc, detF, trFc, F33, dvol, Jac, fact, dvol0, pbar;
    double  stre[9], streDev[9], cc[4][4], param[2];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double dt   = myTime.dt;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(sss == 1)  // plane stress
        {
          if(finite)
            F33 = 1.0/sqrt(detF);
          else
            F33 = 3.0 - F[0] - F[3];
        }
        else if(sss == 2)    // plane strain
        {
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        // evaluate pressure at the quadrature points
        pres = SolnData->var3[elenum];

        if(SolnData->TRULY_INCOMPRESSIBLE)
        {
          stre[0] += pres;
          stre[1] += pres;
          stre[2] += pres;
        }
        else
        {
          pbar = (stre[0]+stre[1]+stre[2])/3.0;

          streDev[0] = stre[0] - pbar;
          streDev[1] = stre[1] - pbar;
          streDev[2] = stre[2] - pbar;

          stre[0] = streDev[0] + pres;
          stre[1] = streDev[1] + pres;
          stre[2] = streDev[2] + pres;
        }

        if(varindex < 9)
          outval[gp] = stre[varindex];
        else if(varindex == 9)
          outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*stre[1]*stre[1])/2.0);
        else if(varindex == 10)
          outval[gp] = (stre[0]+stre[4]+stre[8])/3.0;
          //outval[gp] = pres;
        else
          outval[gp] = 0.0;
    }//gp

    return;
}



void BernsteinElem2DElecMechTri210::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void BernsteinElem2DElecMechTri210::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}



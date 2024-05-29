
/* Node ordering for the quadratic quadrilateral element
 4-----7-----3
 |     |     |
 |     |     |
 8-----9-----6
 |     |     |
 |     |     |
 1-----5-----2
*/

#include "MagnetoMech2DQua21.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "BasisFunctionsBernstein.h"
#include "utilitiesmaterial.h"


using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



MagnetoMech2DQua21::MagnetoMech2DQua21()
{
  ELEM_SHAPE = ELEM_SHAPE_QUAD_BERNSTEIN;

  ndim   = 2;
  degree = 2;
  npElem = 9;
  nlbfU  = 9;
  nlbfP  = 4;
  ndof   = 2;
  nsize  = npElem*ndof;
}


MagnetoMech2DQua21::~MagnetoMech2DQua21()
{
}


void MagnetoMech2DQua21::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double MagnetoMech2DQua21::computeVolume(bool init)
{
  double  dvol, Jac, param[2];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

  int   ii, gp;

  double xNode[9], yNode[9], xx, yy;
  for(ii=0;ii<npElem;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
  }

  vector<double>  gausspoints1, gausspoints2, gaussweights;
  getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    elemVol=0.0;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          if(init)
            GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          else
            GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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



int MagnetoMech2DQua21::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != nsize) || (Mlocal.cols() != nsize) )
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = elmDat[5]*computeVolume(true)/9.0;

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

    double xNode[9], yNode[9], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int nGPt=16;
    getGaussPointsQuad(nGPt, gausspoints1, gausspoints2, gaussweights);

    elemVolOrig=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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
    //cout << " elemVolOrig = " << elemVolOrig << endl;
    //printMatrix(Klocal);  printf("\n\n\n");
  }

  return 0;
}



int MagnetoMech2DQua21::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != 4) || (Mlocal.cols() != 4) )
      Mlocal.resize(4, 4);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = computeVolume(true)/4.0;

      for(int ii=0; ii<4; ii++)
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

    double xNode[9], yNode[9], xx, yy;
    for(ii=0;ii<3;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int nGPt=9;
    getGaussPointsQuad(nGPt, gausspoints1, gausspoints2, gaussweights);

    elemVolOrig=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, 1, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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

        for(ii=0;ii<4;ii++)
        {
          bb1 = dvol0*N[ii];

          for(jj=0; jj<4; jj++)
          {
            Mlocal(ii, jj)  += bb1*N[jj] ;
          }
        }
    } //gp
    //cout << " elemVolOrig = " << elemVolOrig << endl;
    //printMatrix(Klocal);  printf("\n\n\n");
  }

  return 0;
}


double  MagnetoMech2DQua21::calcCriticalTimeStep(bool flag)
{
    int  ii;

    double xNode[9], yNode[9];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosNew[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosNew[nodeNums[ii]][1];
    }

    double  xd, yd, charlen=1.0e10;
    int  edge_nodes[4][2] = { {0,1}, {1,2}, {2,3}, {3,0}};

    for(ii=0; ii<4; ii++)
    {
      xd = xNode[edge_nodes[ii][0]] - xNode[edge_nodes[ii][1]];
      yd = yNode[edge_nodes[ii][0]] - yNode[edge_nodes[ii][1]]; 

      charlen = min(charlen, xd*xd+yd*yd);
    }
    charlen = 0.5*sqrt(charlen);

    double K    = matDat[0] ;
    double mu   = matDat[1] ;
    double rho  = elmDat[5] ;
    double wave_speed;
    
    //if(SolnData->TRULY_INCOMPRESSIBLE)
      wave_speed = sqrt(mu/rho);
    //else
      //wave_speed = sqrt((K+4.0*mu/3.0)/rho);

    return  (charlen/wave_speed);
}







int MagnetoMech2DQua21::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], Fe[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, volstrain, pbar;
    double  bb1, bb2, bb3, bb4, bb5, dUdJ, d2UdJ2, Rp;
    double  stre[4], streDev[4], cc[4][4], param[2], bforce[2], force[2], Np[4];

    double Jg = 1.0 + elmDat[10]*timeFunction[0].prop ;

    double BULK = matDat[0];
    double eps  = 0.0; if(!SolnData->TRULY_INCOMPRESSIBLE) eps=1.0/BULK;
    double mu   = matDat[1];
    double LAMBDA = BULK-2.0*mu/3.0;
    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;
    double tPrev = myTime.cur - dt;

    double xNode[9], yNode[9], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal1.rows() < nsize)   Flocal1.resize(nsize);
    Flocal1.setZero();
    if(Flocal2.rows() < 4)   Flocal2.resize(4);
    Flocal2.setZero();
    if(Kpu.rows() < 4 || Kpu.cols() < nsize)   Kpu.resize(4, nsize);
    Kpu.setZero();
    if(Kup.rows() < nsize || Kup.cols() < 4)   Kup.resize(nsize, 4);
    Kup.setZero();
    if(Kpp.rows() < 4 || Kpp.cols() < 4)   Kpp.resize(4, 4);
    Kpp.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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

        fact = pow(Jg, -1.0/3.0);

        Fe[0] = F[0]*fact;
        Fe[1] = F[1]*fact;
        Fe[2] = F[2]*fact;
        Fe[3] = F[3]*fact;

        F33 = fact;

        //matlib2d_(matDat, Fe, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)
          dvol *= F33;

        // basis functions for pressure dof
        BernsteinBasisFunsQuad(1, param[0], param[1], Np);

        // evaluate pressure at the quadrature points
        pres = 0.0;
        for(ii=0; ii<4; ii++)
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);

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
        force[0] = 0.0;
        force[1] = 0.0;

        //force[0] = analy.computeForce(0, xx, yy, 0.0, tPrev);
        //force[1] = analy.computeForce(1, xx, yy, 0.0, tPrev);

        //double  xF = 0.5, yF = 0.5;
        //force[0] = -(2.0*tPrev/0.0005)*exp(-(tPrev-0.1)*(tPrev-0.1)/0.0005)*100.0*exp(-((xx-xF)*(xx-xF)+(yy-yF)*(yy-yF))/0.0005);
        //force[0] = 100.0*(2.0*tPrev/0.005)*exp(-(tPrev-0.2)*(tPrev-0.2)/0.005)*exp(-((xx-xF)*(xx-xF)+(yy-yF)*(yy-yF))/0.005);
        //force[0] = -exp(-((xx-xF)*(xx-xF)+(yy-yF)*(yy-yF))/0.0005);
        //force[1] = 0.0;
        //force[1] = force[0];

        //cout << " force = " << force[0] << '\t' << force[1] << endl;

        if(!finite)
        {
          d2UdJ2 = 1.0;

          Rp = volstrain-pres*eps;
        }

        if(SolnData->TRULY_INCOMPRESSIBLE)
        {
          stre[0] += pres;
          stre[1] += pres;
          stre[2] += pres;

          if(finite)
          {
            d2UdJ2 = 1.0;

            Rp = detF/Jg-1.0-pres*eps;
          }
        }
        else
        {
          fact = pres - (stre[0]+stre[1]+stre[2])/3.0;

          stre[0] += fact;
          stre[1] += fact;
          stre[2] += fact;

          if(finite)
          {
            if( (int) matDat[2] == 6 )
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

            Rp = (dUdJ/Jg-pres)*eps;
          }
        }

        for(ii=0;ii<4;ii++)
        {
          stre[ii] /= Jg;
        }


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

          bb1 /= Jg;
          bb2 /= Jg;

          for(jj=0; jj<4; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
          }
        }

        for(ii=0; ii<4; ii++)
        {
          bb4 = dvol*Np[ii]*d2UdJ2;
          bb5 = dvol0*Np[ii];

          Flocal2(ii) -= bb5*Rp;

          bb4 /= Jg;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Kpu(ii, TJ)    += bb4*dN_dx[jj];
            Kpu(ii, TJp1)  += bb4*dN_dy[jj];
          }

          bb5 *= eps;
          for(jj=0; jj<4; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

    return 0;
}



int MagnetoMech2DQua21::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ii, jj, ll, gp, TI, TIp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, pbar, volstrain, r2d3 = 2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[4], cc[4][4], param[2], bforce[2], force[2], Np[4];

    double BULK = matDat[0] ;
    double eps  = 1.0/BULK;
    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;
    double tPrev = myTime.cur - dt;

    //LinearElasticTransient2D  analy(rho0, matDat[0], matDat[1], elmDat[5]);

    double xNode[9], yNode[9], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

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

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)
          dvol *= F33;

        // basis functions for pressure dof
        LagrangeBasisFunsQuad(1, param[0], param[1], Np);

        // evaluate pressure at the quadrature points
        pres = 0.0;
        for(ii=0; ii<4; ii++)
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);

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


int MagnetoMech2DQua21::calcResidualPressure(VectorXd& Flocal)
{
    int gp, ii;
    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, volstrain, param[2], dvol0, dvol;
    double  fact, Jac, r2d3 = 2.0/3.0, dUdJ, Np[4];
    double  BULK = matDat[0];
    double mu   = matDat[1];
    double LAMBDA = BULK-r2d3*mu;

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() != 4)
      Flocal.resize(4);
    Flocal.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);

        if(finite)
        {
          F[0] = computeValue(0, dN_dx) + 1.0;
          F[2] = computeValue(0, dN_dy);
          F[1] = computeValue(1, dN_dx);
          F[3] = computeValue(1, dN_dy) + 1.0;

          detF = F[0]*F[3] - F[1]*F[2];

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
          volstrain = computeValue(0, dN_dx) + computeValue(1, dN_dy);

          fact = BULK*volstrain;
        }

        // basis functions for pressure dof
        LagrangeBasisFunsQuad(1, param[0], param[1], Np);

        fact *= dvol0;
        for(ii=0; ii<4; ii++)
        {
          Flocal(ii) += (fact*Np[ii]);
        }
    }

    return 0;
}


int MagnetoMech2DQua21::solveForPressure()
{
    cout << "'MagnetoMech2DQua21::solveForPressure()' ... not implemented yet " << endl;

    return 0;
}



int MagnetoMech2DQua21::calcLoadVector(VectorXd& Flocal)
{
  if( NeumannBCs.size() == 0)
    return 0;

  int side, dir, ii, jj, nn, TI, TIp1, TJ, TJp1, gp, nGP=3;
  int face_node_map[4][3]={{0,1,4},{1,2,5},{2,3,6},{3,0,7}};
  double  specVal, xNode[9], yNode[9], xx, yy, param[2], trac[2], normal[2], tang[2], N[3], dNdx[3], dNdy[3];
  double  Jac, dvol, curvature;

  vector<double>  gausspoints1, gaussweights;
  getGaussPoints1D(nGP, gausspoints1, gaussweights);

  // side 0 - 1-5-2
  // side 1 - 2-6-3
  // side 2 - 3-7-4
  // side 2 - 4-8-1

  // loop over all the sides and compute force vector due to applied traction
  for(nn=0; nn<NeumannBCs.size(); nn++)
  {
    side    = (int) NeumannBCs[nn][0];                      // edge/face number
    dir     = (int) NeumannBCs[nn][1];                     // DOF value
    specVal = NeumannBCs[nn][2]; // specified value
    specVal *= timeFunction[0].prop; // load factor

    for(ii=0; ii<3; ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[face_node_map[side][ii]]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[face_node_map[side][ii]]][1];
    }

    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
      param[0] = gausspoints1[gp];

      BernsteinBasisFunsEdge2D(degree, param, xNode, yNode, N, dNdx, dNdy, normal, curvature, Jac);

      dvol = gaussweights[gp]*(Jac*thick);
      elemVol += dvol;

      //cout << " side = " << side << '\t' << dir << '\t' << normal[0] << '\t' << normal[1] << '\t' << dvol << '\t' << curvature << endl;

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
        else if(dir == 1)
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



int MagnetoMech2DQua21::calcInternalForces()
{
  return 0;
}




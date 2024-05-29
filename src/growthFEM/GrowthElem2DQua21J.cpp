
/* Node ordering for the quadratic quadrilateral element
 4-----7-----3
 |     |     |
 |     |     |
 8-----9-----6
 |     |     |
 |     |     |
 1-----5-----2
*/

#include "GrowthElem2DQua21J.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "BasisFunctionsBernstein.h"
#include "growthfunctions.h"
#include "utilitiesmaterialHyperelastic.h"


using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



GrowthElem2DQua21J::GrowthElem2DQua21J()
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


GrowthElem2DQua21J::~GrowthElem2DQua21J()
{
}


void GrowthElem2DQua21J::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double GrowthElem2DQua21J::computeVolume(bool init)
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



int GrowthElem2DQua21J::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
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



int GrowthElem2DQua21J::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
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


double  GrowthElem2DQua21J::calcCriticalTimeStep(bool flag)
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




int GrowthElem2DQua21J::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, detFn, fact, fact1, dvol, dvol0, Jac, volstrain, volstrainnew, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Rp, Jhat, thetahat;
    double  param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2], Np[4];

    MatrixXd  Fg(3,3), FgInv(3,3), F(3,3), Fe(3,3), Fn(3,3);
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    Fg.setZero();
    FgInv.setZero();

    double  timeFactor = timeFunction[0].prop, Jg = 1.0;
    double  growthFactor = elmDat[10]*timeFunction[0].prop ;

    double mu   = matDat[0];      // shear modulus
    double eps  = matDat[1];      // inverse of bulk modulus
    double BULK = 1.0/eps;
    int  Utype  = int(matDat[2]); // volumetric energy function

    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;
    double tCur = myTime.cur;

    double xNode[9], yNode[9], geom[3];
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
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();
    if(Kup.rows() < nsize || Kup.cols() < 4)   Kup.resize(nsize, 4);
    Kup.setZero();
    if(Kpu.rows() < 4 || Kpu.cols() < nsize)   Kpu.resize(4, nsize);
    Kpu.setZero();
    if(Kpp.rows() < 4 || Kpp.cols() < 4)   Kpp.resize(4, 4);
    Kpp.setZero();

    MatrixXd  Cmat(9,9), Gc(3,9);
    VectorXd  stre(9);
    Gc.setZero();
    stre.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    //cout << elmType << '\t' <<  matType << '\t' <<  secType << endl;
    //printVector(MatlData->matData);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        geom[0] = geom[1] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += N[ii]*xNode[ii];
          geom[1] += N[ii]*yNode[ii];
        }

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;


        Fn(0,0) = computeValuePrev(0, dN_dx) + 1.0;            // eps_xx
        Fn(0,1) = computeValuePrev(0, dN_dy);                  // eps_xy
        Fn(1,0) = computeValuePrev(1, dN_dx);                  // eps_yx
        Fn(1,1) = computeValuePrev(1, dN_dy) + 1.0;            // eps_yy

        detFn = Fn(0,0)*Fn(1,1) - Fn(0,1)*Fn(1,0);

        compute_constants_volumetric_functions(finite, Utype, detFn, BULK, Jhat, thetahat);


        F(0,0) = computeValueCur(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dN_dy);                  // eps_xy
        F(1,0) = computeValueCur(1, dN_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        volstrain = F(0,0) + F(1,1) - 2.0;

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);

        if(detF < 0.0)   return 1;

        // ADJUST F(2,2) for 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        if(sss == 1)  // plane stress
        {
          if(finite)
            F(2,2) = 1.0/sqrt(detF);
          else
            F(2,2) = 3.0 - F(0,0) - F(1,1);
        }
        else if(sss == 2)    // plane strain
        {
          F(2,2) = 1.0;
        }

        //growthfunction_beam_singlelayer(timeFactor,  geom,  Fg);
        //growthfunction_beam_multilayer(timeFactor, geom,  Fg);
        //growthfunction_isotropic(1, ndim,  growthFactor,  Fg);
        //FgInv = Fg.inverse();
        //Jg = Fg.determinant();

        // elastic part of the deformation gradient
        //Fe = F*FgInv;
        Fe = F;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        //if(finite)  dvol *= F(2,2);

        // evaluate pressure at the quadrature points
        BernsteinBasisFunsQuad(1, param[0], param[1], Np);
        volstrainnew = 0.0;
        for(ii=0; ii<4; ii++)
          volstrainnew += (Np[ii]*SolnData->var2[nodeNums[ii]]);

        MatlData->computeStressAndTangent(true, sss, Fn, Fe, volstrainnew, stre, Cmat);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;


        // Calculate Stiffness and Residual
        //==============================================

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          Gc(0,0) = bb1 * Cmat(0,0) + bb2 * Cmat(3,0);
          Gc(0,1) = bb1 * Cmat(0,1) + bb2 * Cmat(3,1);
          Gc(0,3) = bb1 * Cmat(0,3) + bb2 * Cmat(3,3);
          Gc(0,4) = bb1 * Cmat(0,4) + bb2 * Cmat(3,4);

          Gc(1,0) = bb1 * Cmat(1,0) + bb2 * Cmat(4,0);
          Gc(1,1) = bb1 * Cmat(1,1) + bb2 * Cmat(4,1);
          Gc(1,3) = bb1 * Cmat(1,3) + bb2 * Cmat(4,3);
          Gc(1,4) = bb1 * Cmat(1,4) + bb2 * Cmat(4,4);

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[1] + bb2*stre[4];

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

            Kuu(TI,   TJ)    +=  af*(Gc(0,0) * cc1 + Gc(0,3) * cc2) ;
            Kuu(TI,   TJp1)  +=  af*(Gc(0,1) * cc1 + Gc(0,4) * cc2) ;
            Kuu(TIp1, TJ)    +=  af*(Gc(1,0) * cc1 + Gc(1,3) * cc2) ;
            Kuu(TIp1, TJp1)  +=  af*(Gc(1,1) * cc1 + Gc(1,4) * cc2) ;
          }

          bb1 *= BULK;
          bb2 *= BULK;

          for(jj=0; jj<4; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
          }
        }

        if(finite)
        {
          //Rp = detF/Jg - 1.0 - pres*eps;

          Rp = detF - Jg*(Jhat+thetahat*pres);
        }
        else
        {
          Rp = volstrain-volstrainnew;
        }

        for(ii=0; ii<4; ii++)
        {
          //bb4 = dvol*BULK*Np[ii]*af;
          //bb5 = dvol0*BULK*Np[ii];

          bb4 = dvol*Np[ii]*af;
          bb5 = dvol0*Np[ii];

          Flocal2(ii) -= bb5*Rp;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Kpu(ii, TJ)    += (bb4*dN_dx[jj]);
            Kpu(ii, TJp1)  += (bb4*dN_dy[jj]);
          }

          for(jj=0; jj<4; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

    //printMatrix(Kuu);
    //printMatrix(Kpu);
    //printMatrix(Kup);
    //printMatrix(Kpp);
    //printVector(Flocal1);
    //printVector(Flocal2);

    return 0;
}






int GrowthElem2DQua21J::calcForceVectorGrowthModel(VectorXd& Flocal1, VectorXd& Flocal2)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2], Np[4];

    double  gp1 = 1.0 + elmDat[10]*timeFunction[0].prop ;
    double  Jg  = gp1*gp1;
    MatrixXd  FgInv(3,3), F(3,3), Fn(3,3), Fg(3,3), Fe(3,3);
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    Fg.setZero();
    FgInv.setZero();

    FgInv(0,0) = 1.0/gp1;    FgInv(0,1) = 0.0;        FgInv(0,2) = 0.0;
    FgInv(1,0) = 0.0;        FgInv(1,1) = 1.0/gp1;    FgInv(1,2) = 0.0;
    FgInv(2,0) = 0.0;        FgInv(2,1) = 0.0;        FgInv(2,2) = 1.0;

    double  FgH = -2.0/gp1;

    double rho  = elmDat[5] ;
    double rho0 = rho ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

    double xNode[4], yNode[4], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Flocal1.rows() != nsize)
    {
      Flocal1.resize(nsize);
    }
    Flocal1.setZero();

    if(Flocal2.rows() != 4)
    {
      Flocal2.resize(nsize);
    }
    Flocal2.setZero();


    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    MatrixXd  Cmat(9,9), Gc(3,9);
    VectorXd  stre(9);
    Gc.setZero();
    stre.setZero();


    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
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

        dvol  = dvol0;

        F(0,0) = computeValueCur(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dN_dy);                  // eps_xy
        F(1,0) = computeValueCur(1, dN_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        // ADJUST F(2,2) for 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
        if(sss == 1)  // plane stress
        {
          if(finite)
            F(2,2) = 1.0/sqrt(detF);
          else
            F(2,2) = 3.0 - F(0,0) - F(1,1);
        }
        else if(sss == 2)    // plane strain
        {
          F(2,2) = 1.0;
        }

        // elastic part of the deformation gradient
        Fe = F*FgInv;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(thick*Jac);
        }

        // evaluate pressure at the quadrature points
        BernsteinBasisFunsQuad(1, param[0], param[1], Np);
        pres = 0.0;
        for(ii=0; ii<4; ii++)
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);

        MatlData->computeStressAndTangent(true, sss, Fn, Fe, pres, stre, Cmat);
        count++;
        ll += nivGP;

        Cmat /=  Jg;
        stre /=  Jg;


        if(finite)  dvol *= F(2,2);

        // multiply with the factor corresponding to the derivative of F wrt to growth factor
        stre *= FgH;

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[1] + bb2*stre[4];

          Flocal1(TI)   += (sig[0]) ;
          Flocal1(TIp1) += (sig[1]) ;
        }


        Rp(0) = (-detF/Jg/Jg)*(2.0*gp1);

        for(ii=0; ii<4; ii++)
        {
          bb5 = dvol0*Np[ii];

          Flocal2(ii) -= bb5*Rp(0);
        }

    } //gp

    //printVector(Flocal);  printf("\n\n\n");

    return 0;
}




int GrowthElem2DQua21J::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
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



int GrowthElem2DQua21J::calcResidual(VectorXd& Flocal)
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


int GrowthElem2DQua21J::calcResidualPressure(VectorXd& Flocal)
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


int GrowthElem2DQua21J::solveForPressure()
{
    cout << "'GrowthElem2DQua21J::solveForPressure()' ... not implemented yet " << endl;

    return 0;
}



int GrowthElem2DQua21J::calcLoadVector(VectorXd& Flocal)
{
  if( NeumannBCs.size() == 0)
    return 0;

  int side, dir, ii, jj, nn, TI, TIp1, TJ, TJp1, gp, nGP=3;
  int face_node_map[4][3]={{0,4,1},{1,5,2},{2,6,3},{3,7,0}};
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

      BernsteinBasisFunsEdge2D(degree, param, xNode, yNode, N, dNdx, dNdy, normal, curvature, Jac);

      dvol = gaussweights[gp]*(Jac*thick);
      elemVol += dvol;

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



int GrowthElem2DQua21J::calcInternalForces()
{
  return 0;
}



void GrowthElem2DQua21J::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'GrowthElem2DQua21J::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
    {
      vals2project[0] += outval[ii];
    }

    vals2project[0] /= nGP;

    return;
}


void GrowthElem2DQua21J::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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
    stressrecovery_extrapolate_Quadrilateral(degree, outval, &vals2project[0]);

    //printVector(vals2project);

    return;
}


void GrowthElem2DQua21J::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    int   err,  isw,  count,  count1, ll = 0, ii, jj, gp;

    double  detF, dvol, Jac, fact, dvol0, param[2], Np[4], volstrainnew;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), stre(9);

    double  timeFact = elmDat[10]*timeFunction[0].prop, Jg ;
    MatrixXd  Fg(3,3), FgInv(3,3), F(3,3), Fe(3,3);
    F.setZero();
    double dt   = myTime.dt;

    double xNode[9], yNode[9], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        geom[0] = geom[1] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += N[ii]*xNode[ii];
          geom[1] += N[ii]*yNode[ii];
        }

        F(0,0) = computeValue(0, dN_dx) + 1.0;
        F(0,1) = computeValue(0, dN_dy);
        F(1,0) = computeValue(1, dN_dx);
        F(1,1) = computeValue(1, dN_dy) + 1.0;

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        if(sss == 1)  // plane stress
        {
          if(finite)
            F(2,2) = 1.0/sqrt(detF);
          else
            F(2,2) = 3.0 - F(0,0) - F(1,1);
        }
        else if(sss == 2)    // plane strain
        {
          F(2,2) = 1.0;
        }

        growthfunction_beam_multilayer(timeFact, geom,  Fg);
        FgInv = Fg.inverse();
        Jg = Fg.determinant();

        //Fe = F*FgInv;
        Fe = F;

        // basis functions for pressure dof
        BernsteinBasisFunsQuad(1, param[0], param[1], Np);

        // evaluate pressure at the quadrature points
        volstrainnew = 0.0;
        for(ii=0; ii<4; ii++)
          volstrainnew += (Np[ii]*SolnData->var2[nodeNums[ii]]);

        MatlData->computeMechanicalStress(Fe, volstrainnew, stre);
        count++;
        ll += nivGP;

        if(varindex < 9)
           outval[gp] = stre[varindex];
        else if(varindex == 9)
           outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6* stre[1]*stre[1])/2);
        else if(varindex == 10)
          outval[gp] = (stre[0]+stre[4]+stre[8])/3.0;

        //outval[gp] = detF;

    }//gp

    return;
}



void GrowthElem2DQua21J::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void GrowthElem2DQua21J::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}





int  GrowthElem2DQua21J::calcError(int index)
{
    // compute error
    ///////////////////////////////////////////////////////////

  int   err,  isw,  count,  count1, ll, ii, jj, kk, gp;

  double  param[2], Jac, dvol0, dvol, rad, theta, val, fact;
  double  dispExact[2], dispNum[2], streDev[4], pbar;
  double  F[4], detF, F33, streExac[4], streNum[4], cc[4][4], Np[4];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

  double xNode[9], yNode[9], geom[3];
  for(ii=0;ii<nlbfU;ii++)
  {
    xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
    yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
  }

  //ThickCylinder  analy(sss, matDat[1], matDat[2]);
  //ElasticityFiniteStrain  analy(matDat[1], matDat[2]);
  //ThickCylinder  analy(sss, matDat[1], matDat[2]);
  double mu   = matDat[0] ;
  double rho0 = elmDat[5] ;
  double dt   = myTime.dt;
  double tCur = myTime.cur;

  double  lambda10=1.0, lambda30=1.0, beamlength=1.0, lambda1, lambda3, lambda11, alpha;
  lambda10 = 1.0;
  lambda11 = PI;
  lambda30 = 1.0;
  alpha = pi;


  vector<double>  gausspoints1, gausspoints2, gaussweights;
  getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

  elemError = 0.0;
  if(index == 0) // L2 norm
  {
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp] * Jac;

          geom[0] = geom[1] = 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            geom[0] += xNode[ii]*N[ii];
            geom[1] += yNode[ii]*N[ii];
          }

          dispNum[0] = computeValue(0, N);
          dispNum[1] = computeValue(1, N);


          lambda1 = lambda10 + lambda11*geom[1]/beamlength;
          lambda3 = lambda30;

          rad   = lambda30*geom[1] + beamlength*lambda10/alpha;

          dispExact[0] = rad*sin(geom[0]*alpha/beamlength) - geom[0];
          dispExact[1] = ( rad*cos(geom[0]*alpha/beamlength) - beamlength*lambda10/alpha ) - geom[1];

          //cout << dispExact[0] << '\t' << dispNum[0] << endl;
          //cout << dispExact[1] << '\t' << dispNum[1] << endl;

          dispExact[0] -= dispNum[0];
          dispExact[1] -= dispNum[1];

          fact = dispExact[0]*dispExact[0] + dispExact[1]*dispExact[1];

          elemError += ( fact * dvol );
    }//gp
  }
  else if(index == 1)
  {
    /*
    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);

          xx = yy= 0.0;
          for(ii=0;ii<nlbf;ii++)
          {
            xx += N[ii]*xNode[ii];
            yy += N[ii]*yNode[ii];
          }

          F[0] = computeValueCur(0, dN_dx) + 1.0;
          F[2] = computeValueCur(0, dN_dy);
          F[1] = computeValueCur(1, dN_dx);
          F[3] = computeValueCur(1, dN_dy) + 1.0;

          detF =  F[0]*F[3] - F[1]*F[2];

          // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

          if(sss == 1)  // plane stress
          {
            if(finite)
              F33 = 1.0/sqrt(detF);
            else
              F33 = 3.0 - F[0] - F[3];
          }
          else if(sss == 2)    // plane strain
            F33 = 1.0;

          matlib2d_(matDat, F, &F33, streNum, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
          count++;
          ll += nivGP;

          pbar = (streNum[0]+streNum[1]+streNum[2])/3.0;

          streDev[0] = streNum[0] - pbar;
          streDev[1] = streNum[1] - pbar;
          streDev[2] = streNum[2] - pbar;
          streDev[3] = streNum[3];

          streNum[0] = streDev[0] + pres;
          streNum[1] = streDev[1] + pres;
          streNum[2] = streDev[2] + pres;

          rad   = sqrt(xx*xx + yy*yy) ;
          theta = atan2(yy, xx);

          for(ii=0; ii<4; ii++)
            streExac[ii] -= streNum[ii];

          //fact = streExac[0]*streExac[0] + streExac[1]*streExac[1] + streExac[3]*streExac[3];
          fact = streExac[0]*streExac[0] + streExac[1]*streExac[1] + streExac[2]*streExac[2] + streExac[3]*streExac[3];

          elemError += ( fact * dvol );
    }//gp
    */
  }
  else if(index == 2)
  {
    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];

          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);

          geom[0] = geom[1] = 0.0;
          for(ii=0;ii<nlbfU;ii++)
          {
            geom[0] += xNode[ii]*N[ii];
            geom[1] += yNode[ii]*N[ii];
          }

          // evaluate pressure at the quadrature points
          BernsteinBasisFunsQuad(1, param[0], param[1], Np);
          pres = 0.0;
          for(ii=0; ii<4; ii++)
            pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);

          pbar = 0.0;

          val = pres - pbar;
          fact = val*val;

          elemError += ( fact * dvol );
    }//gp
  }

  return 0;
}


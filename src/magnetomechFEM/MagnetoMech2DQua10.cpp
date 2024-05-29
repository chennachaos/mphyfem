
#include "MagnetoMech2DQua10.h"
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



MagnetoMech2DQua10::MagnetoMech2DQua10()
{
  ELEM_SHAPE = ELEM_SHAPE_QUAD_BERNSTEIN;

  ndim   = 2;
  degree = 1;
  npElem = 4;
  nlbfU  = 4;
  nlbfP  = 1;
  ndof   = 2;
  nsize  = npElem*ndof;

  Kup.resize(nsize, 1);
  Kpu.resize(1, nsize);
  Kpp.resize(1, 1);
  Rp.resize(1);

  presDOF.resize(1); presDOF.setZero();
  presDOFprev  = presDOF;
  presDOFprev2 = presDOF;
}


MagnetoMech2DQua10::~MagnetoMech2DQua10()
{
}


void MagnetoMech2DQua10::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double MagnetoMech2DQua10::computeVolume(bool configflag)
{
  double  dvol, Jac, param[2];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

  int   ii, gp;

  double xNode[4], yNode[4], xx, yy;
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

          if(configflag)
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



int MagnetoMech2DQua10::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
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

    double xNode[4], yNode[4], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int nGPt=7;
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
    cout << " elemVolOrig = " << elemVolOrig << endl;
    //printMatrix(Klocal);  printf("\n\n\n");
  }

  return 0;
}


int MagnetoMech2DQua10::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != 1) || (Mlocal.cols() != 1) )
      Mlocal.resize(1, 1);

    Mlocal(0,0) = computeVolume(true);

    return 0;
}


double  MagnetoMech2DQua10::calcCriticalTimeStep(bool flag)
{
    int  ii;

    double xNode[4], yNode[4];
    for(ii=0;ii<npElem;ii++)
    {
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



int MagnetoMech2DQua10::calcStiffnessAndResidual(MatrixXd& Kuu, VectorXd& FlocalU, bool firstIter)
{
    int   err = 0,  ii, jj, kk, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, detFn, Jne, fact, fact1, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Jhat, thetahat;
    double  param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    MatrixXd  F(3,3), Fn(3,3);
    F.setZero();
    Fn.setZero();

    int  Utype  = int(matDat[0]); // volumetric energy function
    double eps  = matDat[1];      // inverse of bulk modulus
    double BULK = 1.0/eps;

    if(Utype == -1)
    {
        cerr << "This element does not support the truly incompressible case" << endl;
    }

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

    double xNode[4], yNode[4], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();

    Kup.setZero();
    Kpu.setZero();
    Kpp.setZero();
    Rp.setZero();

    MatrixXd  Cmat(9,9), Gc(3,9);
    VectorXd  stre(9);
    Gc.setZero();
    stre.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    //cout << elmType << '\t' <<  matType << '\t' <<  secType << endl;
    //printVector(MatlData->matData);

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

        compute_constants_volumetric_functions(finite, Utype, Jne, BULK, Jhat, thetahat);

        F(0,0) = computeValueCur(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dN_dy);                  // eps_xy
        F(1,0) = computeValueCur(1, dN_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        volstrain = F(0,0) + F(1,1) - 2.0;

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

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        //if(finite)  dvol *= F(2,2);

        // evaluate pressure at the quadrature points
        pres = presDOF[0];

        MatlData->computeStressAndTangent(true, sss, Fn, F, pres, stre, Cmat);
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

          FlocalU(TI)   += (bb5*force[0] - sig[0]) ;
          FlocalU(TIp1) += (bb5*force[1] - sig[1]) ;

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

          Kup(TI,   0)  += (bb1);
          Kup(TIp1, 0)  += (bb2);
        }

        if(finite)
        {
          Rp(0) += (detF - (Jhat+thetahat*pres)) *dvol0;
        }
        else
        {
          Rp(0) += (volstrain-pres*eps) *dvol0;
        }

        bb4 = dvol*af;

        for(jj=0; jj<nlbfU; jj++)
        {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Kpu(0, TJ)    += (bb4*dN_dx[jj]);
            Kpu(0, TJp1)  += (bb4*dN_dy[jj]);
        }

        Kpp(0, 0)    -= (dvol0*thetahat);
    }//gp

    // add contribution from the condensed matrix

    Kuu -=  ( Kup*(Kpu/Kpp(0,0)) );
    FlocalU +=  ( Kup*(Rp/Kpp(0,0)) );

    return 0;
}



int MagnetoMech2DQua10::solveForPressure()
{
    double  volstrain = 0.0, volstrainPrev=0.0;
    double  BULK = 1.0/matDat[1];

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
        getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

        Jbar = 0.0;
        pres = 0.0;
        elemVolOrig = 0.0;
        for(gp=0; gp<nGP; gp++)
        {
            param[0] = gausspoints1[gp];
            param[1] = gausspoints2[gp];

            GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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
      VectorXd  dispIncr(nsize);
      int ii, kk, ind1, ind2;
      for(ii=0; ii<npElem; ii++)
      {
        ind1 = ii*2;
        ind2 = nodeNums[ii]*2;

        for(kk=0; kk<2; kk++)
        {
          dispIncr(ind1+kk)  =  SolnData->var1Incr[ind2+kk];
        }
      }

      VectorXd  vtmp = Kpu*dispIncr;

      presDOF[0] += (-Rp(0) - vtmp(0))/Kpp(0,0);
    }

    //cout <<  "pres = " <<  pres <<  endl;

    return 0;
}



int MagnetoMech2DQua10::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
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

    double xNode[4], yNode[4], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

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

        //matlib2d_(matDat, F, &F33, stre, D11[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
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


int MagnetoMech2DQua10::calcResidual(VectorXd& Flocal)
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

    double xNode[4], yNode[4], xx, yy;
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


int MagnetoMech2DQua10::calcResidualPressure(VectorXd& Flocal)
{
    int gp, ii;
    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, volstrain, param[2], dvol0, dvol, fact, Jac, r2d3 = 2.0/3.0, dUdJ;
    double  BULK = matDat[0];
    double mu   = matDat[1];
    double LAMBDA = BULK-r2d3*mu;

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() != 1)
      Flocal.resize(1);
    Flocal.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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






int MagnetoMech2DQua10::calcLoadVector(VectorXd& Flocal)
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



int MagnetoMech2DQua10::calcInternalForces()
{

  return 0;
}



void MagnetoMech2DQua10::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'MagnetoMech2DQua10::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
      vals2project[0] += outval[ii];

    vals2project[0] /= nGP;

    return;
}


void MagnetoMech2DQua10::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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


void MagnetoMech2DQua10::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    int   err = 0, ii, jj, gp;

    double  detF, dvol, Jac, fact, dvol0, pbar, param[3];
    MatrixXd  Cmat(9,9), Bmat(3,9), Amat(3,3);
    VectorXd  stre(9),  elecDisp(3),  elecField(3);
    stre.setZero();

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double dt   = myTime.dt;

    double xNode[4], yNode[4], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    MatrixXd  F(3,3), Fn(3,3);
    F.setZero();
    Fn.setZero();

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

        pres = presDOF[0];
        MatlData->computeMechanicalStress(F, pres, stre);

        if(varindex < 9)
          outval[gp] = stre[varindex];
        else if(varindex == 9)
          outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*stre[1]*stre[1])/2.0);
        else if(varindex == 10)
          outval[gp] = pres;
        else
          //outval[gp] = eivals[0].real();
          outval[gp] = 0.0;
    }//gp

    return;
}



void MagnetoMech2DQua10::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void MagnetoMech2DQua10::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}




int  MagnetoMech2DQua10::calcError(int index)
{
    // compute error
    ///////////////////////////////////////////////////////////

  int   err,  isw,  count,  count1, ll, ii, jj, kk, gp;

  double  param[2], Jac, dvol0, dvol, xx, yy, rad, theta, val, fact;
  double  dispExact[2], dispNum[2], streDev[4], pbar;
  double  detF, streExac[4], streNum[4], cc[4][4], Np[4];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), stre(9);

  double xNode[4], yNode[4], geom[3];
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

    MatrixXd  F(3,3), Fe(3,3), Fg(3,3), FgInv(3,3);
    double  timeFactor = timeFunction[0].prop, Jg;
    double  growthFactor = elmDat[10]*timeFunction[0].prop ;


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

          //cout << dispExact[0] << '\t' << dispNum[0] << '\t' << dispExact[1] << '\t' << dispNum[1] << endl;

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

        F(0,0) = computeValueCur(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dN_dy);                  // eps_xy
        F(1,0) = computeValueCur(1, dN_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        growthfunction_beam_singlelayer(timeFactor,  geom,  Fg);
        //growthfunction_beam_multilayer(timeFactor, geom,  Fg);
        //growthfunction_isotropic(1, ndim,  growthFactor,  Fg);
        FgInv = Fg.inverse();
        Jg = Fg.determinant();

        // elastic part of the deformation gradient
        Fe = F*FgInv;

        pbar = 0.0;

        val = pres - pbar;
        fact = val*val;

        elemError += ( fact * dvol );
    }//gp
  }

  return 0;
}




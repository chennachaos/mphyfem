
#include "BernsteinElem3DElecMechTet220.h"
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



BernsteinElem3DElecMechTet220::BernsteinElem3DElecMechTet220()
{
  ndim   = 3;
  degree = 2;
  npElem = 10;
  nlbfU  = 10;
  nlbfF  = 10;
  nlbfP  = 1;
  ndof   = 3;
  nsize  = npElem*ndof;

  Kup.resize(nsize, 1);
  Kpu.resize(1, nsize);
  Kpp.resize(1, 1);
  Rp.resize(1);

  pres = presPrev = 0.0;
  AlgoType = 2;
}


BernsteinElem3DElecMechTet220::~BernsteinElem3DElecMechTet220()
{
}


void BernsteinElem3DElecMechTet220::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double BernsteinElem3DElecMechTet220::computeVolume(bool configflag)
{
    double  Jac, param[3];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    elemVol=0.0;
    for(int gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        if(configflag)
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
        else
          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        elemVol += (gaussweights[gp]*Jac);
    }//gp

    return elemVol;
}



int BernsteinElem3DElecMechTet220::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != nsize) || (Mlocal.cols() != nsize) )
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = elmDat[5]*computeVolume(false)/10.0;

      for(int ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
    else
    {
      int  ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;
      double  fact, dvol, Jac, bb1, cc1, param[3];
      double  rho0  = elmDat[5] ;

      VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

      vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
      int nGPt=11;
      getGaussPointsTetra(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
  
      for(gp=0; gp<nGPt; gp++)
      {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];
  
          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
  
          dvol = gaussweights[gp]*(thick*Jac);
  
          for(ii=0;ii<nlbfU;ii++)
          {
            bb1 = (dvol*rho0)*N[ii];
  
            TI   = 3*ii;
            TIp1 = TI+1;
            TIp2 = TI+2;
  
            for(jj=0; jj<nlbfU; jj++)
            {
              TJ   = 3*jj;
              TJp1 = TJ+1;
              TJp2 = TJ+2;
  
              fact  = bb1*N[jj];
  
              Mlocal(TI,   TJ)    += fact ;
              Mlocal(TIp1, TJp1)  += fact ;
              Mlocal(TIp2, TJp2)  += fact ;
            }
          }
      } //gp
      //cout << " elemVol = " << elemVol << endl;
      //printMatrix(Klocal);  printf("\n\n\n");
    }

    return 0;
}



int BernsteinElem3DElecMechTet220::calcLoadVector(VectorXd& Flocal)
{
  return 0;
}



double  BernsteinElem3DElecMechTet220::calcCriticalTimeStep(bool flag)
{
    int  ii;

    double xNode[10], yNode[10], zNode[10];
    for(ii=0;ii<npElem;ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

      xNode[ii] = GeomData->NodePosNew[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosNew[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosNew[nodeNums[ii]][2];
    }

    double  xd, yd, zd, charlen=1.0e10;
    int  edge_nodes[6][2] = { {1,0}, {2,0}, {3,0}, {2,1}, {3,1}, {3,2}};

    for(ii=0; ii<6; ii++)
    {
      xd = xNode[edge_nodes[ii][0]] - xNode[edge_nodes[ii][1]];
      yd = yNode[edge_nodes[ii][0]] - yNode[edge_nodes[ii][1]]; 
      zd = zNode[edge_nodes[ii][0]] - zNode[edge_nodes[ii][1]];

      charlen = min(charlen, xd*xd+yd*yd+zd*zd);
    }
    charlen = 0.5*sqrt(charlen);

    double K    = matDat[0] ;
    double mu   = matDat[1] ;
    double rho  = elmDat[5] ;

    double  wave_speed = sqrt((K+4.0*mu/3.0)/rho);

    double  dtCric = charlen/wave_speed;

    return  dtCric;
}


int BernsteinElem3DElecMechTet220::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != 3) || (Mlocal.cols() != 3) )
      Mlocal.resize(3, 3);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = computeVolume(true)/3.0;

      for(int ii=0; ii<3; ii++)
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

    double xNode[3], yNode[3], xx, yy;
    for(ii=0;ii<3;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    int nGPt=3;
    getGaussPointsTetra(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    elemVolOrig=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, 1, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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

        for(ii=0;ii<3;ii++)
        {
          bb1 = dvol0*N[ii];

          for(jj=0; jj<3; jj++)
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


int  BernsteinElem3DElecMechTet220::calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), dNu_dz(nlbfU);
    VectorXd  Nf(nlbfU), dNf_dx(nlbfU), dNf_dy(nlbfU), dNf_dz(nlbfU);

    double  detF, trF, fact, fact1, dvol, dvol0, Jac, volstrain, phi, pbar, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, dUdJ, d2UdJ2;
    double  param[3], bforce[3], force[3], tarr[6],  tarr2[6];
    double  veloCur[3], acceCur[3], sig[3], Np[4];

    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), Bmat(3,9), Amat(3,3),  Gc(3,9);
    Amat.setZero();
    Bmat.setZero();
    Cmat.setZero();
    Gc.setZero();
    F.setZero();
    Fn.setZero();
    VectorXd  stre(9),  elecDisp(3),  elecField(3);
    stre.setZero();

    double BULK = matDat[0];
    double eps  = 0.0; if(!SolnData->TRULY_INCOMPRESSIBLE) eps=1.0/BULK;
    double mu   = matDat[1];
    double LAMBDA = BULK-r2d3*mu;
    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    bforce[2]   = elmDat[8]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

    double xNode[10], yNode[10], zNode[10], xx, yy,  zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

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

    elemVolOrig=0.0;
    elemVolCur=0.0;

    Kup.setZero();
    Kpu.setZero();
    Rp.setZero();

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), &dNu_dz(0), Jac);

        Nf     = Nu;
        dNf_dx = dNu_dx;
        dNf_dy = dNu_dy;
        dNf_dz = dNu_dz;

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F(0,0) = computeValueCur(0, dNu_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dNu_dy);                  // eps_xy
        F(0,2) = computeValueCur(0, dNu_dz);                  // eps_xz
        F(1,0) = computeValueCur(1, dNu_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dNu_dy) + 1.0;            // eps_yy
        F(1,2) = computeValueCur(1, dNu_dz);                  // eps_yz
        F(2,0) = computeValueCur(2, dNu_dx);                  // eps_zx
        F(2,1) = computeValueCur(2, dNu_dy);                  // eps_zy
        F(2,2) = computeValueCur(2, dNu_dz) + 1.0;            // eps_zz

        volstrain = F(0,0)+F(1,1)+F(2,2) - 3.0;

        detF = F.determinant();

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        acceCur[0] = computeValueDotDotCur(0, Nu);
        acceCur[1] = computeValueDotDotCur(1, Nu);
        acceCur[2] = computeValueDotDotCur(2, Nu);

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), &dNu_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;

          dNf_dx = dNu_dx;
          dNf_dy = dNu_dy;
          dNf_dz = dNu_dz;
        }
        elemVolCur  += dvol;


        // evaluate phi at the quadrature points
        phi = 0.0;
        elecField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          phi += (Nf[ii]*SolnData->var3[nodeNums[ii]]);

          elecField[0] -= (dNf_dx[ii]*SolnData->var3[nodeNums[ii]]);
          elecField[1] -= (dNf_dy[ii]*SolnData->var3[nodeNums[ii]]);
          elecField[2] -= (dNf_dz[ii]*SolnData->var3[nodeNums[ii]]);
        }

        // evaluate pressure at the quadrature point
        // pressure is element-wise constant for this element
        // and it is computed as a solution variable.
        // Nothing to do here

        MatlData->computeStressAndTangent(true, sss, Fn, F, elecField, pres, stre, elecDisp, Cmat, Bmat, Amat, ivar, gp, dt);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        // Calculate Stiffness and Residual
        //==============================================

        xx = yy = zz = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += Nu[ii]*xNode[ii];
          yy += Nu[ii]*yNode[ii];
          zz += Nu[ii]*zNode[ii];
        }

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];
        force[2] = rho0*bforce[2];

        //force[0] = analy.computeForce(0, xx, yy, 0.0) ;
        //force[1] = analy.computeForce(1, xx, yy, 0.0) ;

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];
        force[2] -= rho0*acceCur[2];

        d2UdJ2 = 1.0;

        for(ii=0; ii<nlbfU; ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb3 = dvol*dNu_dz[ii];
          bb4 = dvol*Nu[ii];
          bb5 = dvol0*Nu[ii];

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          //  Kuu matrix
          for(jj=0; jj<9; jj++)
          {
            Gc(0,jj) = bb1*Cmat(0,jj) + bb2*Cmat(3,jj) + bb3*Cmat(6,jj);
            Gc(1,jj) = bb1*Cmat(1,jj) + bb2*Cmat(4,jj) + bb3*Cmat(7,jj);
            Gc(2,jj) = bb1*Cmat(2,jj) + bb2*Cmat(5,jj) + bb3*Cmat(8,jj);
          }

          sig[0] = bb1*stre[0] + bb2*stre[3] + bb3*stre[6];
          sig[1] = bb1*stre[1] + bb2*stre[4] + bb3*stre[7];
          sig[2] = bb1*stre[2] + bb2*stre[5] + bb3*stre[8];

          FlocalU(TI)   += (bb5*force[0] - sig[0]) ;
          FlocalU(TIp1) += (bb5*force[1] - sig[1]) ;
          FlocalU(TIp2) += (bb5*force[2] - sig[2]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dNu_dx[jj];
            cc2 = dNu_dy[jj];
            cc3 = dNu_dz[jj];
            cc4 = Nu[jj];

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            // acceleration term
            acceFact2 = acceFact1*cc4*rho0;

            fact  = bb5*acceFact2;

            // material Stiffness contribution
            fact  += af*(sig[0]*cc1+sig[1]*cc2+sig[2]*cc3)*FiniteFact;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;
            Kuu(TIp2, TJp2)  += fact ;

            Kuu(TI,   TJ)    +=  af*(Gc(0,0)*cc1 + Gc(0,3)*cc2 + Gc(0,6)*cc3) ;
            Kuu(TI,   TJp1)  +=  af*(Gc(0,1)*cc1 + Gc(0,4)*cc2 + Gc(0,7)*cc3) ;
            Kuu(TI,   TJp2)  +=  af*(Gc(0,2)*cc1 + Gc(0,5)*cc2 + Gc(0,8)*cc3) ;

            Kuu(TIp1, TJ)    +=  af*(Gc(1,0)*cc1 + Gc(1,3)*cc2 + Gc(1,6)*cc3) ;
            Kuu(TIp1, TJp1)  +=  af*(Gc(1,1)*cc1 + Gc(1,4)*cc2 + Gc(1,7)*cc3) ;
            Kuu(TIp1, TJp2)  +=  af*(Gc(1,2)*cc1 + Gc(1,5)*cc2 + Gc(1,8)*cc3) ;

            Kuu(TIp2, TJ)    +=  af*(Gc(2,0)*cc1 + Gc(2,3)*cc2 + Gc(2,6)*cc3) ;
            Kuu(TIp2, TJp1)  +=  af*(Gc(2,1)*cc1 + Gc(2,4)*cc2 + Gc(2,7)*cc3) ;
            Kuu(TIp2, TJp2)  +=  af*(Gc(2,2)*cc1 + Gc(2,5)*cc2 + Gc(2,8)*cc3) ;
          }
          //printMatrix(Kuu);

          //  Kuf & Kfu matrices
          Gc(0,0) = bb1*Bmat(0,0) + bb2*Bmat(0,3) + bb3*Bmat(0,6);
          Gc(0,1) = bb1*Bmat(1,0) + bb2*Bmat(1,3) + bb3*Bmat(1,6);
          Gc(0,2) = bb1*Bmat(2,0) + bb2*Bmat(2,3) + bb3*Bmat(2,6);

          Gc(1,0) = bb1*Bmat(0,1) + bb2*Bmat(0,4) + bb3*Bmat(0,7);
          Gc(1,1) = bb1*Bmat(1,1) + bb2*Bmat(1,4) + bb3*Bmat(1,7);
          Gc(1,2) = bb1*Bmat(2,1) + bb2*Bmat(2,4) + bb3*Bmat(2,7);

          Gc(2,0) = bb1*Bmat(0,2) + bb2*Bmat(0,5) + bb3*Bmat(0,8);
          Gc(2,1) = bb1*Bmat(1,2) + bb2*Bmat(1,5) + bb3*Bmat(1,8);
          Gc(2,2) = bb1*Bmat(2,2) + bb2*Bmat(2,5) + bb3*Bmat(2,8);

          for(jj=0; jj<nlbfF; jj++)
          {
            tarr2[0] = Gc(0,0)*dNf_dx[jj] + Gc(0,1)*dNf_dy[jj] + Gc(0,2)*dNf_dz[jj];
            tarr2[1] = Gc(1,0)*dNf_dx[jj] + Gc(1,1)*dNf_dy[jj] + Gc(1,2)*dNf_dz[jj];
            tarr2[2] = Gc(2,0)*dNf_dx[jj] + Gc(2,1)*dNf_dy[jj] + Gc(2,2)*dNf_dz[jj];

            Kuf(TI,   jj)  += tarr2[0];
            Kuf(TIp1, jj)  += tarr2[1];
            Kuf(TIp2, jj)  += tarr2[2];

            Kfu(jj, TI)    += tarr2[0];
            Kfu(jj, TIp1)  += tarr2[1];
            Kfu(jj, TIp2)  += tarr2[2];
          }

          // Kup & Kpu matrices
          Kup(TI,   0)  += (af*bb1);
          Kup(TIp1, 0)  += (af*bb2);
          Kup(TIp2, 0)  += (af*bb3);

          Kpu(0, TI  )  += (af*bb1);
          Kpu(0, TIp1)  += (af*bb2);
          Kpu(0, TIp2)  += (af*bb3);
        }

        //  Kff matrix
        for(ii=0; ii<nlbfF; ii++)
        {
          bb1 = dvol*dNf_dx[ii];
          bb2 = dvol*dNf_dy[ii];
          bb3 = dvol*dNf_dz[ii];
          bb4 = dvol*Nf[ii];
          bb5 = dvol0*Nf[ii];

          FlocalF(ii) -= (bb1*elecDisp(0)+bb2*elecDisp(1)+bb3*elecDisp(2));

          tarr[0] = bb1*Amat(0,0) + bb2*Amat(1,0) + bb3*Amat(2,0);
          tarr[1] = bb1*Amat(0,1) + bb2*Amat(1,1) + bb3*Amat(2,1);
          tarr[2] = bb1*Amat(0,2) + bb2*Amat(1,2) + bb3*Amat(2,2);

          for(jj=0; jj<nlbfF; jj++)
          {
            Kff(ii, jj) -= (dNf_dx[jj]*tarr[0] + dNf_dy[jj]*tarr[1] + dNf_dz[jj]*tarr[2]);
          }
        }

        // entries for Kup matrix
        if(finite)
        {
          Rp(0)  += (detF - 1.0 - pres*eps)*dvol0;
        }
        else
        {
          Rp(0) += (volstrain - pres*eps)*dvol0;
        }

    }//gp

    // add contribution from the condensed matrix
    Kpp(0,0) = -af*elemVolOrig/BULK;

    Kuu -=  ( Kup*(Kpu/Kpp(0,0)) );
    FlocalU +=  ( Kup*(Rp(0)/Kpp(0,0)) );

    /*
    printMatrix(Kuu);

    printMatrix(Kuf);
    printMatrix(Kfu);
    printMatrix(Kff);

    printMatrix(Kup);
    printMatrix(Kpu);
    printMatrix(Kpp);

    printVector(FlocalU);
    printVector(FlocalF);
    printVector(FlocalP);
    */

    return 0;
}




int BernsteinElem3DElecMechTet220::solveForPressure()
{
    double  volstrain = 0.0, volstrainPrev=0.0;
    double  BULK = matDat[1];

   // explicit scheme
    if(SolnData->tis >= 400)
    {
      if(finite)
      {
        int gp, ii, ll=0, err, isw, count;
        VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

        double  F[9], detF, param[3], dvol0, dvol, Jac, pbar;
        double  cc[6][6], stre[6], dt, fact;

        vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
        getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

        Jbar = 0.0;
        pres = 0.0;
        elemVolOrig = 0.0;
        for(gp=0; gp<nGP; gp++)
        {
            param[0] = gausspoints1[gp];
            param[1] = gausspoints2[gp];
            param[2] = gausspoints3[gp];

            GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

            dvol0 = gaussweights[gp]*Jac;
            elemVolOrig  += dvol0;

            F[0] = computeValue(0, dN_dx) + 1.0;
            F[3] = computeValue(0, dN_dy);
            F[6] = computeValue(0, dN_dz);

            F[1] = computeValue(1, dN_dx);
            F[4] = computeValue(1, dN_dy) + 1.0;
            F[7] = computeValue(1, dN_dz);

            F[2] = computeValue(2, dN_dx);
            F[5] = computeValue(2, dN_dy);
            F[8] = computeValue(2, dN_dz) + 1.0;

            detF = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[3]*(F[1]*F[8]-F[2]*F[7]) + F[6]*(F[1]*F[5]-F[2]*F[4]);

            pres += (0.5*(detF-1.0/detF));
            //pres += (detF-1.0)*dvol0;
            //pres += 0.5*(detF-1.0/detF)*dvol0;
            //pres += (log(detF)/detF)*dvol0;
        }
        pres = BULK*pres/elemVolOrig;
      }
      else
      {
        volstrain  = computeValue(0, dNc_dx);
        volstrain += computeValue(1, dNc_dy);
        volstrain += computeValue(2, dNc_dz);

        pres = BULK*volstrain;
      }
    }
    else // implicit scheme
    {
      //volstrain  = computeValueIncr(0, dNc_dx);
      //volstrain += computeValueIncr(1, dNc_dy);
      //volstrain += computeValueIncr(2, dNc_dz);

      //pres += BULK*(volstrain+Rp)/elemVolOrig;

      VectorXd  dispIncr(nsize);
      int ii, kk, ind1, ind2;
      for(ii=0; ii<npElem; ii++)
      {
        ind1 = ii*3;
        ind2 = nodeNums[ii]*3;

        for(kk=0; kk<3; kk++)
        {
          dispIncr(ind1+kk)  =  SolnData->var1Incr[ind2+kk];
        }
      }

      presPrev = pres;

      VectorXd  vtmp = Kpu*dispIncr;

      pres += (-Rp(0) - vtmp(0))/Kpp(0,0);
    }

    //cout <<  "pres = " <<  pres <<  endl;

    return 0;
}



int BernsteinElem3DElecMechTet220::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, volstrain, pbar;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  stre[4], streDev[4], cc[4][4], param[2], bforce[2], force[2], Np[3];

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

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal1.rows() <= nsize)   Flocal1.resize(nsize);
    Flocal1.setZero();
    if(Flocal2.rows() <= 3)   Flocal2.resize(3);
    Flocal2.setZero();
    if(Kpu.rows() <= 3 || Kpu.cols() <= nsize)   Kpu.resize(3, nsize);
    Kpu.setZero();
    if(Kpp.rows() <= 3 || Kpp.cols() <= 3)   Kpp.resize(3, 3);
    Kpp.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)
          dvol *= F33;

        // basis functions for pressure dof
        Np[0] = 1.0-param[0]-param[1];    Np[1] = param[0];    Np[2] = param[1];

        // evaluate pressure at the quadrature points
        pres = 0.0;
        for(ii=0; ii<3; ii++)
          pres += (Np[ii]*SolnData->var3[nodeNums[ii]]);

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
        }

        if(finite)
        {
          fact = detF-1.0-pres*eps;
        }
        else
        {
          fact = volstrain-pres*eps;
        }

        for(ii=0; ii<3; ii++)
        {
          bb4 = dvol*Np[ii];
          bb5 = dvol0*Np[ii];

          Flocal2(ii) -= bb4*fact;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Kpu(ii, TJ)    += bb4*dN_dx[jj];
            Kpu(ii, TJp1)  += bb4*dN_dy[jj];
          }

          bb5 *= eps;
          for(jj=0; jj<3; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

    return 0;
}



int BernsteinElem3DElecMechTet220::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TIp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  F[9], detF, F33, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[6], cc[6][6], param[3], bforce[3], force[3];

    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    bforce[2]   = elmDat[8]*timeFunction[0].prop ;
    double dt   = myTime.dt;

    double xNode[10], yNode[10], zNode[10], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);


    elemVol=0.0;
    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0;gp<nGP;gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;
        elemVol += dvol;

        //for(ii=0;ii<nlbfU;ii++)
          //cout << dN_dx(ii) << '\t' << dN_dy(ii) << '\t' << dN_dz(ii) << endl;

        xx = yy = zz = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += xNode[ii]*N[ii];
          yy += yNode[ii]*N[ii];
          zz += zNode[ii]*N[ii];
        }

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[3] = computeValue(0, dN_dy);
        F[6] = computeValue(0, dN_dz);

        F[1] = computeValue(1, dN_dx);
        F[4] = computeValue(1, dN_dy) + 1.0;
        F[7] = computeValue(1, dN_dz);

        F[2] = computeValue(2, dN_dx);
        F[5] = computeValue(2, dN_dy);
        F[8] = computeValue(2, dN_dz) + 1.0;

        //F[0] = computeValuePrev(0, dN_dx) + 1.0;
        //F[3] = computeValuePrev(0, dN_dy);
        //F[6] = computeValuePrev(0, dN_dz);

        //F[1] = computeValuePrev(1, dN_dx);
        //F[4] = computeValuePrev(1, dN_dy) + 1.0;
        //F[7] = computeValuePrev(1, dN_dz);

        //F[2] = computeValuePrev(2, dN_dx);
        //F[5] = computeValuePrev(2, dN_dy);
        //F[8] = computeValuePrev(2, dN_dz) + 1.0;

        //cout << endl;
        //cout << (F[0]-1.0) << '\t' << (F[4]-1.0) << '\t' << (F[8]-1.0) << endl;
        //cout << 0.5*(F[1]+F[3]) << '\t' << 0.5*(F[7]+F[5]) << '\t' << 0.5*(F[6]+F[2]) << endl;
        //cout << endl;

        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

        if(detF < 0.0)   return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol = gaussweights[gp]*Jac;
        }

        //matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;


        // Calculate Residual
        //==============================================

        /*
        fact = -2.0e6;
        double r1 = sqrt(xx*xx+yy*yy);
        double r2 = sqrt(xx*xx+yy*yy+zz*zz);
        double theta = atan(yy/xx);
        if( xx < 0.0 )
          theta += PI;
        if( CompareDoubles(xx,0.0) && CompareDoubles(yy,0.0) )
          theta = 0.0;
        double phi   = atan(r1/zz);
        if( zz < 0.0)
          phi += PI;

        bforce[0]  =  fact*cos(theta)*sin(phi);
        bforce[1]  =  fact*sin(theta)*sin(phi);
        bforce[2]  =  fact*cos(phi);
        */

        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];
        force[2] = rho0*bforce[2];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*dN_dz[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          Flocal(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3] - bb3*stre[5]) ;
          Flocal(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1] - bb3*stre[4]) ;
          Flocal(TIp2) += (bb5*force[2] - bb1*stre[5] - bb2*stre[4] - bb3*stre[2]) ;
        }
    }//gp
    //cout << " elemVol = " << elemVol << endl;
    //printMatrix(Klocal);  printf("\n\n\n");
    //printVector(Flocal);

    return 0;
}



int BernsteinElem3DElecMechTet220::calcInternalForces()
{
  return 0;
}



void BernsteinElem3DElecMechTet220::elementContourplot(int vartype, int varindex, int index)
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


void BernsteinElem3DElecMechTet220::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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
    stressrecovery_extrapolate_Tetrahedron(degree, outval, &vals2project[0]);

    /*
    //else
    //{
       vector<double>  vecTemp;
       vecTemp = vals2project;

       vals2project[0] = vecTemp[0];
       vals2project[1] = vecTemp[1];
       vals2project[2] = vecTemp[2];
       vals2project[3] = vecTemp[3];
       vals2project[4] = 0.5*vecTemp[4] + 0.25*vecTemp[0] + 0.25*vecTemp[1];
       vals2project[5] = 0.5*vecTemp[5] + 0.25*vecTemp[1] + 0.25*vecTemp[2];
       vals2project[6] = 0.5*vecTemp[6] + 0.25*vecTemp[0] + 0.25*vecTemp[2];
       vals2project[7] = 0.5*vecTemp[7] + 0.25*vecTemp[0] + 0.25*vecTemp[3];
       vals2project[8] = 0.5*vecTemp[8] + 0.25*vecTemp[1] + 0.25*vecTemp[3];
       vals2project[9] = 0.5*vecTemp[9] + 0.25*vecTemp[2] + 0.25*vecTemp[3];
    //}
    */

    //printVector(vals2project);

    return;
}


void BernsteinElem3DElecMechTet220::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    double  detF, dvol, Jac, fact, dvol0, pbar, param[3], Np[4];
    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), Bmat(3,9), Amat(3,3);
    VectorXd  stre(9),  elecDisp(3),  elecField(3);
    F.setZero();
    Fn.setZero();
    stre.setZero();

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    int   err,  isw,  count,  count1, ll, ii, jj, gp;

    double dt = myTime.dt;

    double xNode[10], yNode[10], zNode[10], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_TETRA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        F(0,0) = computeValue(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValue(0, dN_dy);                  // eps_xy
        F(0,2) = computeValue(0, dN_dz);                  // eps_xz
        F(1,0) = computeValue(1, dN_dx);                  // eps_yx
        F(1,1) = computeValue(1, dN_dy) + 1.0;            // eps_yy
        F(1,2) = computeValue(1, dN_dz);                  // eps_yz
        F(2,0) = computeValue(2, dN_dx);                  // eps_zx
        F(2,1) = computeValue(2, dN_dy);                  // eps_zy
        F(2,2) = computeValue(2, dN_dz) + 1.0;            // eps_zz

        detF = F.determinant();

        // evaluate pressure at the quadrature point

        MatlData->computeStressAndTangent(false, sss, Fn, F, elecField, pres, stre, elecDisp, Cmat, Bmat, Amat, ivar, gp, dt);
        count++;
        ll += nivGP;

        if(varindex < 9)
           outval[gp] = stre[varindex];
        else if(varindex == 9)
           outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*(stre[1]*stre[1]+stre[2]*stre[2]+stre[5]*stre[5]))/2.0);
        else if(varindex == 10)
           //outval[gp] = (stre[0]+stre[1]+stre[2])/3.0;
           outval[gp] = pres;
        else
        {
           VectorXcd eivals = F.eigenvalues();

           outval[gp] = eivals[0].real();
           //outval[gp] = 0.0;
        }
    } //gp

    return;
}



void BernsteinElem3DElecMechTet220::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void BernsteinElem3DElecMechTet220::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp
    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}



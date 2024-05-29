
#include "MagnetoMech3DHex221.h"
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



MagnetoMech3DHex221::MagnetoMech3DHex221()
{
  ELEM_SHAPE = ELEM_SHAPE_HEXA_BERNSTEIN;

  ndim   = 3;
  degree = 2;
  npElem = 27;
  nlbfU  = 27;
  nlbfF  = 27;
  nlbfP  = 8;
  ndof   = 3;
  nsize  = npElem*ndof;
}


MagnetoMech3DHex221::~MagnetoMech3DHex221()
{
}


void MagnetoMech3DHex221::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double MagnetoMech3DHex221::computeVolume(bool configflag)
{
    double  Jac, param[3];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    elemVol=0.0;
    for(int gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        if(configflag)
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
        else
          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        elemVol += (gaussweights[gp]*Jac);
    }//gp

    return elemVol;
}



int MagnetoMech3DHex221::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
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
      int nGPt=27;
      getGaussPointsHexa(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);

      for(gp=0; gp<nGPt; gp++)
      {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];

          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

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



int MagnetoMech3DHex221::calcLoadVector(VectorXd& Flocal)
{
  return 0;
}



double  MagnetoMech3DHex221::calcCriticalTimeStep(bool flag)
{
    int  ii;

    double xNode[27], yNode[27], zNode[27];
    for(ii=0;ii<npElem;ii++)
    {
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

    //double  wave_speed = sqrt((K+4.0*mu/3.0)/rho);
    double  wave_speed = sqrt(mu/rho);

    double  dtCric = charlen/wave_speed;

    return  dtCric;
}


int MagnetoMech3DHex221::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != 3) || (Mlocal.cols() != 3) )
      Mlocal.resize(3, 3);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = computeVolume(true)/8.0;

      for(int ii=0; ii<8; ii++)
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
    getGaussPointsHexa(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    elemVolOrig=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_HEXA_BERNSTEIN, 1, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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






int MagnetoMech3DHex221::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU,  VectorXd& FlocalP, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), dNu_dz(nlbfU);
    VectorXd  Nf(nlbfU), dNf_dx(nlbfU), dNf_dy(nlbfU), dNf_dz(nlbfU);

    double  detF, trF, fact, fact1, dvol, dvol0, Jac, volstrain, phi, pbar, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, dUdJ, d2UdJ2;
    double  param[3], bforce[3], force[3], tarr[6],  tarr2[6];
    double  veloCur[3], acceCur[3], sig[3], Np[8];

    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), Bmat(3,9), Amat(3,3),  Gc(3,9);
    Amat.setZero();
    Bmat.setZero();
    Cmat.setZero();
    Gc.setZero();
    F.setZero();
    Fn.setZero();

    VectorXd  stre(9),  magnDisp(3),  magnField(3);
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
    double dt   = myTime.dt;

    double xNode[27], yNode[27], zNode[27], xx, yy,  zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(Kup.rows() < nsize || Kup.cols() < nlbfP)     Kup.resize(nsize, nlbfP);
    Kup.setZero();
    if(Kpu.rows() < nlbfP || Kpu.cols() < nsize)     Kpu.resize(nlbfP, nsize);
    Kpu.setZero();
    if(Kpp.rows() < nlbfP || Kpp.cols() < nlbfP)     Kpp.resize(nlbfP, nlbfP);
    Kpp.setZero();

    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalP.rows() < nlbfP)   FlocalP.resize(nlbfP);
    FlocalP.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), &dNu_dz(0), Jac);

        Nf     = Nu;
        dNf_dx = dNu_dx;
        dNf_dy = dNu_dy;
        dNf_dz = dNu_dz;

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F(0,0) = computeValuePrev(0, dNu_dx) + 1.0;            // eps_xx
        F(0,1) = computeValuePrev(0, dNu_dy);                  // eps_xy
        F(0,2) = computeValuePrev(0, dNu_dz);                  // eps_xz
        F(1,0) = computeValuePrev(1, dNu_dx);                  // eps_yx
        F(1,1) = computeValuePrev(1, dNu_dy) + 1.0;            // eps_yy
        F(1,2) = computeValuePrev(1, dNu_dz);                  // eps_yz
        F(2,0) = computeValuePrev(2, dNu_dx);                  // eps_zx
        F(2,1) = computeValuePrev(2, dNu_dy);                  // eps_zy
        F(2,2) = computeValuePrev(2, dNu_dz) + 1.0;            // eps_zz

        volstrain = F(0) + F(4) + F(8) - 3.0;

        detF = F.determinant();

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), &dNu_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;

          dNf_dx = dNu_dx;
          dNf_dy = dNu_dy;
          dNf_dz = dNu_dz;
        }
        elemVolCur  += dvol;


        // evaluate phi at the quadrature points
        magnField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          magnField[0] -= (dNf_dx[ii]*SolnData->var3[nodeNums[ii]]);
          magnField[1] -= (dNf_dy[ii]*SolnData->var3[nodeNums[ii]]);
          magnField[2] -= (dNf_dz[ii]*SolnData->var3[nodeNums[ii]]);
        }

        // evaluate pressure at the quadrature point
        Np[0] = 1.0-param[0]-param[1]-param[2];    Np[1]=param[0];    Np[2]=param[1];    Np[3]=param[2];
        pres = 0.0;
        for(ii=0; ii<nlbfP; ii++)
        {
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);
        }


        MatlData->computeStressAndTangent(true, sss, Fn, F, magnField, pres, stre, magnDisp, Cmat, Bmat, Amat, ivar, gp, dt);
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

          sig[0] = bb1*stre[0] + bb2*stre[3] + bb3*stre[6];
          sig[1] = bb1*stre[1] + bb2*stre[4] + bb3*stre[7];
          sig[2] = bb1*stre[2] + bb2*stre[5] + bb3*stre[8];

          FlocalU(TI)   += (bb5*force[0] - sig[0]) ;
          FlocalU(TIp1) += (bb5*force[1] - sig[1]) ;
          FlocalU(TIp2) += (bb5*force[2] - sig[2]) ;

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
            Kup(TIp2, jj)  += (bb3*Np[jj]);
          }
        }

        if(finite)
        {
          //fact = (dUdJ-pres)*eps;
          fact = detF - 1.0 - pres*eps;
        }
        else
        {
          fact = volstrain-pres*eps;
        }

        for(ii=0; ii<nlbfP; ii++)
        {
          bb4 = dvol*Np[ii]*d2UdJ2;
          bb5 = dvol0*Np[ii];

          FlocalP(ii) -= bb5*fact;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            Kpu(ii, TJ)    += (bb4*dNu_dx[jj]);
            Kpu(ii, TJp1)  += (bb4*dNu_dy[jj]);
            Kpu(ii, TJp2)  += (bb4*dNu_dz[jj]);
          }

          bb5 *= eps;
          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

    return 0;
}




int  MagnetoMech3DHex221::calcStiffnessAndResidualElecField(MatrixXd& Kff, VectorXd& FlocalF, bool firstIter)
{
    int   count,  count1, index, ii, jj, kk, gp;

    VectorXd  Nf(nlbfU), dNf_dx(nlbfU), dNf_dy(nlbfU), dNf_dz(nlbfU);

    double  detF, trF, fact, fact1, dvol, dvol0, Jac, volstrain, phi, pbar, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, dUdJ, d2UdJ2;
    double  param[3], bforce[3], force[3], tarr[6];

    MatrixXd  Amat(3,3);
    Amat.setZero();
    VectorXd  magnDisp(3),  magnField(3);

    double xNode[27], yNode[27], zNode[27], xx, yy,  zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(Kff.rows() < nlbfF || Kff.cols() < nlbfF)     Kff.resize(nlbfF, nlbfF);
    Kff.setZero();

    if(FlocalF.rows() < nlbfF)   FlocalF.resize(nlbfF);
    FlocalF.setZero();

    elemVolCur=0.0;

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];


        GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &Nf(0), &dNf_dx(0), &dNf_dy(0), &dNf_dz(0), Jac);

        dvol = gaussweights[gp]*Jac;

        elemVolCur  += dvol;

        // evaluate phi at the quadrature points
        magnField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          magnField[0] -= (dNf_dx[ii]*SolnData->var3[nodeNums[ii]]);
          magnField[1] -= (dNf_dy[ii]*SolnData->var3[nodeNums[ii]]);
          magnField[2] -= (dNf_dz[ii]*SolnData->var3[nodeNums[ii]]);
        }

        MatlData->computeElectricComponents(magnField, magnDisp, Amat);

        //cout <<  "Amat ......" <<  endl;
        //printMatrix(Amat);
        //printVector(magnField);

        // Calculate Stiffness and Residual
        //==============================================

        //  Kff matrix
        for(ii=0; ii<nlbfF; ii++)
        {
          bb1 = dvol*dNf_dx[ii];
          bb2 = dvol*dNf_dy[ii];
          bb3 = dvol*dNf_dz[ii];
          bb4 = dvol*Nf[ii];
          bb5 = dvol0*Nf[ii];

          FlocalF(ii) -= (bb1*magnDisp(0)+bb2*magnDisp(1)+bb3*magnDisp(2));

          tarr[0] = bb1*Amat(0,0) + bb2*Amat(1,0) + bb3*Amat(2,0);
          tarr[1] = bb1*Amat(0,1) + bb2*Amat(1,1) + bb3*Amat(2,1);
          tarr[2] = bb1*Amat(0,2) + bb2*Amat(1,2) + bb3*Amat(2,2);

          for(jj=0; jj<nlbfF; jj++)
          {
            Kff(ii, jj) -= (dNf_dx[jj]*tarr[0] + dNf_dy[jj]*tarr[1] + dNf_dz[jj]*tarr[2]);
          }
        }
    }//gp

    //printMatrix(Kff);
    //printVector(FlocalF);

    return 0;
}




int MagnetoMech3DHex221::calcResidual(VectorXd& Flocal)
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

    double xNode[27], yNode[27], zNode[27], xx, yy, zz;
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
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);


    elemVol=0.0;
    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0;gp<nGP;gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

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
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

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



int MagnetoMech3DHex221::calcInternalForces()
{
  return 0;
}


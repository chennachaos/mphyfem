
#include "MagnetoMech2DTri221.h"
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



MagnetoMech2DTri221::MagnetoMech2DTri221()
{
  ndim   = 2;
  degree = 2;
  npElem = 6;
  nlbfU  = 6;
  nlbfF  = 6;
  nlbfP  = 3;
  ndof   = 2;
  nsize  = npElem*ndof;

  dNc_dx.resize(nlbfU);
  dNc_dy.resize(nlbfU);

  pres = presPrev = 0.0;
}


MagnetoMech2DTri221::~MagnetoMech2DTri221()
{
}


void MagnetoMech2DTri221::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}




double MagnetoMech2DTri221::computeVolume(bool init)
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

  return elemVol;
}



int MagnetoMech2DTri221::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
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
    //cout << " elemVolOrig = " << elemVolOrig << endl;
    //printMatrix(Klocal);  printf("\n\n\n");
  }

  return 0;
}



int MagnetoMech2DTri221::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
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

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int nGPt=3;
    getGaussPointsTriangle(nGPt, gausspoints1, gausspoints2, gaussweights);

    elemVolOrig=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, 1, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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


double  MagnetoMech2DTri221::calcCriticalTimeStep(bool flag)
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






int  MagnetoMech2DTri221::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU,  VectorXd& FlocalP, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU);
    VectorXd  Nf(nlbfU), dNf_dx(nlbfU), dNf_dy(nlbfU);

    double  detF, trF, fact, fact1, dvol, dvol0, Jac, volstrain, phi, pbar, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, dUdJ, d2UdJ2;
    double  param[3], bforce[3], force[3], tarr[6],  tarr2[6];
    double  veloCur[3], acceCur[3], sig[3], Np[3];
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

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F(0,0) = computeValuePrev(0, dNu_dx) + 1.0;            // eps_xx
        F(0,1) = computeValuePrev(1, dNu_dx);                  // eps_yx
        F(1,0) = computeValuePrev(0, dNu_dy);                  // eps_xy
        F(1,1) = computeValuePrev(1, dNu_dy) + 1.0;            // eps_yy

        volstrain = F(0,0) + F(1,1) - 2.0;

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
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
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);

          Nf     = Nu;
          dNf_dx = dNu_dx;
          dNf_dy = dNu_dy;
        }
        elemVolCur  += dvol;


        if(finite)
          dvol *= F(2,2);

        // evaluate phi at the quadrature point

        elecField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          elecField[0] -= (dNf_dx[ii]*SolnData->var3[nodeNums[ii]]);
          elecField[1] -= (dNf_dy[ii]*SolnData->var3[nodeNums[ii]]);
        }

        // evaluate pressure at the quadrature point
        Np[0] = 1.0-param[0]-param[1];    Np[1]=param[0];    Np[2]=param[1];
        pres = 0.0;
        for(ii=0; ii<nlbfP; ii++)
        {
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);
        }

        MatlData->computeStressAndTangent(true, sss, Fn, F, elecField, pres, stre, elecDisp, Cmat, Bmat, Amat, ivar, gp, dt);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;


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
        //force[0] = 0.0;
        //force[1] = 0.0;

        //force[0] = analy.computeForce(0, xx, yy, 0.0) ;
        //force[1] = analy.computeForce(1, xx, yy, 0.0) ;

        d2UdJ2 = 1.0;

        for(ii=0; ii<nlbfU; ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb4 = dvol*Nu[ii];
          bb5 = dvol0*Nu[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[1] + bb2*stre[4];

          FlocalU(TI)   += (bb5*force[0] - sig[0]) ;
          FlocalU(TIp1) += (bb5*force[1] - sig[1]) ;

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
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
          bb4 = dvol*Np[ii];
          bb5 = dvol0*Np[ii];

          FlocalP(ii) -= bb5*fact;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Kpu(ii, TJ)    += (bb4*dNu_dx[jj]);
            Kpu(ii, TJp1)  += (bb4*dNu_dy[jj]);
          }

          bb5 *= eps;
          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

    /*
    printMatrix(Kup);
    printMatrix(Kpu);
    printMatrix(Kpp);

    printVector(FlocalU);
    printVector(FlocalP);
    */

    return 0;
}




int  MagnetoMech2DTri221::calcStiffnessAndResidualElecField(MatrixXd& Kff, VectorXd& FlocalF, bool firstIter)
{
    int   count,  count1, index, ii, jj, kk, ll, gp;

    VectorXd  Nf(nlbfF), dNf_dx(nlbfF), dNf_dy(nlbfF);

    double  fact, fact1, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[3], tarr[6];
    MatrixXd  Amat(3,3);
    Amat.setZero();
    VectorXd  F(9), stre(9),  elecDisp(3),  elecField(3);

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
    if(Kff.rows() < nlbfF || Kff.cols() < nlbfF)     Kff.resize(nlbfF, nlbfF);
    Kff.setZero();

    if(FlocalF.rows() < nlbfF)   FlocalF.resize(nlbfF);
    FlocalF.setZero();

    elemVolCur=0.0;

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nf(0), &dNf_dx(0), &dNf_dy(0), Jac);
        dvol = gaussweights[gp]*(Jac*thick);

        elemVolCur  += dvol;

        // evaluate phi at the quadrature point

        elecField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          elecField[0] -= (dNf_dx[ii]*SolnData->var3[nodeNums[ii]]);
          elecField[1] -= (dNf_dy[ii]*SolnData->var3[nodeNums[ii]]);
        }

        MatlData->computeElectricComponents(elecField, elecDisp, Amat);

        //cout <<  "Amat ......" <<  endl;
        //printMatrix(Amat);
        //printVector(elecField);

        // Calculate Stiffness and Residual
        //==============================================

        for(ii=0; ii<nlbfF; ii++)
        {
          bb1 = dvol*dNf_dx[ii];
          bb2 = dvol*dNf_dy[ii];
          bb4 = dvol*Nf[ii];
          bb5 = dvol0*Nf[ii];

          FlocalF(ii) -= (bb1*elecDisp(0)+bb2*elecDisp(1));

          tarr[0] = bb1 * Amat(0,0) + bb2 * Amat(1,0);
          tarr[1] = bb1 * Amat(0,1) + bb2 * Amat(1,1);

          for(jj=0; jj<nlbfF; jj++)
          {
            Kff(ii, jj) -= (dNf_dx[jj]*tarr[0] + dNf_dy[jj]*tarr[1]);
          }
        }
    }//gp

    //printMatrix(Kff);
    //printVector(FlocalF);

    return 0;
}








int MagnetoMech2DTri221::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ii, jj, ll, gp, TI, TIp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, pbar, volstrain, r2d3 = 2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[4], cc[4][4], param[2], bforce[2], force[2], Np[3];

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



int MagnetoMech2DTri221::calcResidualPressure(VectorXd& Flocal)
{
    int gp, ii;
    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, volstrain, param[2], dvol0, dvol;
    double  fact, Jac, r2d3 = 2.0/3.0, dUdJ, Np[3];
    double  BULK = matDat[0];
    double mu   = matDat[1];
    double LAMBDA = BULK-r2d3*mu;

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() != 3)
      Flocal.resize(3);
    Flocal.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

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
        Np[0] = 1.0-param[0]-param[1];    Np[1] = param[0];    Np[2] = param[1];

        fact *= dvol0;
        for(ii=0; ii<3; ii++)
        {
          Flocal(ii) += (fact*Np[ii]);
        }
    }

    return 0;
}


int MagnetoMech2DTri221::solveForPressure()
{
    cout << "'MagnetoMech2DTri221::solveForPressure()' ... not implemented yet " << endl;

    return 0;
}



int  MagnetoMech2DTri221::solveForPressureTIC(VectorXd& matK1, double beta, double dt, VectorXd& Flocal)
{
    cout << "'MagnetoMech2DTri221::solveForPressureTIC()' ... not implemented yet " << endl;

    return 0;
}



int MagnetoMech2DTri221::calcLoadVector(VectorXd& Flocal)
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



int MagnetoMech2DTri221::calcInternalForces()
{
  return 0;
}




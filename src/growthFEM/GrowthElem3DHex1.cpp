
#include "GrowthElem3DHex1.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "growthfunctions.h"


using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



GrowthElem3DHex1::GrowthElem3DHex1()
{
  ELEM_SHAPE = ELEM_SHAPE_HEXA_BERNSTEIN;

  ndim   = 3;
  degree = 1;
  npElem = 8;
  nlbfU  = 8;
  nlbfP  = 0;
  ndof   = 3;
  nsize  = nlbfU*ndof;
}

GrowthElem3DHex1::~GrowthElem3DHex1()
{
}


void GrowthElem3DHex1::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



int GrowthElem3DHex1::calcLoadVector(VectorXd& Flocal)
{
  return 0;
}



double GrowthElem3DHex1::computeVolume(bool init)
{
    double  dvol, Jac, param[3];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    int   ii, gp;

    double xNode[8], yNode[8], zNode[8], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    elemVol=0.0;
    for(gp=0;gp<nGP;gp++)
    {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];

          if(init)
            GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
          else
            GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol = gaussweights[gp]*Jac;

      elemVol += dvol;
  }//gp

  return elemVol;
}



int GrowthElem3DHex1::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
  // mass lumping - row-wise sum
  if(MassLumping)
  {
      double rho0  = elmDat[5] ;
      double fact  = rho0*computeVolume(true)/8.0;

      if(Mlocal.rows() != nsize)
        Mlocal.resize(nsize, nsize);
      Mlocal.setZero();

      for(int ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
  }
  else
  {
    int  ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  fact, dvol, Jac, bb1, cc1, param[3];

    double rho0  = elmDat[5] ;
    double rho  = rho0 ;

    double xNode[10], yNode[10], zNode[10], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    if(Mlocal.rows() != nsize)
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    int nGPt=9;
    getGaussPointsHexa(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    elemVol=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol = gaussweights[gp]*Jac;
        elemVol += dvol;

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


double  GrowthElem3DHex1::calcCriticalTimeStep(bool flag)
{
  double  dtCric;

  return  dtCric;
}


/*
int GrowthElem3DHex1::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  detF, fact, fact1, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[3], bforce[3], force[3];
    double  veloCur[3], acceCur[3], sig[3];

    MatrixXd  Fg(3,3), FgInv(3,3), F(3,3), Fn(3,3), Fe(3,3);
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    Fg.setZero();
    FgInv.setZero();

    double  timeFact = elmDat[10]*timeFunction[0].prop, Jg ;

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
    double tCur = myTime.cur;

    double xNode[8], yNode[8], zNode[8], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() <= nsize)   Flocal.resize(nsize);
    Flocal.setZero();
    if(Klocal.rows() < nsize || Klocal.cols() < nsize)   Klocal.resize(nsize, nsize);
    Klocal.setZero();

    MatrixXd  Cmat(9,9), Gc(3,9);
    VectorXd  stre(9);
    Gc.setZero();
    stre.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0;gp<nGP;gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = Jac * gaussweights[gp];
        dvol = dvol0;
        elemVolOrig += dvol0;

        geom[0] = geom[1] = geom[2] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += N[ii]*xNode[ii];
          geom[1] += N[ii]*yNode[ii];
          geom[2] += N[ii]*zNode[ii];
        }

        F(0,0) = computeValueCur(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dN_dy);                  // eps_xy
        F(0,2) = computeValueCur(0, dN_dz);                  // eps_xz
        F(1,0) = computeValueCur(1, dN_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dN_dy) + 1.0;            // eps_yy
        F(1,2) = computeValueCur(1, dN_dz);                  // eps_yz
        F(2,0) = computeValueCur(2, dN_dx);                  // eps_zx
        F(2,1) = computeValueCur(2, dN_dy);                  // eps_zy
        F(2,2) = computeValueCur(2, dN_dz) + 1.0;            // eps_zz

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);
        acceCur[2] = computeValueDotDotCur(2, N);

        detF = F.determinant();

        if(detF < 0.0)
        {
          cerr << " Negative Jacobian in the element " << elenum << "'\t' at quadrature point " << gp << endl;
          return 1;
        }

        growthfunction_circularplate_singlelayer(timeFact, geom, Fg);
        FgInv = Fg.inverse();
        Jg = Fg.determinant();

        // elastic part of the deformation gradient
        Fe = F*FgInv;

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0),  &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }
        elemVolCur  += dvol;

        pres = 0.0;
        MatlData->computeStressAndTangent(true, sss, Fn, Fe, pres, stre, Cmat);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        Cmat /=  Jg;
        stre /=  Jg;

        //if(elenum == 0)  printMatrix(Cmat);

        // Calculate Stiffness and Residual
        //==============================================

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];
        force[2] = rho0*bforce[2];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];
        force[2] -= rho0*acceCur[2];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*dN_dz[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          for(kk=0; kk<9; kk++)
          {
            Gc(0,kk) = bb1 * Cmat(0,kk) + bb2 * Cmat(3,kk) + bb3 * Cmat(6,kk);
            Gc(1,kk) = bb1 * Cmat(1,kk) + bb2 * Cmat(4,kk) + bb3 * Cmat(7,kk);
            Gc(2,kk) = bb1 * Cmat(2,kk) + bb2 * Cmat(5,kk) + bb3 * Cmat(8,kk);
          }

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          sig[0] = bb1*stre[0] + bb2*stre[3] + bb3*stre[6] ;
          sig[1] = bb1*stre[1] + bb2*stre[4] + bb3*stre[7] ;
          sig[2] = bb1*stre[2] + bb2*stre[5] + bb3*stre[8] ;

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;
          Flocal(TIp2) += (bb5*force[2] - sig[2]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = dN_dz[jj];
            cc5 = N[jj];

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            // acceleration term
            acceFact2 = acceFact1*cc5*rho0;

            fact  = bb5*acceFact2;

            // consistent mass matrix
            Klocal(TI,   TJ)    +=  fact;
            Klocal(TIp1, TJp1)  +=  fact;
            Klocal(TIp2, TJp2)  +=  fact;

            fact  = bb5*cc5*rho0;

            //Mlocal(TI,   TJ)    += fact ;
            //Mlocal(TIp1, TJp1)  += fact ;
            //Mlocal(TIp2, TJp2)  += fact ;

            // lumped mass matrix
            //Mlocal(TI,   TI)    +=  fact;
            //Mlocal(TIp1, TIp1)  +=  fact;
            //Mlocal(TIp2, TIp2)  +=  fact;

            // material Stiffness contribution
            fact = af*(sig[0]*cc1 + sig[1]*cc2 + sig[2]*cc3)*FiniteFact;

            Klocal(TI,   TJ)    +=  fact;
            Klocal(TIp1, TJp1)  +=  fact;
            Klocal(TIp2, TJp2)  +=  fact;

            Klocal(TI,   TJ)    +=  af*(Gc(0,0) * cc1 + Gc(0,3) * cc2 + Gc(0,6) * cc3) ;
            Klocal(TI,   TJp1)  +=  af*(Gc(0,1) * cc1 + Gc(0,4) * cc2 + Gc(0,7) * cc3) ;
            Klocal(TI,   TJp2)  +=  af*(Gc(0,2) * cc1 + Gc(0,5) * cc2 + Gc(0,8) * cc3) ;

            Klocal(TIp1, TJ)    +=  af*(Gc(1,0) * cc1 + Gc(1,3) * cc2 + Gc(1,6) * cc3) ;
            Klocal(TIp1, TJp1)  +=  af*(Gc(1,1) * cc1 + Gc(1,4) * cc2 + Gc(1,7) * cc3) ;
            Klocal(TIp1, TJp2)  +=  af*(Gc(1,2) * cc1 + Gc(1,5) * cc2 + Gc(1,8) * cc3) ;

            Klocal(TIp2, TJ)    +=  af*(Gc(2,0) * cc1 + Gc(2,3) * cc2 + Gc(2,6) * cc3) ;
            Klocal(TIp2, TJp1)  +=  af*(Gc(2,1) * cc1 + Gc(2,4) * cc2 + Gc(2,7) * cc3) ;
            Klocal(TIp2, TJp2)  +=  af*(Gc(2,2) * cc1 + Gc(2,5) * cc2 + Gc(2,8) * cc3) ;
          }
        }
  }//gp
  //cout << " elemVol = " << elemVol << endl;
  //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);
  //printVector(Flocal);

    //
    //bool MassLumping=false;
    bool MassLumping=true;
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

    Flocal -= Mlocal*accC;
    Klocal += acceFact1*Mlocal;
    //

    return 0;
}
*/



int GrowthElem3DHex1::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  detF, fact, fact1, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[3], bforce[3], force[3];
    double  veloCur[3], acceCur[3], sig[3], growth[3], Psi;

    MatrixXd  Fg(3,3), FgInv(3,3), F(3,3), Fn(3,3), Fe(3,3);
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    Fg.setZero();
    FgInv.setZero();

    double  timeFact = elmDat[10]*timeFunction[0].prop, Jg ;

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
    double tCur = myTime.cur;

    double xNode[8], yNode[8], zNode[8], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() <= nsize)   Flocal.resize(nsize);
    Flocal.setZero();
    if(Klocal.rows() < nsize || Klocal.cols() < nsize)   Klocal.resize(nsize, nsize);
    Klocal.setZero();

    MatrixXd  Gc(3,3);
    Gc.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    //cout << " lllllllll " << endl;
    //printVector(SolnData->dispTarget);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0;gp<nGP;gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = Jac * gaussweights[gp];
        dvol = dvol0;
        elemVolOrig += dvol0;

        //cout << " eeeeeeeeeeeeee " << endl;
        //cout << " Jac = " << Jac << endl;

        geom[0] = geom[1] = geom[2] = 0.0;
        F.setZero();
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += N[ii]*xNode[ii];
          geom[1] += N[ii]*yNode[ii];
          geom[2] += N[ii]*zNode[ii];

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          F(0,0)  +=  dN_dx[ii]*SolnData->dispTarget[TI];
          F(0,1)  +=  dN_dy[ii]*SolnData->dispTarget[TI];
          F(0,2)  +=  dN_dz[ii]*SolnData->dispTarget[TI];

          F(1,0)  +=  dN_dx[ii]*SolnData->dispTarget[TIp1];
          F(1,1)  +=  dN_dy[ii]*SolnData->dispTarget[TIp1];
          F(1,2)  +=  dN_dz[ii]*SolnData->dispTarget[TIp1];

          F(2,0)  +=  dN_dx[ii]*SolnData->dispTarget[TIp2];
          F(2,1)  +=  dN_dy[ii]*SolnData->dispTarget[TIp2];
          F(2,2)  +=  dN_dz[ii]*SolnData->dispTarget[TIp2];
        }

        cout << " eeeeeeeeeeeeee " << endl;

        F(0,0) += 1.0;
        F(1,1) += 1.0;
        F(2,2) += 1.0;
        
        printMatrix(F);

        growth[0] = computeValueCur(0, N);
        growth[1] = computeValueCur(1, N);
        growth[2] = computeValueCur(2, N);

        detF = F.determinant();

        if(detF < 0.0)
        {
          cerr << " Negative Jacobian in the element " << elenum << "'\t' at quadrature point " << gp << endl;
          return 1;
        }

        Fg(0,0) = 1.0+growth[0];
        Fg(1,1) = 1.0+growth[1];
        Fg(2,2) = 1.0+growth[2];

        printMatrix(Fg);

        FgInv = Fg.inverse();
        Jg = Fg.determinant();

        // elastic part of the deformation gradient
        Fe = F*FgInv;

        //if(finite)
        //{
          //GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0),  &dN_dz(0), Jac);
          //dvol = gaussweights[gp]*Jac;
        //}
        //elemVolCur  += dvol;

        //MatlData->computeStressAndTangent(true, sss, Fn, Fe, pres, stre, Cmat);
        Psi = MatlData->computeValue(sss, Fe);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        cout << " Psi = " << Psi << endl;

        //if(elenum == 0)  printMatrix(Cmat);

        // Calculate Stiffness and Residual
        //==============================================

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];
        force[2] = rho0*bforce[2];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*dN_dz[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii]*Psi;

          Gc(0,0) = 0.0;
          Gc(0,1) = Fg(2,2);
          Gc(0,2) = Fg(1,1);

          Gc(1,0) = Fg(2,2);
          Gc(1,1) = 0.0;
          Gc(1,2) = Fg(0,0);

          Gc(2,0) = Fg(1,1);
          Gc(2,1) = Fg(0,0);
          Gc(2,2) = 0.0;
          
          Gc *= bb5;

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          sig[0] = bb5*Fg(1,1)*Fg(2,2) ;
          sig[1] = bb5*Fg(0,0)*Fg(2,2) ;
          sig[2] = bb5*Fg(0,0)*Fg(1,1) ;

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;
          Flocal(TIp2) += (bb5*force[2] - sig[2]) ;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = dN_dz[jj];
            cc5 = N[jj];

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;
            
            //fact = bb5*cc5;

            //Klocal(TI,   TJ)    +=  fact;
            //Klocal(TIp1, TJp1)  +=  fact;
            //Klocal(TIp2, TJp2)  +=  fact;

            Klocal(TI,   TJ)    +=  (Gc(0,0) * cc5) ;
            Klocal(TI,   TJp1)  +=  (Gc(0,1) * cc5) ;
            Klocal(TI,   TJp2)  +=  (Gc(0,2) * cc5) ;

            Klocal(TIp1, TJ)    +=  (Gc(1,0) * cc5) ;
            Klocal(TIp1, TJp1)  +=  (Gc(1,1) * cc5) ;
            Klocal(TIp1, TJp2)  +=  (Gc(1,2) * cc5) ;

            Klocal(TIp2, TJ)    +=  (Gc(2,0) * cc5) ;
            Klocal(TIp2, TJp1)  +=  (Gc(2,1) * cc5) ;
            Klocal(TIp2, TJp2)  +=  (Gc(2,2) * cc5) ;
          }
        }
  }//gp
  //cout << " elemVol = " << elemVol << endl;
  printMatrix(Klocal);  //printf("\n\n\n");  printVector(Flocal);
  //printVector(Flocal);

    return 0;
}



int GrowthElem3DHex1::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp;
    int   TI, TIp1, TIp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  F[9], detF, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[6], cc[6][6], param[3], bforce[3], force[3];

    double rho  = elmDat[5] ;
    double rho0 = rho ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    bforce[2]   = elmDat[8]*timeFunction[0].prop ;
    double dt   = myTime.dt;

    double xNode[8], yNode[8], zNode[8], xx, yy, zz;
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

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;
        elemVol += dvol;

        F[0] = computeValuePrev(0, dN_dx) + 1.0;
        F[3] = computeValuePrev(0, dN_dy);
        F[6] = computeValuePrev(0, dN_dz);

        F[1] = computeValuePrev(1, dN_dx);
        F[4] = computeValuePrev(1, dN_dy) + 1.0;
        F[7] = computeValuePrev(1, dN_dz);

        F[2] = computeValuePrev(2, dN_dx);
        F[5] = computeValuePrev(2, dN_dy);
        F[8] = computeValuePrev(2, dN_dz) + 1.0;

        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

        if(detF < 0.0)   return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }

        //matlib3d__(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;


        // Calculate Residual
        //==============================================

        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];
        force[2] = rho0*bforce[2];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*dN_dz[ii];
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
  //printVector(Flocal);

  return 0;
}





int GrowthElem3DHex1::calcInternalForces()
{
  return 0;
}



void GrowthElem3DHex1::elementContourplot(int vartype, int varindex, int index)
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


void GrowthElem3DHex1::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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
    stressrecovery_extrapolate_Hexahedron(degree, outval, &vals2project[0]);
    //printVector(vals2project);

    return;
}


void GrowthElem3DHex1::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    double  F[9], detF, F33, fact, dvol, dvol0;
    double  Jac, bb1, bb2, bb3, cc1, cc2, cc3;
    double  stre[6], cc[6][6], param[3];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    int   err,  isw,  count,  count1, ll, ii, jj, gp;

    double dt = myTime.dt;

    double xNode[8], yNode[8], zNode[8], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[3] = computeValueCur(0, dN_dy);
        F[6] = computeValueCur(0, dN_dz);

        F[1] = computeValueCur(1, dN_dx);
        F[4] = computeValueCur(1, dN_dy) + 1.0;
        F[7] = computeValueCur(1, dN_dz);

        F[2] = computeValueCur(2, dN_dx);
        F[5] = computeValueCur(2, dN_dy);
        F[8] = computeValueCur(2, dN_dz) + 1.0;

        detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

        //matlib3d__(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);

        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", stre[0], stre[1], stre[2], stre[3]);

        count++;
        ll += nivGP;

        if(varindex < 6)
           outval[gp] = stre[varindex];
        else if(varindex == 6)
           outval[gp] = sqrt((pow(stre[0]-stre[1],2.0) + pow(stre[1]-stre[2], 2.0) + pow(stre[2]-stre[0], 2.0) + 6* stre[3]*stre[3])/2);
        else if(varindex == 7)
           outval[gp] = (stre[0]+stre[1]+stre[2])/3.0;
    } //gp

    return;
}






void GrowthElem3DHex1::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void GrowthElem3DHex1::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}


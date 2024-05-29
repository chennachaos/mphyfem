
#include "GrowthElem3DMixedBase.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "growthfunctions.h"
#include "utilitiesmaterialHyperelastic.h"

using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



GrowthElem3DMixedBase::GrowthElem3DMixedBase()
{
}



GrowthElem3DMixedBase::~GrowthElem3DMixedBase()
{
}


int GrowthElem3DMixedBase::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU), Np(nlbfP);

    double  detF, detFn, Jne, fact, fact1, dvol, dvol0, Jac, volstrain, r1d3=1.0/3.0, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Rp, dUdJ, d2UdJ2, Jhat, thetahat;
    double  param[3], bforce[3], force[3];
    double  veloCur[3], acceCur[3], sig[3];

    MatrixXd  Fg(3,3), FgInv(3,3), F(3,3), Fn(3,3), Fe(3,3), Fen(3,3);
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    Fen.setZero();
    Fg.setZero();
    FgInv.setZero();

    //cout << " iiiiiiiiiiii " << endl;

    //double  timeFactor = elmDat[10]*timeFunction[0].prop, timeFactorPrev = timeFactor-myTime.dt, Jg ;
    double  timeFactor = timeFunction[0].prop, timeFactorPrev = timeFactor-myTime.dt, Jg ;

    int  Utype  = MatlData->getUtype();                     // volumetric energy function
    double BULK = 1.0/MatlData->getKinv();

    double rho0 = MatlData->getDensity();
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

    double xNode[nlbfU], yNode[nlbfU], zNode[nlbfU], geom[3];
    for(ii=0;ii<nlbfU;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;

    if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
      getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);
    else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
      getGaussPointsWedge(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);
    else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
      getGaussPointsHexa(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);
    else
    {
        throw runtime_error("Element type no defined...");
    }


    // resize local matrices and initialise them to zero
    if(Flocal1.rows() <= nsize)   Flocal1.resize(nsize);
    Flocal1.setZero();
    if(Flocal2.rows() < nlbfP)   Flocal2.resize(nlbfP);
    Flocal2.setZero();
    if(Kup.rows() < nsize || Kup.cols() < nlbfP)       Kup.resize(nsize, nlbfP);
    Kup.setZero();
    if(Kpu.rows() < nlbfP || Kpu.cols() < nsize)       Kpu.resize(nlbfP, nsize);
    Kpu.setZero();
    if(Kpp.rows() < nlbfP || Kpp.cols() < nlbfP)       Kpp.resize(nlbfP, nlbfP);
    Kpp.setZero();
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)       Kuu.resize(nsize, nsize);
    Kuu.setZero();

    MatrixXd  Cmat(9,9), Gc(3,9);
    VectorXd  stre(9);
    Gc.setZero();
    stre.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    //cout << " iiiiiiiiiiii " << endl;
    
    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
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

        Fn(0,0) = computeValuePrev(0, dN_dx) + 1.0;            // eps_xx
        Fn(0,1) = computeValuePrev(0, dN_dy);                  // eps_xy
        Fn(0,2) = computeValuePrev(0, dN_dz);                  // eps_xz
        Fn(1,0) = computeValuePrev(1, dN_dx);                  // eps_yx
        Fn(1,1) = computeValuePrev(1, dN_dy) + 1.0;            // eps_yy
        Fn(1,2) = computeValuePrev(1, dN_dz);                  // eps_yz
        Fn(2,0) = computeValuePrev(2, dN_dx);                  // eps_zx
        Fn(2,1) = computeValuePrev(2, dN_dy);                  // eps_zy
        Fn(2,2) = computeValuePrev(2, dN_dz) + 1.0;            // eps_zz

        detFn = Fn.determinant();

        //growthfunction_beam_singlelayer(timeFactorPrev,  geom,  Fg);
        growthfunction_isotropic(3, ndim, timeFactorPrev,  Fg);

        FgInv = Fg.inverse();
        Fen   = Fn*FgInv;

        Jg = Fg.determinant();
        Jne = detFn / Jg;

        compute_constants_volumetric_functions(finite, Utype, Jne, BULK, Jhat, thetahat);


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

        volstrain = F(0,0)+F(1,1)+F(2,2) - 3.0;

        detF = F.determinant();

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        //cout << " jjjjjjjjjjjj " << endl;

        //growthfunction_beam_singlelayer(timeFactor,  geom,  Fg);
        //growthfunction_beam_singlelayer3D(timeFactor,  geom,  Fg);
        //growthfunction_circularplate_singlelayer(timeFactor, geom, Fg);
        //if(geom[2] > 0.2) fact = timeFactor;
        //else fact = 0.0;
        //fact = (geom[0]*sin(2.0*PI*geom[2])+geom[1]*cos(2.0*PI*geom[2]))*timeFactor;
        //fact = timeFactor*geom[2];
        //fact = timeFactor;
        growthfunction_isotropic(3, ndim, timeFactor,  Fg);
        //growthfunction_rod(1, ndim, timeFact, geom,  Fg);

        FgInv = Fg.inverse();
        Jg = Fg.determinant();

        // elastic part of the deformation gradient
        Fe = F*FgInv;

        //printMatrix(F);
        //printMatrix(Fg);
        //printMatrix(Fe);
        //cout << "Jg = " << Jg << endl;

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0),  &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }
        elemVolCur  += dvol;

        // basis functions for pressure dof

        if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
          BernsteinBasisFunsTetra(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
          BernsteinBasisFunsWedge(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
          BernsteinBasisFunsHexa(1, param[0], param[1], param[2], &Np(0));


        // evaluate pressure at the quadrature points
        pres = computeValue2Cur(0, Np);

        //cout << " jjjjjjjjjjjj " << endl;

        MatlData->computeStressAndTangent(true, sss, Fen, Fe, pres, stre, Cmat, ivar, gp, dt);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;


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

          Flocal1(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal1(TIp1) += (bb5*force[1] - sig[1]) ;
          Flocal1(TIp2) += (bb5*force[2] - sig[2]) ;

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
            Kuu(TI,   TJ)    +=  fact;
            Kuu(TIp1, TJp1)  +=  fact;
            Kuu(TIp2, TJp2)  +=  fact;

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

            Kuu(TI,   TJ)    +=  fact;
            Kuu(TIp1, TJp1)  +=  fact;
            Kuu(TIp2, TJp2)  +=  fact;

            Kuu(TI,   TJ)    +=  af*(Gc(0,0) * cc1 + Gc(0,3) * cc2 + Gc(0,6) * cc3) ;
            Kuu(TI,   TJp1)  +=  af*(Gc(0,1) * cc1 + Gc(0,4) * cc2 + Gc(0,7) * cc3) ;
            Kuu(TI,   TJp2)  +=  af*(Gc(0,2) * cc1 + Gc(0,5) * cc2 + Gc(0,8) * cc3) ;

            Kuu(TIp1, TJ)    +=  af*(Gc(1,0) * cc1 + Gc(1,3) * cc2 + Gc(1,6) * cc3) ;
            Kuu(TIp1, TJp1)  +=  af*(Gc(1,1) * cc1 + Gc(1,4) * cc2 + Gc(1,7) * cc3) ;
            Kuu(TIp1, TJp2)  +=  af*(Gc(1,2) * cc1 + Gc(1,5) * cc2 + Gc(1,8) * cc3) ;

            Kuu(TIp2, TJ)    +=  af*(Gc(2,0) * cc1 + Gc(2,3) * cc2 + Gc(2,6) * cc3) ;
            Kuu(TIp2, TJp1)  +=  af*(Gc(2,1) * cc1 + Gc(2,4) * cc2 + Gc(2,7) * cc3) ;
            Kuu(TIp2, TJp2)  +=  af*(Gc(2,2) * cc1 + Gc(2,5) * cc2 + Gc(2,8) * cc3) ;
          }

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*af*Np[jj]);
            Kup(TIp1, jj)  += (bb2*af*Np[jj]);
            Kup(TIp2, jj)  += (bb3*af*Np[jj]);
          }
        }

        if(finite)
        {
          Rp = detF - Jg*(Jhat+thetahat*pres);
        }
        else
        {
          Rp = volstrain - thetahat*pres;
        }


        for(ii=0; ii<nlbfP; ii++)
        {
          bb4 = dvol*Np[ii]*af;
          bb5 = dvol0*Np[ii];

          Flocal2(ii) -= bb5*Rp;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            Kpu(ii, TJ)    += bb4*dN_dx[jj];
            Kpu(ii, TJp1)  += bb4*dN_dy[jj];
            Kpu(ii, TJp2)  += bb4*dN_dz[jj];
          }

          bb5 *= (Jg*thetahat);
          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp
    //cout << "elemVolCur = " << elemVolCur << endl;

    /*
    //bool MassLumping=false;
    bool MassLumping=true;

    if(MassLumping)
    {
      fact = rho0*elemVolOrig/10.0;

      if(Mlocal.rows() != nsize)
        Mlocal.resize(nsize, nsize);
      Mlocal.setZero();

      for(ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
    else
    {
      calcMassMatrix(Mlocal, false);
    }

    Flocal1 -= Mlocal*accC;
    Kuu += acceFact1*Mlocal;
    */

    return 0;
}




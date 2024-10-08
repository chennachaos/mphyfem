#include "ElementBase.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "BasisFunctionsLagrange.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "utilitiesmaterialHyperelastic.h"
#include "utilitiesmaterial.h"


extern MyTime           myTime;
extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern bool debug;

using namespace std;



int  ElementBase::calcStiffnessAndResidualMixed2D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalP)
{
    int   err, index, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), Np(nlbfP);

    double  detF, detFn, fact, fact1, dvol, dvol0, Jac, volstrain, resip;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Jhat, thetahat;
    double  param[2], bforce[2]={0.0,0.0}, force[2];
    double  veloCur[2], acceCur[2], sig[2];

    int    sss   = ElemTypeData->getModeltypeNum();
    double thick = ElemTypeData->getThickness();
    bool   axsy  = (sss == 3);
    bool  finite = MatlData->isFiniteStrain();

    int  Utype  = MatlData->getUtype(); // volumetric energy function
    double BULK = 1.0/MatlData->getKinv();

    double rho0 = MatlData->getDensity() ;
    double rho  = rho0 ;

    //bforce[0]   = elmDat[6]*timeFunction[0]->getFactor() ;
    //bforce[1]   = elmDat[7]*timeFunction[0]->getFactor() ;

    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;
    double tCur = myTime.cur;

    double xNode[nlbfU], yNode[nlbfU], geom[2];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    // resize local matrices and initialise them to zero
    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalP.rows() < nlbfP)   FlocalP.resize(nlbfP);
    FlocalP.setZero();
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();
    if(Kup.rows() < nsize || Kup.cols() < nlbfP)   Kup.resize(nsize, nlbfP);
    Kup.setZero();
    if(Kpu.rows() < nlbfP || Kpu.cols() < nsize)   Kpu.resize(nlbfP, nsize);
    Kpu.setZero();
    if(Kpp.rows() < nlbfP || Kpp.cols() < nlbfP)   Kpp.resize(nlbfP, nlbfP);
    Kpp.setZero();

    MatrixXd  Cmat(9,9), Gc(3,9), F(3,3), Fn(3,3);
    VectorXd  stre(9);
    Gc.setZero();
    F.setZero();
    Fn.setZero();
    stre.setZero();


    vector<double>  gausspoints1, gausspoints2, gaussweights;

    nGP = getGaussPoints2D(npElem, gausspoints1, gausspoints2, gaussweights);


    elemVolOrig=0.0;
    elemVolCur=0.0;

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

        geom[0] = geom[1] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += Nu[ii]*xNode[ii];
          geom[1] += Nu[ii]*yNode[ii];
        }

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        computeDefGrad2DPrev(dNu_dx, dNu_dy, Fn);

        detFn = Fn(0,0)*Fn(1,1) - Fn(0,1)*Fn(1,0);

        compute_constants_volumetric_functions(finite, Utype, detFn, BULK, Jhat, thetahat);


        computeDefGrad2DCur(dNu_dx, dNu_dy, F);

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        volstrain = F(0,0) + F(1,1) - 2.0;

        acceCur[0] = computeAccelerationCur(0, Nu);
        acceCur[1] = computeAccelerationCur(1, Nu);

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        adjust_deformationgradient_2D(sss, finite, F);


        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)  dvol *= F(2,2);

        // evaluate pressure at the quadrature points

        if(nlbfP == 1)
        {
          Np[0] = 1.0;
          pres = presDOF[0];
        }
        else
        {
          if(nlbfP == 3)
            LagrangeBasisFunsTria(nlbfP, param[0], param[1], &Np(0));
          else if(nlbfP == 4)
            LagrangeBasisFunsQuad(nlbfP, param[0], param[1], &Np(0));

          pres = 0.0;
          for(ii=0; ii<nlbfP; ii++)
            pres += (Np[ii]*SolnData->presCur[nodeNumsPres[ii]]);
        }

        MatlData->computeStressAndTangent(true, sss, Fn, F, pres, stre, Cmat, ivar, gp, dt);
        //if(err !=0)          return 1;


        // Calculate Stiffness and Residual
        //==============================================

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb4 = dvol*Nu[ii];
          bb5 = dvol0*Nu[ii];

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
            cc1 = dNu_dx[jj];
            cc2 = dNu_dy[jj];
            cc3 = Nu[jj];

            TJ   = 2*jj;
            TJp1 = TJ+1;

            // acceleration term
            acceFact2 = acceFact1*cc3*rho0;

            fact  = bb5*acceFact2;

            // material Stiffness contribution
            //fact  += af*(sig[0]*cc1+sig[1]*cc2)*FiniteFact;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;

            Kuu(TI,   TJ)    +=  af*(Gc(0,0) * cc1 + Gc(0,3) * cc2) ;
            Kuu(TI,   TJp1)  +=  af*(Gc(0,1) * cc1 + Gc(0,4) * cc2) ;
            Kuu(TIp1, TJ)    +=  af*(Gc(1,0) * cc1 + Gc(1,3) * cc2) ;
            Kuu(TIp1, TJp1)  +=  af*(Gc(1,1) * cc1 + Gc(1,4) * cc2) ;
          }

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
          }
        }

        if(finite)
        {
          resip = detF-Jhat-thetahat*pres;
        }
        else
        {
          resip = volstrain-thetahat*pres;
        }

        for(ii=0; ii<nlbfP; ii++)
        {
          bb4 = dvol*Np[ii]*af;
          bb5 = dvol0*Np[ii];

          FlocalP(ii) -= bb5*resip;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Kpu(ii, TJ)    += (bb4*dNu_dx[jj]);
            Kpu(ii, TJp1)  += (bb4*dNu_dy[jj]);
          }

          bb5 *= thetahat*af;
          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

    //printMatrix(Kuu);
    //printMatrix(Kpu);
    //printMatrix(Kup);
    //printMatrix(Kpp);
    //printVector(FlocalU);
    //printVector(FlocalP);

    return 0;
}




int  ElementBase::calcStiffnessAndResidualMixed3D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalP)
{
    int   err = 0,  ii, jj, kk, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  detF, detFn, fact, fact1, dvol, dvol0, Jac, volstrain, Np[nlbfP];
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Jhat, thetahat, resip;
    double  param[3], bforce[3]={0.0,0.0,0.0}, force[3];
    double  veloCur[3], acceCur[3], sig[3];

    MatrixXd  F(3,3), Fn(3,3);
    F.setZero();
    Fn.setZero();

    bool  finite = MatlData->isFiniteStrain();

    int  Utype  = MatlData->getUtype();                     // volumetric energy function
    double BULK = 1.0/MatlData->getKinv();

    double rho0 = MatlData->getDensity();
    double rho  = rho0 ;
    //bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    //bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    //bforce[2]   = elmDat[8]*timeFunction[0].prop ;

    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;
    double tCur = myTime.cur;

    double xNode[npElem], yNode[npElem], zNode[npElem], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    nGP = getGaussPoints3D(npElem, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalP.rows() < nlbfP)   FlocalP.resize(nlbfP);
    FlocalP.setZero();

    if(Kuu.rows() < nsize || Kuu.cols() < nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();
    if(Kup.rows() < nsize || Kup.cols() < nlbfP)   Kup.resize(nsize, nlbfP);
    Kup.setZero();
    if(Kpu.rows() < nlbfP || Kpu.cols() < nsize)   Kpu.resize(nlbfP, nsize);
    Kpu.setZero();
    if(Kpp.rows() < nlbfP || Kpp.cols() < nlbfP)   Kpp.resize(nlbfP, nlbfP);
    Kpp.setZero();


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
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;
        elemVolOrig += dvol0;


        geom[0] = geom[1] = geom[2] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += N[ii]*xNode[ii];
          geom[1] += N[ii]*yNode[ii];
          geom[2] += N[ii]*zNode[ii];
        }

        computeDefGradPrev(dN_dx, dN_dy, dN_dz, Fn);

        detFn = Fn.determinant();

        compute_constants_volumetric_functions(finite, Utype, detFn, BULK, Jhat, thetahat);

        computeDefGradCur(dN_dx, dN_dy, dN_dz, F);

        volstrain = F.trace() - 3.0;

        detF = F.determinant();

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        acceCur[0] = computeAccelerationCur(0, N);
        acceCur[1] = computeAccelerationCur(1, N);
        acceCur[2] = computeAccelerationCur(2, N);

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }
        elemVolCur  += dvol;


        // evaluate pressure at the quadrature points

        if(nlbfP == 1)
        {
          Np[0] = 1.0;
          pres = presDOF[0];
        }
        else
        {
          if(nlbfP == 4)
            LagrangeBasisFunsTetra(nlbfP, param[0], param[1], param[2], Np);
          else if(nlbfP == 8)
            LagrangeBasisFunsHexa(nlbfP, param[0], param[1], param[2], Np);

          pres = 0.0;
          for(ii=0; ii<nlbfP; ii++)
            pres += (Np[ii]*SolnData->presCur[nodeNumsPres[ii]]);
        }


        MatlData->computeStressAndTangent(true, sss, Fn, F, pres, stre, Cmat, ivar, gp, dt);
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

          FlocalU(TI)   += (bb5*force[0] - sig[0]) ;
          FlocalU(TIp1) += (bb5*force[1] - sig[1]) ;
          FlocalU(TIp2) += (bb5*force[2] - sig[2]) ;

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

            // consistent mass matrix
            fact  = bb5*acceFact2;

            // material Stiffness contribution
            //fact = af*(sig[0]*cc1 + sig[1]*cc2 + sig[2]*cc3)*FiniteFact;

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
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
            Kup(TIp2, jj)  += (bb3*Np[jj]);
          }
        }

        if(finite)
        {
          resip = detF-Jhat-thetahat*pres;
        }
        else
        {
          resip = volstrain-thetahat*pres;
        }

        for(ii=0; ii<nlbfP; ii++)
        {
          bb4 = dvol*Np[ii]*af;
          bb5 = dvol0*Np[ii];

          FlocalP(ii) -= bb5*resip;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            Kpu(ii, TJ)    += bb4*dN_dx[jj];
            Kpu(ii, TJp1)  += bb4*dN_dy[jj];
            Kpu(ii, TJp2)  += bb4*dN_dz[jj];
          }

          bb5 *= thetahat*af;
          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }

    }//gp

    //printMatrix(Kuu);
    //printMatrix(Kpu);
    //printMatrix(Kup);
    //printMatrix(Kpp);
    //printVector(FlocalU);
    //printVector(FlocalP);

    return 0;
}

int ElementBase::toComputeInfSupCondition(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
{
  if(ndim == 2)
    ElementBase::toComputeInfSupCondition2D(Kuu, Kup, Kpp);
  else
    ElementBase::toComputeInfSupCondition3D(Kuu, Kup, Kpp);

  return 0;
}




int ElementBase::toComputeInfSupCondition2D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
{
    // to compute the inf-sup constant

    int   err, index, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), Np(nlbfP);

    double  fact, fact1, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[2];

    double xNode[nlbfU], yNode[nlbfU], geom[2];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    // resize local matrices and initialise them to zero
    if(Kuu.rows() != nsize)
    {
      Kuu.resize(nsize, nsize);
      Kup.resize(nsize, nlbfP);
      Kpp.resize(nlbfP, nlbfP);
    }
    Kuu.setZero();
    Kup.setZero();
    Kpp.setZero();


    vector<double>  gausspoints1, gausspoints2, gaussweights;
    nGP = getGaussPoints2D(npElem, gausspoints1, gausspoints2, gaussweights);

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;

        // evaluate pressure at the quadrature points

        if(nlbfP == 1)
        {
          Np[0] = 1.0;
          pres = presDOF[0];
        }
        else
        {
          if(nlbfP == 3)
            LagrangeBasisFunsTria(nlbfP, param[0], param[1], &Np(0));
          else if(nlbfP == 4)
            LagrangeBasisFunsQuad(nlbfP, param[0], param[1], &Np(0));
        }

        // Calculate Stiffness and Residual
        //==============================================

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dNu_dx[jj];
            cc2 = dNu_dy[jj];
            cc3 = Nu[jj];

            TJ   = 2*jj;
            TJp1 = TJ+1;

            fact = bb1*dNu_dx[jj] + bb2*dNu_dy[jj] ;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;
          }

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
          }
        }

        for(ii=0; ii<nlbfP; ii++)
        {
          bb4 = dvol*Np[ii];

          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)  += bb4*Np[jj];
          }
        }
    }//gp

    return 0;
}




int ElementBase::toComputeInfSupCondition3D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
{
    // to compute the inf-sup constant

    int   err, index, ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), dNu_dz(nlbfU), Np(nlbfP);

    double  fact, fact1, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[3];

    double xNode[nlbfU], yNode[nlbfU], zNode[nlbfU];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    // resize local matrices and initialise them to zero
    if(Kuu.rows() != nsize)
    {
      Kuu.resize(nsize, nsize);
      Kup.resize(nsize, nlbfP);
      Kpp.resize(nlbfP, nlbfP);
    }
    Kuu.setZero();
    Kup.setZero();
    Kpp.setZero();


    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    nGP = getGaussPoints3D(npElem, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), &dNu_dz(0), Jac);

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;

        // evaluate pressure at the quadrature points

        if(nlbfP == 1)
        {
          Np[0] = 1.0;
          pres = presDOF[0];
        }
        else
        {
          if(nlbfP == 4)
            LagrangeBasisFunsTetra(nlbfP, param[0], param[1], param[2], &Np(0));
          else if(nlbfP == 8)
            LagrangeBasisFunsHexa(nlbfP, param[0], param[1], param[2], &Np(0));
        }

        // Calculate Stiffness and Residual
        //==============================================

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb3 = dvol*dNu_dz[ii];

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          for(jj=0; jj<nlbfU; jj++)
          {
            cc1 = dNu_dx[jj];
            cc2 = dNu_dy[jj];
            cc3 = dNu_dz[jj];
            cc4 = Nu[jj];

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            fact = bb1*dNu_dx[jj] + bb2*dNu_dy[jj] + bb3*dNu_dz[jj] ;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;
            Kuu(TIp2, TJp2)  += fact ;
          }

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
            Kup(TIp2, jj)  += (bb3*Np[jj]);
          }
        }

        for(ii=0; ii<nlbfP; ii++)
        {
          bb4 = dvol*Np[ii];

          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)  += bb4*Np[jj];
          }
        }
    }//gp

    return 0;
}









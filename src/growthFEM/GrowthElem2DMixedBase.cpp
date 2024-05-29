
#include "GrowthElem2DMixedBase.h"
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



GrowthElem2DMixedBase::GrowthElem2DMixedBase()
{
}



GrowthElem2DMixedBase::~GrowthElem2DMixedBase()
{
}


int GrowthElem2DMixedBase::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, detFn, Jne, fact, fact1, dvol, dvol0, Jac, volstrain, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Rp, Jhat, thetahat;
    double  param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2], Np[nlbfP];

    double  timeFactor = timeFunction[0].prop, timeFactorPrev = timeFactor-myTime.dt, Jg;
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

    double xNode[nlbfU], yNode[nlbfU], geom[3];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
      getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);
    else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
      getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);
    else
    {
        throw runtime_error("Element type no defined...");
    }

    // resize local matrices and initialise them to zero
    if(Flocal1.rows() < nsize)   Flocal1.resize(nsize);
    Flocal1.setZero();
    if(Flocal2.rows() < nlbfP)   Flocal2.resize(nlbfP);
    Flocal2.setZero();
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();
    if(Kup.rows() < nsize || Kup.cols() < nlbfP)   Kup.resize(nsize, nlbfP);
    Kup.setZero();
    if(Kpu.rows() < nlbfP || Kpu.cols() < nsize)   Kpu.resize(nlbfP, nsize);
    Kpu.setZero();
    if(Kpp.rows() < nlbfP || Kpp.cols() < nlbfP)   Kpp.resize(nlbfP, nlbfP);
    Kpp.setZero();

    MatrixXd  Cmat(9,9), Gc(3,9), F(3,3), Fn(3,3), Fe(3,3), Fg(3,3), FgInv(3,3);
    VectorXd  stre(9);
    Gc.setZero();
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    stre.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    //cout << elmType << '\t' <<  matType << '\t' <<  secType << '\t' << growthFactor << endl;
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

        growthfunction_beam_singlelayer(timeFactorPrev,  geom,  Fg);

        Jg = Fg.determinant();
        Jne = detFn / Jg;

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

        //growthfunction_beam_singlelayer(timeFactor,  geom,  Fg);
        //growthfunction_beam_multilayer(timeFactor, geom,  Fg);
        growthfunction_isotropic(1, ndim,  growthFactor,  Fg);
        FgInv = Fg.inverse();
        Jg = Fg.determinant();

        // elastic part of the deformation gradient
        Fe = F*FgInv;


        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        //if(finite)  dvol *= F(2,2);

        // evaluate pressure at the quadrature points

        if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
            BernsteinBasisFunsQuad(1, param[0], param[1], Np);
        else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
            BernsteinBasisFunsTria(1, param[0], param[1], Np);

        pres = 0.0;
        for(ii=0; ii<nlbfP; ii++)
          //pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);
          pres += (Np[ii]*SolnData->var2[forAssyVecPres[ii]]);

        MatlData->computeStressAndTangent(true, sss, Fn, Fe, pres, stre, Cmat);
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

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
          }
        }

        if(finite)
        {
          Rp = detF - Jg*(Jhat+thetahat*pres);
        }
        else
        {
          Rp = volstrain-pres*eps;
        }

        for(ii=0; ii<nlbfP; ii++)
        {
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

          bb5 *= (Jg*thetahat);
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
    //printVector(Flocal1);
    //printVector(Flocal2);

    return 0;
}






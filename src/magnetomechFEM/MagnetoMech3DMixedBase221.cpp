
#include "MagnetoMech3DMixedBase221.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"


using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



MagnetoMech3DMixedBase221::MagnetoMech3DMixedBase221()
{
}



MagnetoMech3DMixedBase221::~MagnetoMech3DMixedBase221()
{
}




int MagnetoMech3DMixedBase221::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != nsize) || (Mlocal.cols() != nsize) )
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = elmDat[5]*computeVolume(true)/nlbfU;

      for(int ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
    else
    {
      int  nGPt, ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;
      double  fact, dvol0, Jac, bb1, cc1, param[3];
      double rho0 = MatlData->getDensity();

      VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

      vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
      if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
      {
        nGPt=10;
        getGaussPointsTetra(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
      }
      else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
      {
        nGPt=9;
        getGaussPointsWedge(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
      }
      else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
      {
        nGPt=27;
        getGaussPointsHexa(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
      }
      else
      {
        throw runtime_error("Element type not defined...");
      }

      for(gp=0; gp<nGPt; gp++)
      {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];

          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol0 = gaussweights[gp]*Jac;

          for(ii=0;ii<nlbfU;ii++)
          {
            bb1 = (dvol0*rho0)*N[ii];

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






int  MagnetoMech3DMixedBase221::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalF, VectorXd& FlocalP, bool firstIter)
{
    int   err = 0, ii, jj, kk, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU), Np(nlbfP);

    double  detF, detFn, trF, fact, fact1, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Rp, Jhat, thetahat;
    double  param[3], bforce[3], force[3], tarr[6],  tarr2[6];
    double  veloCur[3], acceCur[3], sig[3];

    MatrixXd  Cmat(9,9), Bmat(9,3), Amat(3,3),  Gc(3,9), F(3,3), Fn(3,3);
    Amat.setZero();
    Bmat.setZero();
    Cmat.setZero();
    Gc.setZero();
    VectorXd  stre(9),  magnDisp(3),  magnField(3);
    F.setZero();
    stre.setZero();

    double  timeFactor = timeFunction[0].prop, timeFactorPrev = timeFactor-myTime.dt;

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

    double xNode[27], yNode[27], zNode[27], xx, yy,  zz;
    for(ii=0;ii<npElem;ii++)
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
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)     Kuu.resize(nsize, nsize);
    Kuu.setZero();

    if(Kuf.rows() < nsize || Kuf.cols() < nlbfF)     Kuf.resize(nsize, nlbfF);
    Kuf.setZero();
    if(Kfu.rows() < nlbfF || Kfu.cols() < nsize)     Kfu.resize(nlbfF, nsize);
    Kfu.setZero();
    if(Kff.rows() < nlbfF || Kff.cols() < nlbfF)     Kff.resize(nlbfF, nlbfF);
    Kff.setZero();

    if(Kup.rows() < nsize || Kup.cols() < nlbfP)     Kup.resize(nsize, nlbfP);
    Kup.setZero();
    if(Kpu.rows() < nlbfP || Kpu.cols() < nsize)     Kpu.resize(nlbfP, nsize);
    Kpu.setZero();
    if(Kpp.rows() < nlbfP || Kpp.cols() < nlbfP)     Kpp.resize(nlbfP, nlbfP);
    Kpp.setZero();

    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalF.rows() < nlbfF)   FlocalF.resize(nlbfF);
    FlocalF.setZero();
    if(FlocalP.rows() < nlbfP)   FlocalP.resize(nlbfP);
    FlocalP.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;


    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = gaussweights[gp]*Jac;
        dvol  = dvol0;
        elemVolOrig += dvol0;

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

        acceCur[0] = computeValueDotDotCur(0, N);
        acceCur[1] = computeValueDotDotCur(1, N);
        acceCur[2] = computeValueDotDotCur(2, N);

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }
        elemVolCur  += dvol;


        // evaluate phi at the quadrature points
        magnField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          magnField[0] -= (dN_dx[ii]*SolnData->var3Cur[nodeNums[ii]]);
          magnField[1] -= (dN_dy[ii]*SolnData->var3Cur[nodeNums[ii]]);
          magnField[2] -= (dN_dz[ii]*SolnData->var3Cur[nodeNums[ii]]);
        }
        //printVector(SolnData->var3Cur);
        //printVector(dN_dy);
        //printVector(magnField);

        // evaluate pressure at the quadrature point

        if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
          BernsteinBasisFunsTetra(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
          BernsteinBasisFunsWedge(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
          BernsteinBasisFunsHexa(1, param[0], param[1], param[2], &Np(0));

        pres = computeValue2(0, Np);

        MatlData->computeStressAndTangent(true, sss, Fn, F, magnField, pres, stre, magnDisp, Cmat, Bmat, Amat, ivar, gp, dt);
        if(err !=0)          return 1;

        //printMatrix(Bmat);


        // Calculate Stiffness and Residual
        //==============================================

        xx = yy = zz = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
          zz += N[ii]*zNode[ii];
        }

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];
        force[2] = rho0*bforce[2];

        // contribution for the acceleration is added after the GP loop
        //force[0] -= rho0*acceCur[0];
        //force[1] -= rho0*acceCur[1];
        //force[2] -= rho0*acceCur[2];

        for(ii=0; ii<nlbfU; ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*dN_dz[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

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
            cc1 = dN_dx[jj];
            cc2 = dN_dy[jj];
            cc3 = dN_dz[jj];
            cc4 = N[jj];

            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            // acceleration term
            acceFact2 = acceFact1*cc4*rho0;

            fact  = bb5*acceFact2;
            // contribution for the acceleration is added after the GP loop
            fact  = 0.0;

            // material Stiffness contribution
            //fact  += af*(sig[0]*cc1+sig[1]*cc2+sig[2]*cc3)*FiniteFact;

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
          Gc(0,0) = bb1*Bmat(0,0) + bb2*Bmat(3,0) + bb3*Bmat(6,0);
          Gc(0,1) = bb1*Bmat(0,1) + bb2*Bmat(3,1) + bb3*Bmat(6,1);
          Gc(0,2) = bb1*Bmat(0,2) + bb2*Bmat(3,2) + bb3*Bmat(6,2);

          Gc(1,0) = bb1*Bmat(1,0) + bb2*Bmat(4,0) + bb3*Bmat(7,0);
          Gc(1,1) = bb1*Bmat(1,1) + bb2*Bmat(4,1) + bb3*Bmat(7,1);
          Gc(1,2) = bb1*Bmat(1,2) + bb2*Bmat(4,2) + bb3*Bmat(7,2);

          Gc(2,0) = bb1*Bmat(2,0) + bb2*Bmat(5,0) + bb3*Bmat(8,0);
          Gc(2,1) = bb1*Bmat(2,1) + bb2*Bmat(5,1) + bb3*Bmat(8,1);
          Gc(2,2) = bb1*Bmat(2,2) + bb2*Bmat(5,2) + bb3*Bmat(8,2);

          for(jj=0; jj<nlbfF; jj++)
          {
            tarr2[0] = af*( Gc(0,0)*dN_dx[jj] + Gc(0,1)*dN_dy[jj] + Gc(0,2)*dN_dz[jj] );
            tarr2[1] = af*( Gc(1,0)*dN_dx[jj] + Gc(1,1)*dN_dy[jj] + Gc(1,2)*dN_dz[jj] );
            tarr2[2] = af*( Gc(2,0)*dN_dx[jj] + Gc(2,1)*dN_dy[jj] + Gc(2,2)*dN_dz[jj] );

            Kuf(TI,   jj)  += tarr2[0];
            Kuf(TIp1, jj)  += tarr2[1];
            Kuf(TIp2, jj)  += tarr2[2];

            Kfu(jj, TI)    += tarr2[0];
            Kfu(jj, TIp1)  += tarr2[1];
            Kfu(jj, TIp2)  += tarr2[2];
          }

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*af*Np[jj]);
            Kup(TIp1, jj)  += (bb2*af*Np[jj]);
            Kup(TIp2, jj)  += (bb3*af*Np[jj]);
          }
        }

        //  Kff matrix
        for(ii=0; ii<nlbfF; ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*dN_dz[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          FlocalF(ii) -= (bb1*magnDisp(0)+bb2*magnDisp(1)+bb3*magnDisp(2));

          tarr[0] = af*( bb1*Amat(0,0) + bb2*Amat(1,0) + bb3*Amat(2,0) );
          tarr[1] = af*( bb1*Amat(0,1) + bb2*Amat(1,1) + bb3*Amat(2,1) );
          tarr[2] = af*( bb1*Amat(0,2) + bb2*Amat(1,2) + bb3*Amat(2,2) );

          for(jj=0; jj<nlbfF; jj++)
          {
            Kff(ii, jj) -= (dN_dx[jj]*tarr[0] + dN_dy[jj]*tarr[1] + dN_dz[jj]*tarr[2]);
          }
        }


        if(finite)
        {
          Rp = detF - (Jhat+thetahat*pres);
        }
        else
        {
          Rp = volstrain-thetahat*pres;
        }

        for(ii=0; ii<nlbfP; ii++)
        {
          bb4 = dvol*Np[ii]*af;
          bb5 = dvol0*Np[ii];

          FlocalP(ii) -= bb5*Rp;

          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 3*jj;
            TJp1 = TJ+1;
            TJp2 = TJ+2;

            Kpu(ii, TJ)    += (bb4*dN_dx[jj]);
            Kpu(ii, TJp1)  += (bb4*dN_dy[jj]);
            Kpu(ii, TJp2)  += (bb4*dN_dz[jj]);
          }

          bb5 *= thetahat;
          for(jj=0; jj<nlbfP; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

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

    if(SolnData->tis > 0)
    {
        VectorXd  accC(nsize);

        for(ii=0;ii<npElem;ii++)
        {
          int ind1 = ii*ndim;
          int ind2 = nodeNums[ii]*ndim;

          for(kk=0;kk<ndim;kk++)
          {
            accC(ind1+kk)   =  SolnData->var1DotDotCur[ind2+kk];
          }
        }

        Cmat.resize(nsize, nsize);

        calcMassMatrix(Cmat, false);

        Kuu     +=  acceFact1*Cmat;
        FlocalU -=  Cmat*accC;
    }

    return 0;
}




/*
int MagnetoMech3DMixedBase221::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  detF, detFn, Jne, fact, fact1, dvol, dvol0, Jac, volstrain, r1d3=1.0/3.0, r2d3=2.0/3.0;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Rp, dUdJ, d2UdJ2, Jhat, thetahat;
    double  param[3], bforce[3], force[3];
    double  veloCur[3], acceCur[3], sig[3], Np[nlbfP];

    MatrixXd  Fg(3,3), FgInv(3,3);
    double  timeFactor = elmDat[10]*timeFunction[0].prop, timeFactorPrev = timeFactor-myTime.dt, Jg ;

    double mu   = matDat[0];      // shear modulus
    double eps  = matDat[1];      // inverse of bulk modulus
    double BULK = 1.0/eps;
    int  Utype  = int(matDat[2]); // volumetric energy function

    double rho0 = elmDat[5];
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
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();

    MatrixXd  Cmat(9,9), Gc(3,9), F(3,3), Fe(3,3);
    VectorXd  stre(9);
    Gc.setZero();
    F.setZero();
    stre.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

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

        F(0,0) = computeValuePrev(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValuePrev(0, dN_dy);                  // eps_xy
        F(0,2) = computeValuePrev(0, dN_dz);                  // eps_xz
        F(1,0) = computeValuePrev(1, dN_dx);                  // eps_yx
        F(1,1) = computeValuePrev(1, dN_dy) + 1.0;            // eps_yy
        F(1,2) = computeValuePrev(1, dN_dz);                  // eps_yz
        F(2,0) = computeValuePrev(2, dN_dx);                  // eps_zx
        F(2,1) = computeValuePrev(2, dN_dy);                  // eps_zy
        F(2,2) = computeValuePrev(2, dN_dz) + 1.0;            // eps_zz

        detFn = F.determinant();

        growthfunction_beam_singlelayer(timeFactorPrev,  geom,  Fg);

        Jg = Fg.determinant();
        Jne = detFn / Jg;


        if(finite)
        {
          if(Utype == -1)
          {
            dUdJ   = 1.0;
            d2UdJ2 = 1.0;
            Jhat   = 1.0;
            thetahat = 0.0;
          }
          else
          {
            dUdJ   = volumetricfunctions_getFirstDerivative(Utype, BULK, Jne);
            d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, BULK, Jne);

            Jhat = Jne - dUdJ/d2UdJ2;
            thetahat = 1.0/d2UdJ2;
          }
        }
        else
        {
          thetahat = 1.0/BULK;
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

        volstrain = F(0,0)+F(1,1)+F(2,2) - 3.0;

        detF = F.determinant();

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        //growthfunction_beam_singlelayer(timeFactor,  geom,  Fg);
        //growthfunction_beam_singlelayer3D(timeFactor,  geom,  Fg);
        //growthfunction_circularplate_singlelayer(timeFactor, geom, Fg);
        //if(geom[2] > 0.2) fact = timeFactor;
        //else fact = 0.0;
        //fact = (geom[0]*sin(2.0*PI*geom[2])+geom[1]*cos(2.0*PI*geom[2]))*timeFactor;
        //fact = timeFactor*geom[2];
        fact = timeFactor;
        growthfunction_isotropic(3, ndim, fact,  Fg);
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
          BernsteinBasisFunsTetra(1, param[0], param[1], param[2], Np);
        else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
          BernsteinBasisFunsWedge(1, param[0], param[1], param[2], Np);
        else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
          BernsteinBasisFunsHexa(1, param[0], param[1], param[2], Np);


        // evaluate pressure at the quadrature points
        pres = 0.0; presPrev = 0.0;
        for(ii=0; ii<nlbfP; ii++)
        {
          //pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);
          //presPrev += (Np[ii]*SolnData->var2Prev[nodeNums[ii]]);
          pres += (Np[ii]*SolnData->var2[forAssyVecPres[ii]]);
        }
        //pres = af*pres + (1.0-af)*presPrev;

        MatlData->computeStressAndTangent(true, sss, Fe, pres, stre, Cmat);
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
          //Rp = detF/Jg - 1.0 - pres*Kinv;

          Rp = detF - Jg*(Jhat+thetahat*pres);
        }
        else
        {
          Rp = volstrain - pres*eps;
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


    return 0;
}
*/







void MagnetoMech3DMixedBase221::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'MagnetoMech3DMixedBase221::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
      vals2project[0] += outval[ii];

    vals2project[0] /= nGP;

    return;
}


void MagnetoMech3DMixedBase221::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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
    if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
      stressrecovery_extrapolate_Tetrahedron(degree, outval, &vals2project[0]);
    else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
      stressrecovery_extrapolate_Wedge(degree, outval, &vals2project[0]);
    else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
      stressrecovery_extrapolate_Hexahedron(degree, outval, &vals2project[0]);

    //printVector(vals2project);

    return;
}


void MagnetoMech3DMixedBase221::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    double  detF, dvol, Jac, fact, dvol0, pbar, param[3];
    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), Bmat(9,3), Amat(3,3);
    VectorXd  stre(9),  magnDisp(3),  magnField(3);
    stre.setZero();
    F.setZero();
    Fn.setZero();

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU), Np(nlbfP);

    int   err = 0;

    double dt = myTime.dt;

    double xNode[npElem], yNode[npElem], zNode[npElem], xx, yy, zz;
    for(int ii=0;ii<npElem;ii++)
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


    for(int gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE_HEXA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        computeDefGradPrev(dN_dx, dN_dy, dN_dz, Fn);
        computeDefGrad(dN_dx, dN_dy, dN_dz, F);


        // evaluate pressure at the quadrature point
        if(ELEM_SHAPE == ELEM_SHAPE_TETRA_BERNSTEIN)
          BernsteinBasisFunsTetra(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_WEDGE_BERNSTEIN)
          BernsteinBasisFunsWedge(1, param[0], param[1], param[2], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_HEXA_BERNSTEIN)
          BernsteinBasisFunsHexa(1, param[0], param[1], param[2], &Np(0));

        pres = computeValue2(0, Np);

        if(varindex != 10)
          MatlData->computeStressAndTangent(false, sss, Fn, F, magnField, pres, stre, magnDisp, Cmat, Bmat, Amat, ivar, gp, dt);

        if(varindex < 9)
           outval[gp] = stre[varindex];
        else if(varindex == 9)
           outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*(stre[1]*stre[1]+stre[2]*stre[2]+stre[5]*stre[5]))/2.0);
        else if(varindex == 10)
          //outval[gp] = (stre[0]+stre[4]+stre[8])/3.0;
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



void MagnetoMech3DMixedBase221::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void MagnetoMech3DMixedBase221::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}



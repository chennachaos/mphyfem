
#include "MagnetoMech2DMixedBase221.h"
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



MagnetoMech2DMixedBase221::MagnetoMech2DMixedBase221()
{
}



MagnetoMech2DMixedBase221::~MagnetoMech2DMixedBase221()
{
}


int  MagnetoMech2DMixedBase221::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalF, VectorXd& FlocalP, bool firstIter)
{
    int  err = 0,  ii, jj, kk, gp, TI, TIp1, TJ, TJp1;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU);

    double  detF, detFn, trF, fact, fact1, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Rp, Jhat, thetahat;
    double  param[3], bforce[3], force[3], tarr[6],  tarr2[6];
    double  veloCur[3], acceCur[3], sig[3], Np[nlbfP];
    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), Bmat(9,3), Amat(3,3),  Gc(3,9);
    Amat.setZero();
    Bmat.setZero();
    Cmat.setZero();
    Gc.setZero();
    F.setZero();
    Fn.setZero();

    VectorXd  stre(9),  elecDisp(3),  elecField(3);
    stre.setZero();

    int  Utype  = int(matDat[0]); // volumetric energy function
    double eps  = matDat[1];      // inverse of bulk modulus
    double BULK = 1.0/eps;

    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

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
        throw runtime_error("Element type not defined...");
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

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_QUAD_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

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

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        adjust_deformationgradient_2D(sss, finite, F);


        acceCur[0] = computeValueDotDotCur(0, Nu);
        acceCur[1] = computeValueDotDotCur(1, Nu);

        //printVector(dNu_dx);
        //printVector(dNu_dy);
        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_QUAD_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)    dvol *= F(2,2);

        // evaluate phi at the quadrature point

        elecField.setZero();
        for(ii=0; ii<nlbfF; ii++)
        {
          elecField[0] -= (dNu_dx[ii]*SolnData->var3[nodeNums[ii]]);
          elecField[1] -= (dNu_dy[ii]*SolnData->var3[nodeNums[ii]]);
        }

        // evaluate pressure at the quadrature point
        if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
            BernsteinBasisFunsQuad(1, param[0], param[1], Np);
        else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
            BernsteinBasisFunsTria(1, param[0], param[1], Np);

        pres = 0.0;
        for(ii=0; ii<nlbfP; ii++)
        {
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);
        }

        MatlData->computeStressAndTangent(true, sss, Fn, F, elecField, pres, stre, elecDisp, Cmat, Bmat, Amat, ivar, gp, dt);
        if(err !=0)          return 1;

        // Calculate Stiffness and Residual
        //==============================================

        // body force and acceleration terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        for(ii=0; ii<nlbfU; ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb4 = dvol*Nu[ii];
          bb5 = dvol0*Nu[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          //  Kuu matrix
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
            cc4 = Nu[jj];

            TJ   = 2*jj;
            TJp1 = TJ+1;

            // acceleration term
            acceFact2 = acceFact1*cc4*rho0;

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
          //printMatrix(Kuu);

          //  Kuf & Kfu matrices
          Gc(0,0) = bb1 * Bmat(0,0) + bb2 * Bmat(3,0);
          Gc(0,1) = bb1 * Bmat(0,1) + bb2 * Bmat(3,1);

          Gc(1,0) = bb1 * Bmat(1,0) + bb2 * Bmat(4,0);
          Gc(1,1) = bb1 * Bmat(1,1) + bb2 * Bmat(4,1);

          for(jj=0; jj<nlbfF; jj++)
          {
            tarr2[0] = Gc(0,0) * dNu_dx[jj] + Gc(0,1) * dNu_dy[jj];
            tarr2[1] = Gc(1,0) * dNu_dx[jj] + Gc(1,1) * dNu_dy[jj];

            Kuf(TI,   jj)  += tarr2[0];
            Kuf(TIp1, jj)  += tarr2[1];

            Kfu(jj, TI)    += tarr2[0];
            Kfu(jj, TIp1)  += tarr2[1];
          }

          for(jj=0; jj<nlbfP; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
          }
        }

        //  Kff matrix
        for(ii=0; ii<nlbfF; ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb4 = dvol*Nu[ii];
          bb5 = dvol0*Nu[ii];

          FlocalF(ii) -= (bb1*elecDisp(0)+bb2*elecDisp(1));

          tarr[0] = bb1 * Amat(0,0) + bb2 * Amat(1,0);
          tarr[1] = bb1 * Amat(0,1) + bb2 * Amat(1,1);

          for(jj=0; jj<nlbfF; jj++)
          {
            Kff(ii, jj) -= (dNu_dx[jj]*tarr[0] + dNu_dy[jj]*tarr[1]);
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
            TJ   = 2*jj;
            TJp1 = TJ+1;

            Kpu(ii, TJ)    += (bb4*dNu_dx[jj]);
            Kpu(ii, TJp1)  += (bb4*dNu_dy[jj]);
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

    //printMatrix(Kuf);
    //printMatrix(Kfu);
    //printMatrix(Kff);

    printMatrix(Kup);
    printMatrix(Kpu);
    printMatrix(Kpp);

    printVector(FlocalU);
    printVector(FlocalF);
    //printVector(FlocalP);
    */

    return 0;
}








void MagnetoMech2DMixedBase221::elementContourplot(int vartype, int varindex, int index)
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


void MagnetoMech2DMixedBase221::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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
    if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
        stressrecovery_extrapolate_Quadrilateral(degree, outval, &vals2project[0]);
    else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
        stressrecovery_extrapolate_Triangle(degree, outval, &vals2project[0]);

    //printVector(vals2project);

    return;
}


void MagnetoMech2DMixedBase221::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    int err = 0, ii, jj, gp;

    double  detF, dvol, Jac, fact, dvol0, param[3], Np[nlbfP];
    MatrixXd  F(3,3), Fn(3,3), Cmat(9,9), Bmat(3,9), Amat(3,3);
    VectorXd  stre(9),  elecDisp(3),  elecField(3);
    F.setZero();
    stre.setZero();

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU);

    double dt   = myTime.dt;

    double xNode[npElem], yNode[npElem], xx, yy;
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


    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_QUAD_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

        computeDefGrad2D(dNu_dx, dNu_dy, F);

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        adjust_deformationgradient_2D(sss, finite, F);


        // basis functions for pressure dof
        if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
            BernsteinBasisFunsQuad(1, param[0], param[1], Np);
        else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
            BernsteinBasisFunsTria(1, param[0], param[1], Np);

        pres = 0.0;
        for(ii=0; ii<nlbfP; ii++)
        {
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);
        }

        MatlData->computeStressAndTangent(false, sss, Fn, F, elecField, pres, stre, elecDisp, Cmat, Bmat, Amat, ivar, gp, dt);

        if(varindex < 9)
          outval[gp] = stre[varindex];
        else if(varindex == 9)
          outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*stre[1]*stre[1])/2.0);
        else if(varindex == 10)
          //outval[gp] = (stre[0]+stre[4]+stre[8])/3.0;
          outval[gp] = pres;
        else
        {
          VectorXcd eivals = F.eigenvalues();
          outval[gp] = eivals[0].real();
        }
    }//gp

    return;
}



void MagnetoMech2DMixedBase221::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void MagnetoMech2DMixedBase221::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}


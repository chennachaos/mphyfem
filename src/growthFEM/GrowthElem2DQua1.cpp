
/* Node ordering for the bilinear quadrilateral element
 4-----------3
 |           |
 |           |
 |           |
 |           |
 |           |
 1-----------2
*/

#include "GrowthElem2DQua1.h"
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


GrowthElem2DQua1::GrowthElem2DQua1()
{
  ELEM_SHAPE = ELEM_SHAPE_QUAD_BERNSTEIN;

  ndim   = 2;
  degree = 1;
  npElem = 4;
  nlbfU  = 4;
  nlbfP  = 0;
  ndof   = 2;
  nsize  = npElem*ndof;
}


GrowthElem2DQua1::~GrowthElem2DQua1()
{
}


void GrowthElem2DQua1::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double GrowthElem2DQua1::computeVolume(bool init)
{
    int  ii, jj, gp;
    double  dvol, Jac, param[2];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double xNode[nlbfU], yNode[nlbfU], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    elemVol = 0.0;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        if(init)
          GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
        else
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

         xx = yy= 0.0;
         for(ii=0;ii<nlbfU;ii++)
         {
           xx += N[ii]*xNode[ii];
           yy += N[ii]*yNode[ii];
         }

         dvol = gaussweights[gp]*(Jac*thick);

         if(axsy)
           dvol *= 2.0*PI*xx;

         elemVol += dvol;
    }

    return  elemVol;
}



int GrowthElem2DQua1::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
  // mass lumping - row-wise sum
  if(MassLumping)
  {
      double rho0  = elmDat[5] ;
      double fact  = rho0*computeVolume(true)/4.0;

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

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  fact, dvol, Jac, bb1, cc1, param[3];
    double rho0  = elmDat[5] ;

    double xNode[nlbfU], yNode[nlbfU], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Mlocal.rows() != nsize)
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int nGPt=9;
    getGaussPointsQuad(nGPt, gausspoints1, gausspoints2, gaussweights);

    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        xx = yy= 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        dvol = gaussweights[gp]*(Jac*thick);

        if(axsy)
          dvol *= 2.0*PI*xx;

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = (dvol*rho0)*N[ii];

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
    //cout << " elemVol = " << elemVol << endl;
    //printMatrix(Klocal);  printf("\n\n\n");
  }

  return 0;
}


double  GrowthElem2DQua1::calcCriticalTimeStep(bool flag)
{
  double  dtCric;

  return  dtCric;
}



int GrowthElem2DQua1::calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    MatrixXd  Fg(3,3), FgInv(3,3), F(3,3), Fn(3,3), Fe(3,3);
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    Fg.setZero();
    FgInv.setZero();

    double  timeFactor = timeFunction[0].prop, Jg;
    double  growthFactor = elmDat[10]*timeFunction[0].prop ;

    double rho  = elmDat[5] ;
    double rho0 = rho ;
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

    if(Klocal.rows() != nsize)
    {
      Klocal.resize(nsize, nsize);
      Flocal.resize(nsize);
    }
    Klocal.setZero();
    Flocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    MatrixXd  Cmat(9,9), Gc(3,9);
    VectorXd  stre(9);
    Gc.setZero();
    stre.setZero();


    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(thick*Jac);

        geom[0] = geom[1] = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          geom[0] += N[ii]*xNode[ii];
          geom[1] += N[ii]*yNode[ii];
        }

        if(axsy)
          dvol0 *= 2.0*PI*geom[1];

        dvol  = dvol0;

        F(0,0) = computeValueCur(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dN_dy);                  // eps_xy
        F(1,0) = computeValueCur(1, dN_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

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

        growthfunction_beam_singlelayer(timeFactor,  geom,  Fg);
        //growthfunction_beam_multilayer(timeFactor, geom,  Fg);
        //growthfunction_isotropic(1, ndim,  growthFactor,  Fg);
        FgInv = Fg.inverse();
        Jg = Fg.determinant();
        Jg = 1.0;

        // elastic part of the deformation gradient
        Fe = F*FgInv;
        Fe = F;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(thick*Jac);
        }

        pres = 0.0;
        MatlData->computeStressAndTangent(true, sss, Fn, Fe, pres, stre, Cmat);
        count++;
        ll += nivGP;


        if(finite)  dvol *= F(2,2);

        // Calculate Stiffness and Residual
        //==============================================

        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        //   part 1. -- material part (not necessarily symmetric!!)
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

          Flocal(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal(TIp1) += (bb5*force[1] - sig[1]) ;

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

              Klocal(TI,   TJ)    += fact ;
              Klocal(TIp1, TJp1)  += fact ;

              Klocal(TI,   TJ)    +=  af*(Gc(0,0) * cc1 + Gc(0,3) * cc2) ;
              Klocal(TI,   TJp1)  +=  af*(Gc(0,1) * cc1 + Gc(0,4) * cc2) ;
              Klocal(TIp1, TJ)    +=  af*(Gc(1,0) * cc1 + Gc(1,3) * cc2) ;
              Klocal(TIp1, TJp1)  +=  af*(Gc(1,1) * cc1 + Gc(1,4) * cc2) ;
          }
        }
    } //gp

    //printMatrix(Klocal);  printf("\n\n\n");
    //printVector(Flocal);  printf("\n\n\n");

    return 0;
}





int GrowthElem2DQua1::calcForceVectorGrowthModel(VectorXd& Flocal1, VectorXd& Flocal2)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5;
    double  param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    double  gp1 = 1.0 + elmDat[10]*timeFunction[0].prop ;
    double  Jg  = gp1*gp1;
    MatrixXd  FgInv(3,3), F(3,3), Fn(3,3), Fe(3,3);
    F.setZero();
    Fn.setZero();
    Fe.setZero();
    FgInv.setZero();


    FgInv(0,0) = 1.0/gp1;    FgInv(0,1) = 0.0;        FgInv(0,2) = 0.0;
    FgInv(1,0) = 0.0;        FgInv(1,1) = 1.0/gp1;    FgInv(1,2) = 0.0;
    FgInv(2,0) = 0.0;        FgInv(2,1) = 0.0;        FgInv(2,2) = 1.0;

    double  FgH = -2.0/gp1;

    double rho  = elmDat[5] ;
    double rho0 = rho ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

    double xNode[nlbfU], yNode[nlbfU], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Flocal1.rows() != nsize)
    {
      Flocal1.resize(nsize);
    }
    Flocal1.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    MatrixXd  Cmat(9,9), Gc(3,9);
    VectorXd  stre(9);
    Gc.setZero();
    stre.setZero();


    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
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

        dvol  = dvol0;

        F(0,0) = computeValueCur(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValueCur(0, dN_dy);                  // eps_xy
        F(1,0) = computeValueCur(1, dN_dx);                  // eps_yx
        F(1,1) = computeValueCur(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

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

        // elastic part of the deformation gradient
        Fe = F*FgInv;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(thick*Jac);
        }

        pres = 0.0;
        MatlData->computeStressAndTangent(true, sss, Fn, Fe, pres, stre, Cmat);
        count++;
        ll += nivGP;

        Cmat /=  Jg;
        stre /=  Jg;


        if(finite)  dvol *= F(2,2);

        // multiply with the factor corresponding to the derivative of F wrt to growth factor
        stre *= FgH;

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb4 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[1] + bb2*stre[4];

          Flocal1(TI)   += (sig[0]) ;
          Flocal1(TIp1) += (sig[1]) ;
        }
    } //gp

    //printVector(Flocal);  printf("\n\n\n");

    return 0;
}



int GrowthElem2DQua1::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  F[4], detF, F33, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[4], cc[4][4], param[2], bforce[2], force[2];

    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;

    double xNode[nlbfU], yNode[nlbfU], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
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
          dvol0 *= 2.0*PI*xx;

        dvol  = dvol0;

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[2] = computeValueCur(0, dN_dy);
        F[1] = computeValueCur(1, dN_dx);
        F[3] = computeValueCur(1, dN_dy) + 1.0;

        detF =  F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(thick*Jac);
        }

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

        //printf("F... \t%20.18f\t%20.18f\t%20.18f\t%20.18f\t%20.18f \n", F[0], F[1], F[2], F[3], detF);

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;


        if(err !=0)    return 1;

        if(finite)  dvol *= F33;

        // Calculate Residual
        //==============================================

        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dN_dx[ii];
          bb2 = dvol*dN_dy[ii];
          bb3 = dvol*N[ii];
          bb5 = dvol0*N[ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          Flocal(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3]) ;
          Flocal(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1]) ;
        }
    } //gp

    return 0;
}




int GrowthElem2DQua1::calcInternalForces()
{
  return 0;
}



void GrowthElem2DQua1::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'GrowthElem2DQua1::elementContourplot'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
      vals2project[0] += outval[ii];

    vals2project[0] /= nGP;

    return;
}


void GrowthElem2DQua1::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
{
   double outval[50];
   //cout << vartype << '\t' << varindex << '\t' << index << endl;

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
    //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", outval[0], outval[1], outval[2], outval[3]);
    stressrecovery_extrapolate_Quadrilateral(degree, outval, &vals2project[0]);
    //printVector(vals2project);

  return;
}


void GrowthElem2DQua1::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    double  detF, fact, dvol, dvol0;
    double  Jac, param[2];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), stre(9);
    MatrixXd  F(3,3), FgInv(3,3), Fe(3,3);
    F.setZero();
    double  gp1 = 1.0 + elmDat[10]*timeFunction[0].prop ;
    FgInv(0,0) = 1.0/gp1;    FgInv(0,1) = 0.0;        FgInv(0,2) = 0.0;
    FgInv(1,0) = 0.0;        FgInv(1,1) = 1.0/gp1;    FgInv(1,2) = 0.0;
    FgInv(2,0) = 0.0;        FgInv(2,1) = 0.0;        FgInv(2,2) = 1.0;

    int   err,  isw,  count,  count1, ll, ii, jj, gp;

    double dt = myTime.dt;

    double xNode[nlbfU], yNode[nlbfU], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsQuad(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        F(0,0) = computeValue(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValue(0, dN_dy);                  // eps_xy
        F(1,0) = computeValue(1, dN_dx);                  // eps_yx
        F(1,1) = computeValue(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

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

        //Fe = F*FgInv;
        Fe = F;

        pres = 0.0;
        MatlData->computeMechanicalStress(Fe, pres, stre);
        count++;
        ll += nivGP;

        if(varindex < 9)
           outval[gp] = stre[varindex];
        else if(varindex == 9)
           outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6* stre[1]*stre[1])/2);
        else if(varindex == 10)
           outval[gp] = (stre[0]+stre[4]+stre[8])/3.0;
    } //gp

    return;
}



void GrowthElem2DQua1::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void GrowthElem2DQua1::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



int GrowthElem2DQua1::calcLoadVector(VectorXd& Flocal)
{
  return 0;
}





int  GrowthElem2DQua1::calcError(int index)
{
    // compute error
    ///////////////////////////////////////////////////////////

  int   err,  isw,  count,  count1, ll, ii, jj, kk, gp;

  double  param[2], Jac, dvol0, dvol, xx, yy, rad, theta, val, fact;
  double  dispExact[2], dispNum[2], streDev[4], pbar;
  double  detF, streExac[4], streNum[4], cc[4][4];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), stre(9);

  double xNode[nlbfU], yNode[nlbfU], geom[3];
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

          //cout << dispExact[0] << '\t' << dispNum[0] << endl;
          //cout << dispExact[1] << '\t' << dispNum[1] << endl;

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

        F(0,0) = computeValue(0, dN_dx) + 1.0;            // eps_xx
        F(0,1) = computeValue(0, dN_dy);                  // eps_xy
        F(1,0) = computeValue(1, dN_dx);                  // eps_yx
        F(1,1) = computeValue(1, dN_dy) + 1.0;            // eps_yy

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        growthfunction_beam_singlelayer(timeFactor,  geom,  Fg);
        //growthfunction_beam_multilayer(timeFactor, geom,  Fg);
        //growthfunction_isotropic(1, ndim,  growthFactor,  Fg);
        FgInv = Fg.inverse();
        Jg = Fg.determinant();

        // elastic part of the deformation gradient
        Fe = F*FgInv;

        pres = 0.0;
        MatlData->computeMechanicalStress(Fe, pres, stre);

        pbar = 0.0;
        pres = (stre(0)+stre(4)+stre(8))/3.0;
        //cout << pbar << '\t' << pres << endl;

        val = pres - pbar;
        fact = val*val;

        elemError += ( fact * dvol );
    }//gp
  }

  return 0;
}



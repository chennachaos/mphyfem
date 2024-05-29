
#include "MagnetoMech2DMixedBase21.h"
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



MagnetoMech2DMixedBase21::MagnetoMech2DMixedBase21()
{
}



MagnetoMech2DMixedBase21::~MagnetoMech2DMixedBase21()
{
}


int MagnetoMech2DMixedBase21::calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err = 0, index, ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU), Np(nlbfP);

    double  detF, detFn, fact, fact1, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, Rp, Jhat, thetahat;
    double  param[2], bforce[2], force[2];
    double  veloCur[2], acceCur[2], sig[2];

    double  timeFactor = timeFunction[0].prop, timeFactorPrev = timeFactor-myTime.dt;

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

    MatrixXd  Cmat(9,9), Gc(3,9), F(3,3), Fn(3,3);
    VectorXd  stre(9);
    Gc.setZero();
    F.setZero();
    Fn.setZero();
    stre.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);

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

        acceCur[0] = computeValueDotDotCur(0, Nu);
        acceCur[1] = computeValueDotDotCur(1, Nu);

        if(detF < 0.0)
        {
          throw runtime_error("Negative Jacobian in the element");
        }

        adjust_deformationgradient_2D(sss, finite, F);


        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur  += dvol;

        if(finite)  dvol *= F(2,2);

        // evaluate pressure at the quadrature points

        if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
            BernsteinBasisFunsQuad(1, param[0], param[1], &Np(0));
        else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
            BernsteinBasisFunsTria(1, param[0], param[1], &Np(0));

        pres = computeValue2(0, Np);

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

          Flocal1(TI)   += (bb5*force[0] - sig[0]) ;
          Flocal1(TIp1) += (bb5*force[1] - sig[1]) ;

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

          Flocal2(ii) -= bb5*Rp;

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

    //printMatrix(Kuu);
    //printMatrix(Kpu);
    //printMatrix(Kup);
    //printMatrix(Kpp);
    //printVector(Flocal1);
    //printVector(Flocal2);

    return 0;
}




void MagnetoMech2DMixedBase21::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'MagnetoMech2DMixedBase21::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
    {
      vals2project[0] += outval[ii];
    }

    vals2project[0] /= nGP;

    return;
}


void MagnetoMech2DMixedBase21::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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


void MagnetoMech2DMixedBase21::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    int  err = 0, ii, jj, gp;

    double  detF, dvol, Jac, fact, dvol0, param[2], Np[nlbfP];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), stre(9);

    double  timeFact = elmDat[10]*timeFunction[0].prop;
    MatrixXd  F(3,3), Fn(3,3);
    F.setZero();
    Fn.setZero();
    double dt   = myTime.dt;

    double xNode[9], yNode[9], geom[3];
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


    for(int gp=0; gp<nGP; gp++)
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

        computeDefGrad2D(dN_dx, dN_dy, F);

        detF = F(0,0)*F(1,1) - F(0,1)*F(1,0);

        adjust_deformationgradient_2D(sss, finite, F);


        // basis functions for pressure dof
        if(ELEM_SHAPE == ELEM_SHAPE_QUAD_BERNSTEIN)
            BernsteinBasisFunsQuad(1, param[0], param[1], Np);
        else if(ELEM_SHAPE == ELEM_SHAPE_TRIA_BERNSTEIN)
            BernsteinBasisFunsTria(1, param[0], param[1], Np);

        // evaluate pressure at the quadrature points
        pres = 0.0;
        for(ii=0; ii<nlbfP; ii++)
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);
          //pres += (Np[ii]*SolnData->var2[forAssyVecPres[ii]]);


        MatlData->computeMechanicalStress(F, pres, stre);

        if(varindex < 9)
           outval[gp] = stre[varindex];
        else if(varindex == 9)
           outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6* stre[1]*stre[1])/2);
        else if(varindex == 10)
          outval[gp] = pres;

        //outval[gp] = detF;

    }//gp

    return;
}



void MagnetoMech2DMixedBase21::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void MagnetoMech2DMixedBase21::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;

    return;
}






int  MagnetoMech2DMixedBase21::calcError(int index)
{
    // compute error
    ///////////////////////////////////////////////////////////

  int   err,  isw,  count,  count1, ll, ii, jj, kk, gp;

  double  param[2], Jac, dvol0, dvol, rad, theta, val, fact;
  double  dispExact[2], dispNum[2], streDev[4], pbar;
  double  F[4], detF, F33, streExac[4], streNum[4], cc[4][4], Np[4];

  VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

  double xNode[9], yNode[9], geom[3];
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

          // evaluate pressure at the quadrature points
          BernsteinBasisFunsQuad(1, param[0], param[1], Np);
          pres = 0.0;
          for(ii=0; ii<4; ii++)
            pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);

          pbar = 0.0;

          val = pres - pbar;
          fact = val*val;

          elemError += ( fact * dvol );
    }//gp
  }

  return 0;
}



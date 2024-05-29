
/* Node ordering for the linear triangular element
 3
 | |
 |   |
 1-----2
*/


#include "BernsteinElem2DElecMechTri11.h"
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



BernsteinElem2DElecMechTri11::BernsteinElem2DElecMechTri11()
{
  ndim   = 2;
  degree = 1;
  npElem = 3;
  nlbfU  = 3;
  nlbfF  = 3;
  ndof   = 2;
  nsize  = nlbfU*ndof;
}

BernsteinElem2DElecMechTri11::~BernsteinElem2DElecMechTri11()
{
}


void BernsteinElem2DElecMechTri11::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double BernsteinElem2DElecMechTri11::computeVolume(bool configflag)
{
    double  dvol, Jac, param[2];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    int   ii, gp;

    double xNode[3], yNode[3], xx, yy;
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

        if(configflag)
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

    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

    return elemVol;
}



int BernsteinElem2DElecMechTri11::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != nsize) || (Mlocal.cols() != nsize) )
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = elmDat[5]*computeVolume(true)/3.0;

      for(int ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
  else
  {
    int  ii, jj, gp, TI, TIp1, TJ, TJp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  fact, dvol, Jac, bb1, cc1, param[2];
    double rho0 = elmDat[5] ;

    double xNode[3], yNode[3], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    int nGPt=3;
    getGaussPointsTriangle(nGPt, gausspoints1, gausspoints2, gaussweights);

    elemVol=0.0;
    for(gp=0; gp<nGPt; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol = gaussweights[gp]*(thick*Jac);

        xx = yy= 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        if(axsy)
          dvol *= 2.0*PI*yy;

        elemVol += dvol;

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
    }

    return 0;
}



double  BernsteinElem2DElecMechTri11::calcCriticalTimeStep(bool flag)
{
    int  ii;

    double xNode[3], yNode[3];
    for(ii=0;ii<npElem;ii++)
    {
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

    double  dtCric = charlen/wave_speed;

    return  dtCric;
}


int BernsteinElem2DElecMechTri11::calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp, TI, TIp1, TJ, TJp1;

    VectorXd  Nu(nlbfU), dNu_dx(nlbfU), dNu_dy(nlbfU);
    VectorXd  Nf(nlbfU), dNf_dx(nlbfU), dNf_dy(nlbfU);

    double  F[4], detF, trF, F33, fact, fact1, dvol, dvol0, Jac, volstrain, phi, gradphi[2];
    double  bb1, bb2, bb3, bb4, bb5, cc1, cc2, cc3, cc4, cc5, dUdJ, d2UdJ2;
    double  stre[4], Cmat[4][4], D11[4][4], bc[2][4], param[2], bforce[2], force[2];
    double  Bmat[2][4], Amat[2][2], tarr[4],  tarr2[4];
    double  veloCur[2], acceCur[2], sig[2];

    double BULK = matDat[0];
    double mu   = matDat[1];
    //double LAMBDA = BULK-r2d3*mu;
    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double af   = SolnData->td(2);
    double dt   = myTime.dt;
    double acceFact1 = SolnData->td(5);
    double acceFact2 = acceFact1;
    double FiniteFact = (finite == true) ? 1.0 : 0.0;

    double xNode[3], yNode[3], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    // resize local matrices and initialise them to zero
    if(Kuu.rows() < nsize || Kuu.cols() < nsize)   Kuu.resize(nsize, nsize);
    Kuu.setZero();

    if(Kuf.rows() < nsize || Kuf.cols() < nlbfF)     Kuf.resize(nsize, nlbfF);
    Kuf.setZero();
    if(Kfu.rows() < nlbfF || Kfu.cols() < nsize)     Kfu.resize(nlbfF, nsize);
    Kfu.setZero();
    if(Kff.rows() < nlbfF || Kff.cols() < nlbfF)     Kff.resize(nlbfF, nlbfF);
    Kff.setZero();

    if(FlocalU.rows() < nsize)   FlocalU.resize(nsize);
    FlocalU.setZero();
    if(FlocalF.rows() < nlbfF)   FlocalF.resize(nlbfF);
    FlocalF.setZero();

    elemVolOrig=0.0;
    elemVolCur=0.0;

    //Stokes2DEx1  analy;
    //Elasticity2DMixedEx1  analy(BULK, mu);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);


        Nf     = Nu;
        dNf_dx = dNu_dx;
        dNf_dy = dNu_dy;

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F[0] = computeValueCur(0, dNu_dx) + 1.0;
        F[2] = computeValueCur(0, dNu_dy);
        F[1] = computeValueCur(1, dNu_dx);
        F[3] = computeValueCur(1, dNu_dy) + 1.0;

        volstrain = F[0] + F[3] - 2.0;

        detF = F[0]*F[3] - F[1]*F[2];

        if(detF < 0.0)   return 1;

        acceCur[0] = computeValueDotDotCur(0, Nu);
        acceCur[1] = computeValueDotDotCur(1, Nu);

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

        ////matlib2d_(matDat, F, &F33, stre, Cmat[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;
        if(err !=0)          return 1;

        Cmat[0][0] = 12.6;  Cmat[0][1] =  5.3;  Cmat[0][2] = 0.0;  Cmat[0][3] = 0.0;
        Cmat[1][0] =  5.3;  Cmat[1][1] = 11.7;  Cmat[1][2] = 0.0;  Cmat[1][3] = 0.0;
        Cmat[2][0] =  0.0;  Cmat[2][1] =  0.0;  Cmat[2][2] = 0.0;  Cmat[2][3] = 0.0;
        Cmat[3][0] =  0.0;  Cmat[3][1] =  0.0;  Cmat[3][2] = 0.0;  Cmat[3][3] = 3.53;
        for(ii=0; ii<4; ii++)
        {
            for(jj=0; jj<4; jj++)
            {
                Cmat[ii][jj] *=  10.0e4;
            }
        }
        Bmat[0][0] =  0.0;  Bmat[0][1] =  0.0;  Bmat[0][2] = 0.0;  Bmat[0][3] = 17.0;
        Bmat[1][0] = -6.5;  Bmat[1][1] = 23.3;  Bmat[1][2] = 0.0;  Bmat[1][3] =  0.0;
        for(ii=0; ii<2; ii++)
        {
            for(jj=0; jj<4; jj++)
            {
                Bmat[ii][jj] *=  10.0e-6;
            }
        }

        Amat[0][0] =  1.51e-11; Amat[0][1] =  0.0;
        Amat[1][0] =  0.0;      Amat[1][1] =  1.3e-11;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &Nu(0), &dNu_dx(0), &dNu_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);

          dNf_dx = dNu_dx;
          dNf_dy = dNu_dy;
        }
        elemVolCur  += dvol;

        if(finite)
          dvol *= F33;

        // evaluate phi at the quadrature points
        phi = 0.0;
        gradphi[0] = gradphi[1] = 0.0;
        for(ii=0; ii<nlbfF; ii++)
        {
          phi += (Nf[ii]*SolnData->var3[nodeNums[ii]]);

          gradphi[0] += (Nf[ii]*SolnData->var3[nodeNums[ii]]);
          gradphi[1] += (Nf[ii]*SolnData->var3[nodeNums[ii]]);
        }


        //for(ii=0;ii<4;ii++)
          //cout << D11[ii][0] << '\t' << D11[ii][1] << '\t' << D11[ii][2] << '\t' << D11[ii][3] << endl;

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

        //force[0] = analy.computeForce(0, xx, yy, 0.0) ;
        //force[1] = analy.computeForce(1, xx, yy, 0.0) ;

        force[0] -= rho0*acceCur[0];
        force[1] -= rho0*acceCur[1];

        if(finite)
        {
          d2UdJ2 = 1.0;
        }
        else
        {
          d2UdJ2 = 1.0;
        }

        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dNu_dx[ii];
          bb2 = dvol*dNu_dy[ii];
          bb4 = dvol*Nu[ii];
          bb5 = dvol0*Nu[ii];

          bc[0][0] = bb1 * Cmat[0][0] + bb2 * Cmat[3][0];
          bc[0][1] = bb1 * Cmat[0][1] + bb2 * Cmat[3][1];
          bc[0][2] = bb1 * Cmat[0][3] + bb2 * Cmat[3][3];

          bc[1][0] = bb2 * Cmat[1][0] + bb1 * Cmat[3][0];
          bc[1][1] = bb2 * Cmat[1][1] + bb1 * Cmat[3][1];
          bc[1][2] = bb2 * Cmat[1][3] + bb1 * Cmat[3][3];

          TI   = 2*ii;
          TIp1 = TI+1;

          sig[0] = bb1*stre[0] + bb2*stre[3];
          sig[1] = bb1*stre[3] + bb2*stre[1];

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
            fact  += af*(sig[0]*cc1+sig[1]*cc2)*FiniteFact;

            Kuu(TI,   TJ)    += fact ;
            Kuu(TIp1, TJp1)  += fact ;

            Kuu(TI,   TJ)   +=  af*(bc[0][0] * cc1 + bc[0][2] * cc2) ;
            Kuu(TI,   TJp1) +=  af*(bc[0][1] * cc2 + bc[0][2] * cc1) ;
            Kuu(TIp1, TJ)   +=  af*(bc[1][0] * cc1 + bc[1][2] * cc2) ;
            Kuu(TIp1, TJp1) +=  af*(bc[1][1] * cc2 + bc[1][2] * cc1) ;
          }
        }

        for(ii=0; ii<nlbfF; ii++)
        {
          bb1 = dvol*dNf_dx[ii];
          bb2 = dvol*dNf_dy[ii];
          bb4 = dvol*Nf[ii];
          bb5 = dvol0*Nf[ii];

          FlocalF(ii) -= bb5*fact;

          tarr[0] = bb1 * Bmat[0][0] + bb2 * Bmat[1][0];
          tarr[1] = bb1 * Bmat[0][1] + bb2 * Bmat[1][1];
          tarr[2] = bb1 * Bmat[0][2] + bb2 * Bmat[1][2];
          tarr[3] = bb1 * Bmat[0][3] + bb2 * Bmat[1][3];


          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            tarr2[0] = tarr[0]*dNu_dx[jj]+tarr[3]*dNu_dy[jj];
            tarr2[1] = tarr[1]*dNu_dy[jj]+tarr[3]*dNu_dx[jj];

            Kfu(ii, TJ)    += tarr2[0];
            Kfu(ii, TJp1)  += tarr2[1];

            Kuf(TJ,   ii)  += tarr2[0];
            Kuf(TJp1, ii)  += tarr2[1];
          }


          tarr[0] = bb1 * Amat[0][0] + bb2 * Amat[1][0];
          tarr[1] = bb1 * Amat[0][1] + bb2 * Amat[1][1];

          for(jj=0; jj<nlbfF; jj++)
          {
            Kff(ii, jj)    += (dNf_dx[jj]*tarr[0] + dNf_dy[jj]*tarr[1]);
          }
        }
    }//gp

    printMatrix(Kuu);
    printMatrix(Kfu);
    printMatrix(Kff);

    return 0;
}



int BernsteinElem2DElecMechTri11::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ll, ii, jj, gp, TI, TIp1;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    double  detF, fact, dvol, dvol0, Jac;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[4], cc[4][4], param[2], bforce[2], force[2];
    //double  F[3][3] = {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
    double  F[4], F33;

    double mu   = matDat[1] ;
    double rho0 = elmDat[5] ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    double dt   = myTime.dt;
    double tPrev = myTime.cur ;

    //LinearElasticTransient2D  analy(rho0, matDat[0], matDat[1]);
    //LinearElasticTransientFIC2D  analy(rho0, mu);
    //LinearElastic2DEx1  analy(matDat[0], matDat[1]);

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    elemVolOrig=0.0;
    elemVolCur=0.0;
    for(gp=0;gp<nGP;gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        dvol0 = gaussweights[gp]*(Jac*thick);
        dvol  = dvol0;
        elemVolOrig += dvol;

        //cout << dN_dx[0] << '\t' << dN_dx[1] << '\t' << dN_dx[2] << '\t' << dN_dx[3] << '\t' << dN_dx[4] << '\t' << dN_dx[5] << endl;
        //cout << dN_dy[0] << '\t' << dN_dy[1] << '\t' << dN_dy[2] << '\t' << dN_dy[3] << '\t' << dN_dy[4] << '\t' << dN_dy[5] << endl;  cout << endl;

        xx = yy= 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += N[ii]*xNode[ii];
          yy += N[ii]*yNode[ii];
        }

        //if(axsy)
          //dvol *= 2.0*PI*yy;

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        detF = F[0]*F[3] - F[1]*F[2];

        //F[0] = computeValuePrev(0, dN_dx) + 1.0;
        //F[2] = computeValuePrev(0, dN_dy);
        //F[1] = computeValuePrev(1, dN_dx);
        //F[3] = computeValuePrev(1, dN_dy) + 1.0;

        if(detF < 0.0)   return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions2D(CONFIG_DEFORMED, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);
          dvol = gaussweights[gp]*(Jac*thick);
        }
        elemVolCur += dvol;

        // ADJUST F33 fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
          {
            //F[2][2] = 1.0/sqrt(detF);
            F33 = 1.0/sqrt(detF);
          }
        }
        else if(sss == 2)    // plane strain
        {
          //F[2][2] = 1.0;
          F33 = 1.0;
        }

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        //ComputeStressAndTangent(matId, matDat, F, stre, cc, intVar1, intVar2, dt, nivGP, finite, isw, gp, sss);
        count++;
        ll += nivGP;

        //printf("strains ... %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", F[0]-1.0, F[1], F[2], F[3]-1.0);
        //printf("stress  ... %5d \t %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", elenum, stre[0], stre[1], stre[3], stre[2]);

        if(err !=0)    return 1;

        if(finite)
          dvol *= F33;

        // body force terms
        force[0] = rho0*bforce[0];
        force[1] = rho0*bforce[1];

        //force[0] = analy.computeForce(0, xx, yy, 0.0, tPrev);
        //force[1] = analy.computeForce(1, xx, yy, 0.0, tPrev);

        // Calculate Residual
        //==============================================
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

    //printVector(Flocal);

    return 0;
}




int BernsteinElem2DElecMechTri11::calcLoadVector(VectorXd& Flocal)
{
  if( NeumannBCs.size() == 0)
    return 0;


  int side, dir, ii, jj, nn, TI, TIp1, TJ, TJp1, gp, nGP=3;
  int face_node_map[3][3]={{0,3,1},{1,4,2},{2,5,0}};
  double  specVal, xNode[3], yNode[3], xx, yy, param[2], trac[2], normal[2], tang[2], N[3];
  double  Jac, dvol;

  vector<double>  gausspoints1, gaussweights;
  getGaussPoints1D(nGP, gausspoints1, gaussweights);

  //printVector(Flocal);  printf("\n\n\n");

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

    for(ii=0; ii<3; ii++)
    {
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

        //cout << " trac = " << trac[0] << '\t' << trac[1] << endl;

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

  //printVector(Flocal);  printf("\n\n\n");

  return 0;
}




int BernsteinElem2DElecMechTri11::calcInternalForces()
{
  return 0;
}



void BernsteinElem2DElecMechTri11::elementContourplot(int vartype, int varindex, int index)
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

              cout  << " Invalid Variable Type to project in 'BernsteinElem2DElecMechTri11::projectToNodes'" << endl;
              break;
    }

    assert(vals2project.size() >= 1);

    vals2project[0] = 0.0;
    for(int ii=0; ii<nGP; ii++)
      vals2project[0] += outval[ii];

    vals2project[0] /= nGP;

    return;
}



void BernsteinElem2DElecMechTri11::projectToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
{
   double outval[50];

   //printVector(forAssyVec);

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

              cout  << " Invalid Variable Type to project in 'BernsteinElem2DElecMechTri11::projectToNodes'" << endl;
              break;
    }


    assert(vals2project.size() == npElem);

    // project from quadrature points to nodes/control points
    stressrecovery_extrapolate_Triangle(degree, outval, &vals2project[0]);

      //vector<double>  valtemp;
      //valtemp = vals2project;

      //vals2project[0] = valtemp[0];
      //vals2project[1] = valtemp[1];
      //vals2project[2] = valtemp[2];
      //vals2project[3] = 2.0*( valtemp[3] - 0.25*(valtemp[0]+valtemp[1]) );
      //vals2project[4] = 2.0*( valtemp[4] - 0.25*(valtemp[1]+valtemp[2]) );
      //vals2project[5] = 2.0*( valtemp[5] - 0.25*(valtemp[2]+valtemp[0]) );

      //vals2project[3] = 0.5*valtemp[3] + 0.25*(valtemp[0]+valtemp[1]);
      //vals2project[4] = 0.5*valtemp[4] + 0.25*(valtemp[1]+valtemp[2]);
      //vals2project[5] = 0.5*valtemp[5] + 0.25*(valtemp[2]+valtemp[0]);

    //cout << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;
    //cout << vals2project[0] << '\t' << vals2project[1] << '\t' << vals2project[2] << '\t' << vals2project[3] << '\t' << vals2project[4] << '\t' << vals2project[5] << endl; cout << endl;

    //printVector(vals2project);

    //printVector(forAssyVec);

    return;
}


void BernsteinElem2DElecMechTri11::projectStress(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    double  detF, dvol, Jac, param[2];
    double  stre[9], cc[4][4];
    double  F[4], F33;

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU);

    int   err,  isw,  count,  count1, ll = 0, ii, jj, gp;

    double dt   = myTime.dt;

    double xNode[6], yNode[6], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
    }

    vector<double>  gausspoints1, gausspoints2, gaussweights;
    getGaussPointsTriangle(nGP, gausspoints1, gausspoints2, gaussweights);

    count = 1;   ll = 0;   err = 0;   isw = 3;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];

        GeomData->computeBasisFunctions2D(CONFIG_ORIGINAL, ELEM_SHAPE_TRIA_BERNSTEIN, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), Jac);

        F[0] = computeValue(0, dN_dx) + 1.0;
        F[2] = computeValue(0, dN_dy);
        F[1] = computeValue(1, dN_dx);
        F[3] = computeValue(1, dN_dy) + 1.0;

        //detF = F[0][0]*F[1][1] - F[0][1]*F[1][0];
        detF = F[0]*F[3] - F[1]*F[2];

        // ADJUST F[2][2] fOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC

        if(sss == 1)  // plane stress
        {
          if(finite)
            //F[2][2] = 1.0/sqrt(detF);
            F33 = 1.0/sqrt(detF);
        }
        else if(sss == 2)    // plane strain
        {
          //F[2][2] = 1.0;
          F33 = 1.0;
        }

        //cout << " nivGP = " << nivGP << endl;
        //for(int i=0; i<nGP*nivGP; i++)  cout << '\t' << i << '\t' << intVar1[i] << '\t' << intVar2[i] << endl;

        //matlib2d_(matDat, F, &F33, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &sss, &isw, &err, &count, NULL);
        //ComputeStressAndTangent(matId, matDat, F, stre, cc, intVar1, intVar2, dt, nivGP, finite, isw, gp, sss);
        count++;
        ll += nivGP;

        //printf("stress  %14.12f \t %14.12f \t %14.12f \t %14.12f \t %14.12f \t %14.12f \n ", stre[0], stre[1], stre[2], stre[3], stre[4], stre[5]);

        if(varindex < 9)
          outval[gp] = stre[varindex];
        else if(varindex == 9)
          outval[gp] = sqrt((pow(stre[0]-stre[4],2.0) + pow(stre[4]-stre[8], 2.0) + pow(stre[8]-stre[0], 2.0) + 6.0*stre[1]*stre[1])/2.0);
        else if(varindex == 10)
          outval[gp] = (stre[0]+stre[4]+stre[8])/3.0;
          //outval[gp] = pres;
        else
          outval[gp] = 0.0;
    }//gp

    return;
}



void BernsteinElem2DElecMechTri11::projectStrain(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
  return;
}



void BernsteinElem2DElecMechTri11::projectInternalVariable(bool extrapolateFlag, int vartype, int varindex, int index, double* outval)
{
    assert( (ivar.var.rows() > 0) && (varindex < nivGP) );

    for(int gp=0; gp<nGP; gp++)
    {
       outval[gp] = ivar.var(varindex, gp);
    }//gp

    //cout << varindex << '\t' << outval[0] << '\t' << outval[1] << '\t' << outval[2] << endl;
    return;
}






/* Node ordering for the quadratic triangular element with bubble function
 3
 | |
 |   |
 6     5
 |   7   |
 |         |
 1-----4-----2
*/


#include "BernsteinElem2DINSTriaP2bP1dc.h"
#include "BasisFunctionsBernstein.h"
#include "headersBasic.h"
#include "KimMoinFlow.h"
#include "elementutilitiescfd.h"
#include "KimMoinFlow.h"
#include "ViscoElasticModels.h"
#include "QuadratureUtil.h"
#include "BasisFunctionsLagrange.h"


using namespace std;


BernsteinElem2DINSTriaP2bP1dc::BernsteinElem2DINSTriaP2bP1dc()
{
  ELEM_TYPE = ELEM_TRIA_BERNSTEIN_P2bP1dc;

  degree  = 2;
  npElem  = 7;

  nlbfU   = 7;
  ndof    = 2;
  nsize   = nlbfU*ndof;
  nlbfP   = 3;

  nGP = 7;
}


BernsteinElem2DINSTriaP2bP1dc::~BernsteinElem2DINSTriaP2bP1dc()
{

}




void BernsteinElem2DINSTriaP2bP1dc::findLocalCoordinates(vector<myPoint>& nodeCoords, myPoint& coords_global, myPoint& coords_local)
{
    //cout << "BernsteinElem2DINSTriaP2bP1dc::findLocalCoordinates" << endl;
    //printVector(nodeNums);

    double xNode[3], yNode[3];
    for(int ii=0; ii<3; ii++)
    {
      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    //cout << "global coordinates ... " << '\t' << coords_global[0] << '\t' <<  coords_global[1] << endl;

    pointInverseTria3node(xNode, yNode, &coords_global[0], &coords_local[0]);

    //cout << "local coordinates ... " << '\t' << coords_local[0] << '\t' <<  coords_local[1] << endl;

    return;
}





void BernsteinElem2DINSTriaP2bP1dc::prepareElemData(vector<myPoint>& nodeCoords)
{
    // compute Volume and basis function derivatives

    double  dvol, Jac, param[2];

    int   ii, gp;

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      //xNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][0];
      //yNode[ii] = GeomData->NodePosOrig[SolnData->node_map_new_to_old[nodeNums[ii]]][1];

      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    // Gauss point coordinates and weights
    vector<double>  gpts1(nGP), gpts2(nGP), gwts(nGP);
    getGaussPointsTriangle(nGP, gpts1, gpts2, gwts);

    Nv.resize(nGP);
    dNvdx.resize(nGP);
    dNvdy.resize(nGP);

    Np.resize(nGP);
    dNpdx.resize(nGP);
    dNpdy.resize(nGP);

    for(int ii=0; ii<nGP; ii++)
    {
      Nv[ii].resize(nlbfU);
      dNvdx[ii].resize(nlbfU);
      dNvdy[ii].resize(nlbfU);

      Np[ii].resize(nlbfU);
      dNpdx[ii].resize(nlbfU);
      dNpdy[ii].resize(nlbfU);

      Nv[ii].setZero();
      dNvdx[ii].setZero();
      dNvdy[ii].setZero();

      Np[ii].setZero();
      dNpdx[ii].setZero();
      dNpdy[ii].setZero();
    }

    //printVector(nodeNums);

    elemVolGP.resize(nGP);

    elemVol=0.0;
    for(gp=0; gp<nGP; gp++)
    {
      param[0] = gpts1[gp];
      param[1] = gpts2[gp];

      computeBasisFunctions2D(false, ELEM_TYPE, degree, param, xNode, yNode, &Nv[gp][0], &dNvdx[gp][0], &dNvdy[gp][0], Jac);

      //cout << " Jac = " << Jac << endl;
      if(Jac < 0.0)
      {
        cout << " Negative Jacobian in 'BernsteinElem2DINSTriaP2bP1dc::prepareElemData' at gp " << gp << endl;
        cout << " Jac = " << Jac << endl;
        exit(1);
      }

      dvol = gwts[gp]*Jac;
      elemVolGP[gp] = dvol;

      elemVol += dvol;

      Np[gp][0] = 1.0-param[0]-param[1];  Np[gp][1] = param[0];  Np[gp][2] = param[1];

      //computeBasisFunctions2D(false, ELEM_TYPE, 1, param, xNode, yNode, &Np[gp][0], &dNpdx[gp][0], &dNpdy[gp][0], Jac);
    }//gp

    //cout << " elemVol = " << elemVol << endl;

    //charlen = sqrt(4.0*elemVol/PI);
    //charlen = 0.75*charlen;

    //charlen = (xNode[1]-xNode[0])*(xNode[1]-xNode[0])+(yNode[1]-yNode[0])*(yNode[1]-yNode[0]);
    //charlen = min(charlen, (xNode[2]-xNode[1])*(xNode[2]-xNode[1])+(yNode[2]-yNode[1])*(yNode[2]-yNode[1]));
    //charlen = min(charlen, (xNode[2]-xNode[0])*(xNode[2]-xNode[0])+(yNode[2]-yNode[0])*(yNode[2]-yNode[0]));

    charlen = (xNode[3]-xNode[4])*(xNode[3]-xNode[4])+(yNode[3]-yNode[4])*(yNode[3]-yNode[4]);
    charlen = min(charlen, (xNode[4]-xNode[5])*(xNode[4]-xNode[5])+(yNode[4]-yNode[5])*(yNode[4]-yNode[5]));
    charlen = min(charlen, (xNode[5]-xNode[3])*(xNode[5]-xNode[3])+(yNode[5]-yNode[3])*(yNode[5]-yNode[3]));

    charlen = sqrt(charlen);

    //cout << " charlen = " << charlen << endl;

    //printMatrix(Klocal);  printf("\n\n\n");  printVector(Flocal);

    globalDOFnums.resize(npElem*2);

    for(ii=0; ii<npElem; ii++)
    {
      globalDOFnums[ii*2]   = nodeNums[ii]*2;
      globalDOFnums[ii*2+1] = nodeNums[ii]*2 + 1;
    }

    // pressure DOFs
    nodeNumsPres.resize(3);
    ii = elenum*3;
    nodeNumsPres[0] = ii;
    nodeNumsPres[1] = ii+1;
    nodeNumsPres[2] = ii+2;

    return;
}





double  BernsteinElem2DINSTriaP2bP1dc::ResidualIncNavStokesAlgo1(vector<myPoint>& node_coods, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& veloDotVec, VectorXd& veloDotVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  FlocalVelo, VectorXd&  FlocalPres, double timeCur)
{
    double  b1, b2, b3, b4;
    double  velo[2], veloPrev[2], resi[2], gradTvel[2], pres;
    double  force[2], veloDot[2], divVelo;
    double  xx, yy, fact, fact2;
    double  grad[2][2], stress[2][2];

    int  gp, ii, jj, TI, TIp1;

    // material parameters
    //double  rho = elemData[0];
    //double  mu  = elemData[1];

    double  rho = FluidMatlData->getDensity();
    double  mu  = FluidMatlData->getDynamicViscosity();

    double  bforce[2] = {elemData[2],  elemData[3]};
    double  beta0  = elemData[5];
    double  betaSq = beta0*beta0;
    double  v_conv=-1.0, dtCric=1.0e10;
    double  beta=-1.0;
    double  v_diff = mu/rho/charlen;
    double  dtCric_visco = charlen*charlen*rho/mu/2.0;

    //KimMoinFlow  analy(rho, mu, 0.0);

    //loop over Gauss points and compute element residual
    FlocalVelo.setZero();
    FlocalPres.setZero();
    force[0] = 0.0;    force[1] = 0.0;

    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of velocity
          xx = 0.0; yy = 0.0;
          velo[0] = velo[1] = 0.0;
          grad[0][0] = grad[0][1] = 0.0;
          grad[1][0] = grad[1][1] = 0.0;
          for(ii=0; ii<6; ii++)
          {
            //xx += xNode[ii]*N[ii];
            //yy += yNode[ii]*N[ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];

            velo[0]    +=  b1*Nv[gp][ii];
            velo[1]    +=  b2*Nv[gp][ii];

            //b1 = 1.5*veloVec[TI]-0.5*veloVecPrev[TI];
            //b2 = 1.5*veloVec[TIp1]-0.5*veloVecPrev[TIp1];

            grad[0][0] += b1*dNvdx[gp][ii];
            grad[0][1] += b1*dNvdy[gp][ii];
            grad[1][0] += b2*dNvdx[gp][ii];
            grad[1][1] += b2*dNvdy[gp][ii];

            //veloPrev[0] += veloVecPrev[TI]*Nv[gp][ii];
            //veloPrev[1] += veloVecPrev[TIp1]*Nv[gp][ii];

            //b1 = veloDotVec[TI];
            //b2 = veloDotVec[TIp1];

            //veloDot[0]    +=  b1*Nv[gp][ii];
            //veloDot[1]    +=  b2*Nv[gp][ii];
          }

          pres = 0.0;
          for(ii=0; ii<3; ii++)
          {
            //pres += (1.5*presVec[nodeNums[ii]]-0.5*presVecPrev[nodeNums[ii]])*Np[gp][ii];
            b1 = presVec[nodeNums[ii]];

            pres  += b1*Np[gp][ii];
          }

          //fact = 2.0*mu*(grad(0,0)+grad(1,1))/3.0;
          //stress(0,0) = mu*(grad(0,0)+grad(0,0)) - fact - pres;
          //stress(0,1) = mu*(grad(0,1)+grad(1,0));
          //stress(1,0) = stress(0,1);
          //stress(1,1) = mu*(grad(1,1)+grad(1,1)) - fact - pres;

          stress[0][0] = mu*grad[0][0] - pres;
          stress[0][1] = mu*grad[0][1];
          stress[1][0] = mu*grad[1][0];
          stress[1][1] = mu*grad[1][1] - pres;

          //force[0] = 0.0;
          //force[1] = 0.0;
          //force[0] = analy.computeForce(0, xx, yy, timeCur);
          //force[1] = analy.computeForce(1, xx, yy, timeCur);

          divVelo = grad[0][0]+grad[1][1];

          gradTvel[0] = velo[0]*grad[0][0] + velo[1]*grad[0][1] ;
          gradTvel[1] = velo[0]*grad[1][0] + velo[1]*grad[1][1] ;

          //gradTvel[0] = velo[0]*grad(0,0) + velo[1]*grad(0,1) + velo[0]*divVelo;
          //gradTvel[1] = velo[0]*grad(1,0) + velo[1]*grad(1,1) + velo[1]*divVelo;

          resi[0] = rho*(force[0] - gradTvel[0]) ;
          resi[1] = rho*(force[1] - gradTvel[1]) ;

          for(ii=0; ii<6; ii++)
          {
            TI   = ii*2;
            TIp1 = TI+1;

            b1 = dNvdx[gp][ii]*elemVolGP[gp];
            b2 = dNvdy[gp][ii]*elemVolGP[gp];
            b4 = Nv[gp][ii]*elemVolGP[gp];

            FlocalVelo[TI]   += (b4*resi[0] - b1*stress[0][0] - b2*stress[0][1]);
            FlocalVelo[TIp1] += (b4*resi[1] - b1*stress[1][0] - b2*stress[1][1]);
          }

          v_conv = max(v_conv, velo[0]*velo[0]+velo[1]*velo[1]);

          fact = -elemVolGP[gp]*betaSq*divVelo;

          for(ii=0; ii<3; ii++)
          {
            FlocalPres[ii] += Np[gp][ii]*fact;
          }
    }

    dtCric = charlen/(sqrt(v_conv)+sqrt(v_conv+betaSq));

    return  min(dtCric, dtCric_visco);
}



double  BernsteinElem2DINSTriaP2bP1dc::ResidualIncNavStokesSemiImpl(vector<myPoint>& node_coods, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& veloDotVec, VectorXd& veloDotVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  FlocalVelo, VectorXd&  FlocalPres, double timeCur)
{
    double  b1, b2, b3, b4;
    double  velo[2], veloPrev[2], resi[2], gradTvel[2], pres;
    double  force[2], veloDot[2], divVelo;
    double  xx, yy, fact, fact2;
    double  grad[2][2], stress[2][2];

    int  gp, ii, jj, TI, TIp1;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    double  bforce[2] = {elemData[2],  elemData[3]};
    double  v_conv=-1.0, dtCric=1.0e10;
    double  dtCric_visco = charlen*charlen*rho/mu/4.0;

    //KimMoinFlow  analy(rho, mu, 0.0);

    //loop over Gauss points and compute element residual
    FlocalVelo.setZero();
    FlocalPres.setZero();
    force[0] = 0.0;    force[1] = 0.0;

    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of velocity
          xx = 0.0; yy = 0.0;
          velo[0] = velo[1] = 0.0;
          grad[0][0] = grad[0][1] = 0.0;
          grad[1][0] = grad[1][1] = 0.0;
          for(ii=0; ii<6; ii++)
          {
            //xx += xNode[ii]*N[ii];
            //yy += yNode[ii]*N[ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];

            velo[0]    +=  b1*Nv[gp][ii];
            velo[1]    +=  b2*Nv[gp][ii];

            //b1 = 1.5*veloVec[TI]-0.5*veloVecPrev[TI];
            //b2 = 1.5*veloVec[TIp1]-0.5*veloVecPrev[TIp1];

            grad[0][0] += b1*dNvdx[gp][ii];
            grad[0][1] += b1*dNvdy[gp][ii];
            grad[1][0] += b2*dNvdx[gp][ii];
            grad[1][1] += b2*dNvdy[gp][ii];

            //veloPrev[0] += veloVecPrev[TI]*Nv[gp][ii];
            //veloPrev[1] += veloVecPrev[TIp1]*Nv[gp][ii];

            //b1 = veloDotVec[TI];
            //b2 = veloDotVec[TIp1];

            //veloDot[0]    +=  b1*Nv[gp][ii];
            //veloDot[1]    +=  b2*Nv[gp][ii];
          }

          pres = 0.0;
          for(ii=0; ii<3; ii++)
          {
            //pres += (1.5*presVec[nodeNums[ii]]-0.5*presVecPrev[nodeNums[ii]])*Np[gp][ii];
            b1 = presVec[nodeNums[ii]];

            pres  += b1*Np[gp][ii];
          }

          //stress[0][0] = mu*(grad[0][0]+grad[0][0]) - pres;
          //stress[0][1] = mu*(grad[0][1]+grad[1][0]);
          //stress[1][0] = stress[0][1];
          //stress[1][1] = mu*(grad[1][1]+grad[1][1]) - pres;

          mu  = compute_viscosity_CarreauYasuda(grad);

          stress[0][0] = mu*grad[0][0] - pres;
          stress[0][1] = mu*grad[0][1];
          stress[1][0] = mu*grad[1][0];
          stress[1][1] = mu*grad[1][1] - pres;

          //force[0] = 0.0;
          //force[1] = 0.0;
          //force[0] = analy.computeForce(0, xx, yy, timeCur);
          //force[1] = analy.computeForce(1, xx, yy, timeCur);

          divVelo = grad[0][0]+grad[1][1];

          gradTvel[0] = velo[0]*grad[0][0] + velo[1]*grad[0][1] ;
          gradTvel[1] = velo[0]*grad[1][0] + velo[1]*grad[1][1] ;

          //gradTvel[0] = velo[0]*grad(0,0) + velo[1]*grad(0,1) + velo[0]*divVelo;
          //gradTvel[1] = velo[0]*grad(1,0) + velo[1]*grad(1,1) + velo[1]*divVelo;

          resi[0] = rho*(force[0] - gradTvel[0]) ;
          resi[1] = rho*(force[1] - gradTvel[1]) ;

          for(ii=0; ii<6; ii++)
          {
            TI   = ii*2;
            TIp1 = TI+1;

            b1 = dNvdx[gp][ii]*elemVolGP[gp];
            b2 = dNvdy[gp][ii]*elemVolGP[gp];
            b4 = Nv[gp][ii]*elemVolGP[gp];

            FlocalVelo[TI]   += (b4*resi[0] - b1*stress[0][0] - b2*stress[0][1]);
            FlocalVelo[TIp1] += (b4*resi[1] - b1*stress[1][0] - b2*stress[1][1]);
          }

          v_conv = max(v_conv, velo[0]*velo[0]+velo[1]*velo[1]);

          fact = -elemVolGP[gp]*divVelo;

          for(ii=0; ii<3; ii++)
          {
            FlocalPres[ii] += Np[gp][ii]*fact;
          }
    }

    dtCric = charlen/sqrt(v_conv);

    return  min(dtCric, dtCric_visco);
}




int  BernsteinElem2DINSTriaP2bP1dc::ResidualIncNavStokesAlgo2(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloDotPrev, VectorXd& presPrev, VectorXd&  FlocalPres)
{
      double  velo[2], grad[2][2];
      double  xx, yy, fact, fact2, b1, b2, v_conv;

      int  gp, ii, jj, TI, TIp1;

      // material parameters
      double  rho = elemData[0];
      double  mu  = elemData[1];
      double  v_diff = mu/rho/charlen;

      double  beta2 = elemData[5]*elemData[5];
      double  beta=-1.0;
      double  beta0=elemData[5];

      // time integration parameters
      double  af   = timeData[1];
      double  timefact = timeData[2];

      //loop over Gauss points and compute element residual
      FlocalPres.setZero();

      for(gp=0; gp<nGP; gp++)
      {
          // compute the gradient of displacement first
          xx = 0.0; yy = 0.0;
          velo[0] = velo[1] = 0.0;
          grad[0][0] = grad[0][1] = grad[1][0] = grad[1][1] = 0.0;

          for(ii=0; ii<6; ii++)
          {
            //xx += xNode[ii]*N[ii];
            //yy += yNode[ii]*N[ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloPrev[TI];
            b2 = veloPrev[TIp1];

            velo[0]    +=  b1*Nv[gp][ii];
            velo[1]    +=  b2*Nv[gp][ii];

            grad[0][0] += b1*dNvdx[gp][ii];
            grad[1][1] += b2*dNvdy[gp][ii];
          }

          v_conv = sqrt(velo[0]*velo[0]+velo[1]*velo[1]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta2 = beta*beta;

          fact = -elemVolGP[gp]*beta2*(grad[0][0]+grad[1][1]);

          for(ii=0; ii<3; ii++)
          {
            FlocalPres[ii] += Np[gp][ii]*fact;
          }
      }

  return 0;
}




double  BernsteinElem2DINSTriaP2bP1dc::calcCriticalTimeStep(double* elemData, double* timeData, VectorXd&  veloVec)
{
    int ii, jj, gp, TI, TIp1;

    double  velo[2], b1, b2, v_conv, dtCric=1.0e10;
    double  beta=-1.0;
    double  beta0=elemData[5];
    double  rho=elemData[0];
    double  mu =elemData[1];
    double  v_diff = mu/rho/charlen;
    double  dtCric_visco = charlen*charlen*rho/mu/2.0;

    for(gp=0; gp<nGP; gp++)
    {
          velo[0] = velo[1] = 0.0;
          for(ii=0; ii<6; ii++)
          {
            //xx += xNode[ii]*N[ii];
            //yy += yNode[ii]*N[ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];

            velo[0]  +=  b1*Nv[gp][ii];
            velo[1]  +=  b2*Nv[gp][ii];
          }

          v_conv = sqrt(velo[0]*velo[0]+velo[1]*velo[1]);

          //beta = max(beta0, max(v_conv, v_diff));
          //beta = max(beta0, v_conv);
          beta = beta0;

          dtCric = min(dtCric, charlen/(v_conv+beta));
    }

    return  min(dtCric, dtCric_visco);
}



/*
int  BernsteinElem2DINSTriaP2bP1dc::StiffnessAndResidualFullyImplicit(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloPrev2, VectorXd& veloCur, VectorXd& veloDotCur, VectorXd& presCur, MatrixXd& Kuu, MatrixXd& Kup, VectorXd& Fu, VectorXd& Fp, double timeCur)
{
    // Fully-implicit formulation

    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);
    //KimMoinFlow  analy(rho, mu, 1.0);

    int ii, jj, gp, TI, TIp1, count, TJ, TJp1;

    double  b1, b2, b3, b4, b5, b6, b7, b8;
    double  dvol, pres, fact, fact2, veloDiv, Da, Db;
    double  param[2], vel[2], velPrev[2], velExtp[2], velDot[2], force[2];

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0; ii<npElem; ii++)
    {
      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    VectorXd  gradTvel(2);
    MatrixXd  grad(2,2), gradN(2,2), stress(2,2);


    double rho = elemData[0];
    double mu  = elemData[1];

    double  am = timeData[1];
    double  af = timeData[2];
    double  acceFact = timeData[8];
    double  muTaf = mu*af;
    double  dt = timeData[5]/af;

    Kuu.setZero();
    Kup.setZero();
    Fu.setZero();
    Fp.setZero();

    //KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);

    double  tCur  = timeCur;
    double  tPrev = tCur - dt;

    //cout << " AAAAAAAAAA " << rho << '\t' << mu << endl;
    //cout << nGP << endl;

    for(gp=0; gp<nGP; gp++)
    {
        // compute the gradient of velocity
        xx = 0.0; yy = 0.0;
        vel[0] = vel[1] = 0.0;
        velPrev[0] = velPrev[1] = 0.0;
        velExtp[0] = velExtp[1] = 0.0;
        velDot[0] = velDot[1] = 0.0;
        grad.setZero();

        for(ii=0; ii<npElem; ii++)
        {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloCur[TI];
            b2 = veloCur[TIp1];

            vel[0]     += b1*Nv[gp][ii];
            vel[1]     += b2*Nv[gp][ii];

            grad(0,0)  += b1*dNvdx[gp][ii];
            grad(0,1)  += b1*dNvdy[gp][ii];
            grad(1,0)  += b2*dNvdx[gp][ii];
            grad(1,1)  += b2*dNvdy[gp][ii];

            velDot[0] += veloDotCur[TI]*Nv[gp][ii];
            velDot[1] += veloDotCur[TIp1]*Nv[gp][ii];
        }

        pres = 0.0;
        for(ii=0; ii<3; ii++)
        {
            pres += presCur[nodeNums[ii]]*Np[gp][ii];
        }

        //printMatrix(grad);

        // this is pseudo-stress
        stress = mu*grad;
        stress(0,0) -= pres;
        stress(1,1) -= pres;

        //cout << vel[0] << '\t' << vel[1] << endl;
        //cout << stress(0,0) << '\t' << stress(0,1) << '\t' << stress(1,0) << '\t' << stress(1,1) << endl;

        force[0] = 0.0;
        force[1] = 0.0;
        //force[0] = analy.computeForce(0, xx, yy, tCur);
        //force[1] = analy.computeForce(1, xx, yy, tCur);

        //force[0] = af*analy.computeForce(0, xx, yy, 0.0, tCur)+(1.0-af)*analy.computeForce(0, xx, yy, 0.0, tPrev);
        //force[1] = af*analy.computeForce(1, xx, yy, 0.0, tCur)+(1.0-af)*analy.computeForce(1, xx, yy, 0.0, tPrev);

        gradTvel[0] = grad(0,0)*vel[0] + grad(0,1)*vel[1];
        gradTvel[1] = grad(1,0)*vel[0] + grad(1,1)*vel[1];


        force[0] -= rho*( velDot[0] + gradTvel[0] ) ;
        force[1] -= rho*( velDot[1] + gradTvel[1] ) ;
        veloDiv = grad.trace();

        dvol = elemVolGP[gp];

        //cout << force[0] << '\t' << force[1] << endl;

        for(ii=0; ii<npElem; ii++)
        {
            TI   = 2*ii;
            TIp1 = TI+1;

            b1 = dNvdx[gp][ii]*dvol;
            b2 = dNvdy[gp][ii]*dvol;
            b4 = Nv[gp][ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;

            b8 = af*b4;

            //cout << ii << '\t' << TI << '\t' << TIp1 << endl;

            Fu(TI)   += (b4*force[0] - b1*stress(0,0) - b2*stress(0,1) );
            Fu(TIp1) += (b4*force[1] - b1*stress(1,0) - b2*stress(1,1) );

            for(jj=0; jj<npElem; jj++)
            {
               TJ   = 2*jj;
               TJp1 = TJ+1;

               //cout << jj << '\t' << TJ << '\t' << TJp1 << endl;

               fact2 = rho*acceFact*Nv[gp][jj];

               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += (b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj]);

               Kuu(TI,   TJ)   += fact;
               Kuu(TIp1, TJp1) += fact;

               // convection term

               gradN = grad*(rho*Nv[gp][jj]);

               Db = rho*(vel[0]*dNvdx[gp][jj] + vel[1]*dNvdy[gp][jj]);

               gradN(0,0) += Db;
               gradN(1,1) += Db;

               Kuu(TI,   TJ)   += (b8*gradN(0,0));
               Kuu(TI,   TJp1) += (b8*gradN(0,1));
               Kuu(TIp1, TJ)   += (b8*gradN(1,0));
               Kuu(TIp1, TJp1) += (b8*gradN(1,1));
            }

            for(jj=0; jj<3; jj++)
            {
               // pressure term
               Kup(TI,   jj) -= (b1*af*Np[gp][jj]);
               Kup(TIp1, jj) -= (b2*af*Np[gp][jj]);
            }
        }

        //cout << velPrev[0] << '\t' << velPrev[1] << endl;

        for(ii=0; ii<3; ii++)
        {
            b8 = Np[gp][ii]*dvol;

            Fp(ii) += (b8*veloDiv);
        }

    }//gp

    //printMatrix(Kuu);
    //printVector(Fu);

    return 0;
}
*/



/*
int  BernsteinElem2DINSTriaP2bP1dc::StiffnessAndResidualFullyImplicit(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloPrev2, VectorXd& veloCur, VectorXd& veloDotCur, VectorXd& presCur, MatrixXd& Kuu, MatrixXd& Kup, VectorXd& Fu, VectorXd& Fp, double timeCur)
{
    // Semi-implicit formulation - type A

    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);
    //KimMoinFlow  analy(rho, mu, 1.0);

    int ii, jj, gp, TI, TIp1, count, TJ, TJp1;

    double  b1, b2, b3, b4, b5, b6, b7, b8;
    double  dvol, pres, fact, fact2, veloDiv, Da, Db;
    double  param[2], force[2];

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0; ii<npElem; ii++)
    {
      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    VectorXd  vel(2), velPrev(2), velPrev2(2), velExtp(2), velDot(2), gradTvel(2);
    MatrixXd  grad(2,2), gradPrev(2,2), gradN(2,2), stress(2,2);


    double rho = elemData[0];
    double mu  = elemData[1];

    double  am = timeData[1];
    double  af = timeData[2];
    double  acceFact = timeData[8];
    double  muTaf = mu*af;
    double  dt = timeData[5]/af;

    Kuu.setZero();
    Kup.setZero();
    Fu.setZero();
    Fp.setZero();

    //KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);

    double  tCur  = timeCur;
    double  tPrev = tCur - dt;

    //cout << " AAAAAAAAAA " << rho << '\t' << mu << '\t' << tPrev << '\t' << tCur << endl;
    //cout << nGP << endl;

    for(gp=0; gp<nGP; gp++)
    {
        // compute the gradient of velocity
        xx = 0.0; yy = 0.0;
        vel[0]     = vel[1]     = 0.0;
        velPrev[0] = velPrev[1] = 0.0;
        velPrev2[0] = velPrev2[1] = 0.0;
        velExtp[0] = velExtp[1] = 0.0;
        velDot[0]  = velDot[1]  = 0.0;
        grad.setZero();
        gradPrev.setZero();

        for(ii=0; ii<npElem; ii++)
        {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloCur[TI];
            b2 = veloCur[TIp1];

            vel[0]     += b1*Nv[gp][ii];
            vel[1]     += b2*Nv[gp][ii];

            grad(0,0)  += b1*dNvdx[gp][ii];
            grad(0,1)  += b1*dNvdy[gp][ii];
            grad(1,0)  += b2*dNvdx[gp][ii];
            grad(1,1)  += b2*dNvdy[gp][ii];

            b1 = veloPrev[TI];
            b2 = veloPrev[TIp1];

            velPrev[0]  += b1*Nv[gp][ii];
            velPrev[1]  += b2*Nv[gp][ii];

            velDot[0]   += veloDotCur[TI]*Nv[gp][ii];
            velDot[1]   += veloDotCur[TIp1]*Nv[gp][ii];

            velPrev2[0] += veloPrev2[TI]*Nv[gp][ii];
            velPrev2[1] += veloPrev2[TIp1]*Nv[gp][ii];
        }

        pres = 0.0;
        for(ii=0; ii<3; ii++)
        {
            pres += presCur[nodeNums[ii]]*Np[gp][ii];
        }

        //printMatrix(grad);

        // this is pseudo-stress
        stress = mu*grad;
        stress(0,0) -= pres;
        stress(1,1) -= pres;

        //cout << vel[0] << '\t' << vel[1] << endl;
        //cout << stress(0,0) << '\t' << stress(0,1) << '\t' << stress(1,0) << '\t' << stress(1,1) << endl;

        force[0] = 0.0;
        force[1] = 0.0;
        //force[0] = analy.computeForce(0, xx, yy, tCur);
        //force[1] = analy.computeForce(1, xx, yy, tCur);

        //force[0] = af*analy.computeForce(0, xx, yy, 0.0, tCur)+(1.0-af)*analy.computeForce(0, xx, yy, 0.0, tPrev);
        //force[1] = af*analy.computeForce(1, xx, yy, 0.0, tCur)+(1.0-af)*analy.computeForce(1, xx, yy, 0.0, tPrev);
        //cout << force[0] << '\t' << force[1] << endl;

        velExtp = 2.0*velPrev - velPrev2;
        //velExtp = velPrev;

        velExtp = af*velExtp + (1.0-af)*velPrev;

        gradTvel = grad*velExtp;

        force[0] -= rho*( velDot[0] + gradTvel[0] ) ;
        force[1] -= rho*( velDot[1] + gradTvel[1] ) ;
        veloDiv = grad.trace();

        dvol = elemVolGP[gp];

        //cout << force[0] << '\t' << force[1] << endl;

        for(ii=0; ii<npElem; ii++)
        {
            TI   = 2*ii;
            TIp1 = TI+1;

            b1 = dNvdx[gp][ii]*dvol;
            b2 = dNvdy[gp][ii]*dvol;
            b4 = Nv[gp][ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;

            b8 = af*b4;

            //cout << ii << '\t' << TI << '\t' << TIp1 << endl;

            Fu(TI)   += (b4*force[0] - b1*stress(0,0) - b2*stress(0,1) );
            Fu(TIp1) += (b4*force[1] - b1*stress(1,0) - b2*stress(1,1) );

            for(jj=0; jj<npElem; jj++)
            {
               TJ   = 2*jj;
               TJp1 = TJ+1;

               //cout << jj << '\t' << TJ << '\t' << TJp1 << endl;

               fact2 = rho*acceFact*Nv[gp][jj];

               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += (b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj]);

               Kuu(TI,   TJ)   += fact;
               Kuu(TIp1, TJp1) += fact;

               // convection term

               Db = rho*(velExtp[0]*dNvdx[gp][jj] + velExtp[1]*dNvdy[gp][jj]);

               Kuu(TI,   TJ)   += (b8*Db);
               Kuu(TIp1, TJp1) += (b8*Db);
            }

            for(jj=0; jj<3; jj++)
            {
               // pressure term
               Kup(TI,   jj) -= (b1*af*Np[gp][jj]);
               Kup(TIp1, jj) -= (b2*af*Np[gp][jj]);
            }
        }

        //cout << velPrev[0] << '\t' << velPrev[1] << endl;

        for(ii=0; ii<3; ii++)
        {
            b8 = Np[gp][ii]*dvol;

            Fp(ii) += (b8*veloDiv);
        }

    }//gp

    //printMatrix(Kuu);
    //printVector(Fu);

    return 0;
}
*/


//
int  BernsteinElem2DINSTriaP2bP1dc::StiffnessAndResidualFullyImplicit(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloPrev2, VectorXd& veloCur, VectorXd& veloDotCur, VectorXd& presCur, MatrixXd& Kuu, MatrixXd& Kup, VectorXd& Fu, VectorXd& Fp, double dt, double timeCur)
{
    // Semi-implicit formulation - type B

    //Stokes2DEx1 analy;
    //Kovasznay  analy;
    //analy.SetPressure(0.0);
    //KimMoinFlow  analy(rho, mu, 1.0);

    int ii, jj, gp, TI, TIp1, count, TJ, TJp1;

    double  b1, b2, b3, b4, b5, b6, b7, b8;
    double  dvol, pres, fact, fact2, veloDiv, Da, Db;
    double  param[2], force[2];

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0; ii<npElem; ii++)
    {
      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    VectorXd  vel(2), velPrev(2), velExtp(2), velDot(2), gradTvel(2);
    MatrixXd  grad(2,2), gradPrev(2,2), gradN(2,2), stress(2,2);


    //double rho = elemData[0];
    //double mu  = elemData[1];

    double  rho = FluidMatlData->getDensity();
    double  mu  = FluidMatlData->getDynamicViscosity();

    double  am = timeData[1];
    double  af = timeData[2];
    double  acceFact = timeData[8];
    double  muTaf = mu*af;

    Kuu.setZero();
    Kup.setZero();
    Fu.setZero();
    Fp.setZero();

    //KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);

    double  tCur  = timeCur;
    double  tPrev = tCur - dt;

    //cout << " AAAAAAAAAA " << rho << '\t' << mu << '\t' << am << '\t' << af << '\t' << acceFact << endl;
    //cout << nGP << endl;

    for(gp=0; gp<nGP; gp++)
    {
        // compute the gradient of velocity
        xx = 0.0; yy = 0.0;
        vel[0]     = vel[1]     = 0.0;
        velPrev[0] = velPrev[1] = 0.0;
        velExtp[0] = velExtp[1] = 0.0;
        velDot[0]  = velDot[1]  = 0.0;
        grad.setZero();
        gradPrev.setZero();

        for(ii=0; ii<npElem; ii++)
        {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloCur[TI];
            b2 = veloCur[TIp1];

            vel[0]     += b1*Nv[gp][ii];
            vel[1]     += b2*Nv[gp][ii];

            grad(0,0)  += b1*dNvdx[gp][ii];
            grad(0,1)  += b1*dNvdy[gp][ii];
            grad(1,0)  += b2*dNvdx[gp][ii];
            grad(1,1)  += b2*dNvdy[gp][ii];

            b1 = veloPrev[TI];
            b2 = veloPrev[TIp1];

            velPrev[0]     += b1*Nv[gp][ii];
            velPrev[1]     += b2*Nv[gp][ii];

            gradPrev(0,0)  += b1*dNvdx[gp][ii];
            gradPrev(0,1)  += b1*dNvdy[gp][ii];
            gradPrev(1,0)  += b2*dNvdx[gp][ii];
            gradPrev(1,1)  += b2*dNvdy[gp][ii];

            velDot[0] += veloDotCur[TI]*Nv[gp][ii];
            velDot[1] += veloDotCur[TIp1]*Nv[gp][ii];
        }

        pres = 0.0;
        for(ii=0; ii<3; ii++)
        {
            pres += presCur[nodeNumsPres[ii]]*Np[gp][ii];
        }

        //printMatrix(grad);

        // this is pseudo-stress
        stress = mu*grad;
        stress(0,0) -= pres;
        stress(1,1) -= pres;

        //cout << vel[0] << '\t' << vel[1] << endl;
        //cout << stress(0,0) << '\t' << stress(0,1) << '\t' << stress(1,0) << '\t' << stress(1,1) << endl;

        force[0] = 0.0;
        force[1] = 0.0;
        //force[0] = analy.computeForce(0, xx, yy, tCur);
        //force[1] = analy.computeForce(1, xx, yy, tCur);

        //force[0] = af*analy.computeForce(0, xx, yy, 0.0, tCur)+(1.0-af)*analy.computeForce(0, xx, yy, 0.0, tPrev);
        //force[1] = af*analy.computeForce(1, xx, yy, 0.0, tCur)+(1.0-af)*analy.computeForce(1, xx, yy, 0.0, tPrev);
        //cout << force[0] << '\t' << force[1] << endl;

        gradTvel = gradPrev*(vel-velPrev) + grad*velPrev;

        force[0] -= rho*( velDot[0] + gradTvel[0] ) ;
        force[1] -= rho*( velDot[1] + gradTvel[1] ) ;
        veloDiv = grad.trace();

        dvol = elemVolGP[gp];

        //cout << force[0] << '\t' << force[1] << endl;

        for(ii=0; ii<npElem; ii++)
        {
            TI   = 2*ii;
            TIp1 = TI+1;

            b1 = dNvdx[gp][ii]*dvol;
            b2 = dNvdy[gp][ii]*dvol;
            b4 = Nv[gp][ii]*dvol;

            b5 = muTaf*b1;
            b6 = muTaf*b2;

            b8 = af*b4;

            //cout << ii << '\t' << TI << '\t' << TIp1 << endl;

            Fu(TI)   += (b4*force[0] - b1*stress(0,0) - b2*stress(0,1) );
            Fu(TIp1) += (b4*force[1] - b1*stress(1,0) - b2*stress(1,1) );

            for(jj=0; jj<npElem; jj++)
            {
               TJ   = 2*jj;
               TJp1 = TJ+1;

               //cout << jj << '\t' << TJ << '\t' << TJp1 << endl;

               fact2 = rho*acceFact*Nv[gp][jj];

               // time acceleration term
               fact = b4*fact2 ;

               // diffusion term
               fact += (b5*dNvdx[gp][jj]+b6*dNvdy[gp][jj]);

               Kuu(TI,   TJ)   += fact;
               Kuu(TIp1, TJp1) += fact;

               // convection term

               gradN = gradPrev*(rho*Nv[gp][jj]);

               Db = rho*(velPrev[0]*dNvdx[gp][jj] + velPrev[1]*dNvdy[gp][jj]);

               gradN(0,0) += Db;
               gradN(1,1) += Db;

               Kuu(TI,   TJ)   += (b8*gradN(0,0));
               Kuu(TI,   TJp1) += (b8*gradN(0,1));
               Kuu(TIp1, TJ)   += (b8*gradN(1,0));
               Kuu(TIp1, TJp1) += (b8*gradN(1,1));
            }

            for(jj=0; jj<3; jj++)
            {
               // pressure term
               Kup(TI,   jj) -= (b1*af*Np[gp][jj]);
               Kup(TIp1, jj) -= (b2*af*Np[gp][jj]);
            }
        }

        //cout << velPrev[0] << '\t' << velPrev[1] << endl;

        for(ii=0; ii<3; ii++)
        {
            b8 = Np[gp][ii]*dvol;

            Fp(ii) += (b8*veloDiv);
        }

    }//gp

    //printMatrix(Kuu);
    //printVector(Fu);

    return 0;
}
//



int  BernsteinElem2DINSTriaP2bP1dc::StiffnessForSemiImpl(double* elemData, double* timeData, MatrixXd& Kup)
{
    // Semi-Implicit formulation

    int ii, jj, TI, TIp1;

    double  b1, b2, dvol;

    Kup.setZero();

    for(int gp=0; gp<nGP; gp++)
    {
        dvol = elemVolGP[gp];
        for(ii=0; ii<6; ii++)
        {
            TI   = 2*ii;
            TIp1 = TI+1;

            b1 = dNvdx[gp][ii]*dvol;
            b2 = dNvdy[gp][ii]*dvol;

            for(jj=0; jj<3; jj++)
            {
              // pressure term
              Kup(TI,   jj) -= (b1*Np[gp][jj]);
              Kup(TIp1, jj) -= (b2*Np[gp][jj]);
            }
          }
    }//gp

    //printVector(Flocal);

    return 0;
}



// Mass is assumed to be lumped.
// So, it is stored as a vector of diagonal vector
int BernsteinElem2DINSTriaP2bP1dc::MassMatrices(vector<myPoint>& nodeCoords, double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2)
{
    int  ii;

    //cout << " VOLUME = " << dvol << endl;

    // Mass Matrix for the Velocity DOF
    // elemData[0] --> density
    double fact = elemData[0]*elemVol/6.0;

    for(ii=0; ii<6; ii++)
    {
      Mlocal1[ii] = fact;
    }

    // Mass Matrix for the Pressure DOF
    fact = elemVol/3.0;

    for(ii=0; ii<3; ii++)
    {
      Mlocal2[ii] = fact;
    }

    return 0;
}


int BernsteinElem2DINSTriaP2bP1dc::calcLoadVector(VectorXd& Flocal)
{
    return 0;
}



//
double  BernsteinElem2DINSTriaP2bP1dc::CalculateError(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloCur, VectorXd& veloDotPrev, VectorXd& presCur, double timeCur, int index)
{
    double  b1, b2, b3, b4, fact, elemError;
    double  valNum[3], valExact[3], dp[2], gradNum[4], gradExact[4];

    int  nlbfU=6, gp, ii, jj, TI, TIp1, TIp2;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    KimMoinFlowUnsteadyNavierStokes  analy(rho, mu, 1.0);
    //Stokes2DEx1 analy;
    //KimMoinFlow  analy(rho, mu, 0.0);
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    elemError = 0.0;
    // L2 norm in x-velocity (index=0), y-velocity (index=1) and pressure (index=2)
    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of displacement first
          xx = 0.0; yy = 0.0;
          valNum[0] = valNum[1] = valNum[2] = 0.0;
          gradNum[0] = gradNum[1] = gradNum[2] = gradNum[3] = 0.0;
          for(ii=0; ii<6; ii++)
          {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloCur[TI];
            b2 = veloCur[TIp1];

            valNum[0]  +=  b1*Nv[gp][ii];
            valNum[1]  +=  b2*Nv[gp][ii];

            valNum[2]  +=  presCur[nodeNums[ii]]*Np[gp][ii];

            gradNum[0] += b1*dNvdx[gp][ii];
            gradNum[2] += b1*dNvdy[gp][ii];
            gradNum[1] += b2*dNvdx[gp][ii];
            gradNum[3] += b2*dNvdy[gp][ii];
          }

          //valNum[2] = 0.0;
          //for(ii=0; ii<3; ii++)
            //valNum[2] += veloPrev[nodeNums[ii]*3+2]*Np[gp][ii];

          if(index < 3)
          {
            valExact[index] = analy.computeValue(index, xx, yy, 0.0, timeCur);
            //if(index == 0)
              //valExact[index] = 4.0*yy*(1.0-yy);
            //else if(index == 1)
              //valExact[index] = 0.0;
            //else
              //valExact[index] = 8.0*mu*(2.0-xx);


            fact = valNum[index] - valExact[index];

            //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", valNum[index], valExact[index], fact);

            //if(index == 2)
              //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", valNum[2], valExact[2], fact);

            elemError += ( (fact*fact) * elemVolGP[gp] );
          }
          else
          {
            valExact[0] = analy.computeValue(0, xx, yy, 0.0, timeCur);
            valExact[1] = analy.computeValue(1, xx, yy, 0.0, timeCur);

            analy.computeDerivatives(xx, yy, 0.0, timeCur, gradExact);

            valNum[0] -= valExact[0];
            valNum[1] -= valExact[1];

            gradNum[0] -= gradExact[0];
            gradNum[1] -= gradExact[1];
            gradNum[2] -= gradExact[2];
            gradNum[3] -= gradExact[3];

            fact  = valNum[0]*valNum[0] + valNum[1]*valNum[1];
            fact += (gradNum[0]*gradNum[0]+gradNum[1]*gradNum[1]+gradNum[2]*gradNum[2]+gradNum[3]*gradNum[3]);

            elemError += ( fact * elemVolGP[gp] );
          }

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", valNum[index], valExact, fact);
    }

    return elemError;
}
//


/*
double  BernsteinElem2DINSTriaP2bP1dc::CalculateError(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloDotPrev, VectorXd& presPrev, double timeCur, int index)
{
    double  b1, b2, b3, b4, fact, elemError, valExact[3];
    double  valNum[3], dp[2], gradNum[4], gradExact[4];

    int  gp, ii, jj, TI, TIp1, TIp2;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    double xNode[npElem], yNode[npElem], xx, yy;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    PoissonAnnulus analy;
    //PoissonEx3 analy;
    //Stokes2DEx1 analy;
    //KimMoinFlow  analy(rho, mu, 0.0);
    //Kovasznay  analy;
    //analy.SetPressure(0.0);

    elemError = 0.0;
    // L2 norm in x-velocity (index=0), y-velocity (index=1) and pressure (index=2)
    for(gp=0; gp<nGP; gp++)
    {
          // compute the gradient of displacement first
          xx = 0.0; yy = 0.0;
          valNum[0] = valNum[1] = valNum[2] = 0.0;
          gradNum[0] = gradNum[1] = gradNum[2] = gradNum[3] = 0.0;
          for(ii=0; ii<6; ii++)
          {
            xx += xNode[ii]*Nv[gp][ii];
            yy += yNode[ii]*Nv[gp][ii];

            b1 = veloPrev[nodeNums[ii]];

            valNum[0]  +=  b1*Nv[gp][ii];

            gradNum[0] += b1*dNvdx[gp][ii];
            gradNum[1] += b1*dNvdy[gp][ii];
          }

          valExact[0] = analy.computeValue(index, xx, yy, timeCur);

          //printf(" computed, exact and difference \t %12.8f \t %12.8f \t %12.8f \n", valNum[0], valExact[0], fact);

          valNum[0] -= valExact[0];

          fact = valNum[0]*valNum[0];

          if(index > 0)
          {
            analy.computeDerivatives(xx, yy, gradExact);

            //printf(" computed, exact \t %12.8f \t %12.8f \t %12.8f \t %12.8f  \t %12.8f \t %12.8f \n", valNum[0], valExact[0], gradNum[0], gradNum[1], gradExact[0], gradExact[1]);

            gradNum[0] -= gradExact[0];
            gradNum[1] -= gradExact[1];

            fact  += ( gradNum[0]*gradNum[0]+gradNum[1]*gradNum[1] );
          }

          elemError += ( fact * elemVolGP[gp] );
    }

    cout << " elemError = " << elemError << endl;

    return elemError;
}
*/




int BernsteinElem2DINSTriaP2bP1dc::toComputeInfSupCondition(vector<myPoint>& nodeCoords, double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
{
//
//  printStiffnessMatrix();
//  printf("\n\n");
//  printForceVector();

    // to compute the inf-sup constant

    int  sizep=3, ii, jj, kk, ll, gp, TI, TIp1, TJ, TJp1;
    double  bb1, bb2, bb3, bb4, nsize=12;

    double  fact, fact1, dvol;

    // resize local matrices and initialise them to zero
    if(Kuu.rows() != nsize)
    {
      Kuu.resize(nsize, nsize);
      Kup.resize(nsize, sizep);
      Kpp.resize(sizep, sizep);
    }
    Kuu.setZero();
    Kup.setZero();
    Kpp.setZero();

    for(gp=0; gp<nGP; gp++)
    {
        dvol = elemVolGP[gp];

        // compute Kuu and Kup
        for(ii=0;ii<nlbfU;ii++)
        {
          bb1 = dvol*dNvdx[gp][ii];
          bb2 = dvol*dNvdy[gp][ii];

          TI   = 2*ii;
          TIp1 = TI+1;

          // compute Kuu
          for(jj=0; jj<nlbfU; jj++)
          {
            TJ   = 2*jj;
            TJp1 = TJ+1;

            fact = bb1*dNvdx[gp][jj] + bb2*dNvdy[gp][jj] ;

            Kuu(TI,   TJ)   +=  fact ;
            Kuu(TIp1, TJp1) +=  fact ;
          }

          // compute Kup
          for(jj=0; jj<sizep; jj++)
          {
            Kup(TI,   jj) += bb1*Np[gp][jj];
            Kup(TIp1, jj) += bb2*Np[gp][jj];
          }
        }

        for(ii=0;ii<sizep;ii++)
        {
          bb3 = dvol*Np[gp][ii];

          // compute Kup
          for(jj=0; jj<sizep; jj++)
          {
            Kpp(ii, jj) += bb3*Np[gp][jj];
          }
        }
    }//gp

    return 0;
}




int  BernsteinElem2DINSTriaP2bP1dc::CalculateForces(int side, vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& presVec, VectorXd&  Flocal1)
{
    // side 0 - 1-4-2
    // side 1 - 2-5-3
    // side 2 - 3-6-1
    int face_node_map[3][3]={{0,3,1},{1,4,2},{2,5,0}};


    double  Jac, dvol, b1, b2, b3, b4, pres;
    double  grad[2][2], stress[2][2], trac[2], normal[2], tarpoint[2], Ne[3];
    double  xx, yy, fact, fact2;
    double  xNode[6], yNode[6], xNodeEdge[3], yNodeEdge[3], param[3];

    VectorXd  Nv(6), dNvdx(6), dNvdy(6), Np(3);
    Nv.setZero();    dNvdx.setZero();    dNvdy.setZero();
    Np.setZero();

    int  gp, ii, jj, TI, TIp1, degree=2;

    // material parameters
    double  rho = elemData[0];
    double  mu  = elemData[1];

    for(ii=0; ii<6; ii++)
    {
      xNode[ii] = nodeCoords[nodeNums[ii]][0];
      yNode[ii] = nodeCoords[nodeNums[ii]][1];
    }

    for(ii=0; ii<3; ii++)
    {
      xNodeEdge[ii] = xNode[face_node_map[side][ii]];
      yNodeEdge[ii] = yNode[face_node_map[side][ii]];
    }


    vector<double>  gausspoints(3), gaussweights(3);

    gausspoints[0] = -0.774596669241483;  gaussweights[0] = 5.0/9.0;
    gausspoints[1] =  0.0;                gaussweights[1] = 8.0/9.0;
    gausspoints[2] =  0.774596669241483;  gaussweights[2] = 5.0/9.0;

    //loop over Gauss points and compute element residual

    for(gp=0; gp<3; gp++)
    {
        // compute the normal and the Jacobian
        param[0] = gausspoints[gp];

        //BernsteinBasisFunsEdge2D(degree, param, xNodeEdge, yNodeEdge, Ne, normal, Jac);

        dvol = gaussweights[gp]*Jac;

        //compute the basis functions of calculating the stress

        if(side == 0)
        {
          param[0] = (1.0+param[0])/2.0;
          param[1] = 0.0;
        }
        else if(side == 1)
        {
          param[0] = 0.5;
          param[1] = 0.5;
        }
        else
        {
          param[0] = 0.0;
          param[1] = (1.0+param[0])/2.0;
        }

        tarpoint[0] = tarpoint[1] = 0.0;
        for(ii=0; ii<3; ii++)
        {
          tarpoint[0] += xNodeEdge[ii]*Ne[ii];
          tarpoint[1] += yNodeEdge[ii]*Ne[ii];
        }

        pointInverseTria6node(xNode, yNode, tarpoint, param);

        computeBasisFunctions2D(false, ELEM_TYPE, degree, param, xNode, yNode, &Nv[0], &dNvdx[0], &dNvdy[0], Jac);

        Np[0] = 1.0-param[0]-param[1];  Np[1] = param[0];  Np[2] = param[1];

        // calculate the velocity gradient

        grad[0][0] = grad[0][1] = grad[1][0] = grad[1][1] = 0.0;
        for(ii=0; ii<6; ii++)
        {
            TI   = nodeNums[ii]*2;
            TIp1 = TI+1;

            b1 = veloVec[TI];
            b2 = veloVec[TIp1];

            grad[0][0] += b1*dNvdx[ii];
            grad[0][1] += b1*dNvdy[ii];
            grad[1][0] += b2*dNvdx[ii];
            grad[1][1] += b2*dNvdy[ii];
        }

        pres = 0.0;
        for(ii=0; ii<3; ii++)
        {
          pres += presVec[nodeNums[ii]]*Np[ii];
        }

          fact = 2.0*mu*(grad[0][0]+grad[1][1])/3.0;
          stress[0][0] = mu*(grad[0][0]+grad[0][0]) - fact - pres;
          stress[0][1] = mu*(grad[0][1]+grad[1][0]);
          stress[1][0] = mu*(grad[1][0]+grad[0][1]);
          stress[1][1] = mu*(grad[1][1]+grad[1][1]) - fact - pres;

          //stress[0][0] = mu*grad[0][0] - pres;
          //stress[0][1] = mu*grad[0][1] ;
          //stress[1][0] = mu*grad[1][0] ;
          //stress[1][1] = mu*grad[1][1] - pres;

          trac[0] = stress[0][0]*normal[0] + stress[0][1]*normal[1];
          trac[1] = stress[1][0]*normal[0] + stress[1][1]*normal[1];

          for(ii=0; ii<3; ii++)
          {
            b4 = Ne[ii]*dvol;

            Flocal1[0] += (b4*trac[0]);
            Flocal1[1] += (b4*trac[1]);
          }
    }

    return 0;
}




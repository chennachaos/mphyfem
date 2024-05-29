
#include "MagnetoMech3DTet21.h"
#include "MyTime.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "QuadratureUtil.h"
#include "TimeFunction.h"
#include "stressrecovery.h"
#include "BasisFunctionsBernstein.h"
#include "growthfunctions.h"

using namespace std;

extern MyTime myTime;
extern List<TimeFunction> timeFunction;



MagnetoMech3DTet21::MagnetoMech3DTet21()
{
  ELEM_SHAPE = ELEM_SHAPE_TETRA_BERNSTEIN;

  ndim   = 3;
  degree = 2;
  npElem = 10;
  nlbfU  = 10;
  nlbfP  = 4;
  ndof   = 3;
  nsize  = nlbfU*ndof;
}



MagnetoMech3DTet21::~MagnetoMech3DTet21()
{
}


void MagnetoMech3DTet21::prepareElemData()
{
  ElementBase::prepareElemData();

  return;
}



double MagnetoMech3DTet21::computeVolume(bool configflag)
{
    double  Jac, param[3];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    elemVol=0.0;
    for(int gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        if(configflag)
          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);
        else
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        elemVol += (gaussweights[gp]*Jac);
    }//gp

    return elemVol;
}




/*
int MagnetoMech3DTet21::calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != nsize) || (Mlocal.cols() != nsize) )
      Mlocal.resize(nsize, nsize);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = elmDat[5]*computeVolume(true)/27.0;

      for(int ii=0; ii<nsize; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
    else
    {
      int  ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;
      double  fact, dvol0, Jac, bb1, cc1, param[3];
      double  rho0  = elmDat[5] ;

      VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

      vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
      int nGPt=27;
      getGaussPointsTetra(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
  
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
*/



int MagnetoMech3DTet21::calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
{
    if( (Mlocal.rows() != 8) || (Mlocal.cols() != 8) )
      Mlocal.resize(8, 8);
    Mlocal.setZero();

    // mass lumping - row-wise sum
    if(MassLumping)
    {
      double fact = computeVolume(true)/8.0;

      for(int ii=0; ii<4; ii++)
      {
        Mlocal(ii,ii) = fact;
      }
    }
    else
    {
      int  ii, jj, gp, TI, TIp1, TIp2, TJ, TJp1, TJp2;
      double  fact, dvol0, Jac, bb1, cc1, param[3];
      double  rho0  = elmDat[5] ;

      VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

      vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
      int nGPt=27;
      getGaussPointsTetra(nGPt, gausspoints1, gausspoints2, gausspoints3, gaussweights);
  
      for(gp=0; gp<nGPt; gp++)
      {
          param[0] = gausspoints1[gp];
          param[1] = gausspoints2[gp];
          param[2] = gausspoints3[gp];
  
          GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, 1, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

          dvol0 = gaussweights[gp]*Jac;
  
          for(ii=0;ii<4;ii++)
          {
            bb1 = dvol0*N[ii];
  
            for(jj=0; jj<4; jj++)
            {  
              Mlocal(ii, jj)  += (bb1*N[jj]) ;
            }
          }
      } //gp
      //cout << " elemVol = " << elemVol << endl;
      //printMatrix(Klocal);  printf("\n\n\n");
    }

    return 0;
}
  


double  MagnetoMech3DTet21::calcCriticalTimeStep(bool flag)
{
    int  ii, jj, gp;
/*
    double  F[9], detF, fact, dvol, dvol0, Jac,  param[3];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double BULK = matDat[0];
    double eps  = 0.0; if(!SolnData->TRULY_INCOMPRESSIBLE) eps=1.0/BULK;
    double  mu  = matDat[1], mu_effective = -1.0 ;
    double  C10  = matDat[1], C20 = matDat[2], C01 = matDat[3] ;
    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    double  wave_speed=-1.0, dtCrit, psi1, psi2, Ib;

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = Jac * gaussweights[gp];
        dvol = dvol0;

        F[0] = computeValueCur(0, dN_dx) + 1.0;
        F[3] = computeValueCur(0, dN_dy);
        F[6] = computeValueCur(0, dN_dz);

        F[1] = computeValueCur(1, dN_dx);
        F[4] = computeValueCur(1, dN_dy) + 1.0;
        F[7] = computeValueCur(1, dN_dz);

        F[2] = computeValueCur(2, dN_dx);
        F[5] = computeValueCur(2, dN_dy);
        F[8] = computeValueCur(2, dN_dz) + 1.0;

        detF = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[3]*(F[1]*F[8]-F[2]*F[7]) + F[6]*(F[1]*F[5]-F[2]*F[4]);

        Ib = (F[0]*F[0]+F[3]*F[3]+F[6]*F[6]) + (F[1]*F[1]+F[4]*F[4]+F[7]*F[7]) + (F[2]*F[2]+F[5]*F[5]+F[8]*F[8]);

        // Neo-Hookean
        psi1 = 0.5*mu;
        psi2 = 0.0;

        // Isihara
        psi1 = C10;
        psi2 = 2.0*C20*(Ib-3.0);

        mu_effective = 2.0*(psi1 + Ib*psi2);
        rho = rho0/detF;

        fact = mu_effective/rho0;

        wave_speed = max(wave_speed, fact);
    }
    wave_speed = sqrt(wave_speed);
*/

    double xNode[10], yNode[10], zNode[10];
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosNew[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosNew[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosNew[nodeNums[ii]][2];
    }

    double  xd, yd, zd, charlen=1.0e10;
    int  edge_nodes[6][2] = { {1,0}, {2,0}, {3,0}, {2,1}, {3,1}, {3,2}};

    for(ii=0; ii<6; ii++)
    {
      xd = xNode[edge_nodes[ii][0]] - xNode[edge_nodes[ii][1]];
      yd = yNode[edge_nodes[ii][0]] - yNode[edge_nodes[ii][1]]; 
      zd = zNode[edge_nodes[ii][0]] - zNode[edge_nodes[ii][1]];

      charlen = min(charlen, xd*xd+yd*yd+zd*zd);
    }
    charlen = 0.5*sqrt(charlen);

    //cout << " elemVolCur = " << elemVolCur << endl;
    //charlen = 0.5*pow(6.0*abs(elemVolCur)/PI, 1.0/3.0);
    double BULK = matDat[0];
    double  mu  = matDat[1] ;
    //double  mu  = matDat[6] ;
    //double  mu  = 0.5*(matDat[2]+matDat[3]+matDat[4]);
    double rho0 = elmDat[5] ;

    double  wave_speed;
    if(flag)
      wave_speed = sqrt( mu/rho0 );
    else
      wave_speed = sqrt((BULK+4.0*mu/3.0)/rho0);

    return  (charlen/wave_speed);
}





int MagnetoMech3DTet21::calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter)
{
    int   err,  isw,  count,  count1, index, ii, jj, kk, ll, mm, gp;
    int   TI, TIp1, TIp2, TJ, TJp1, TJp2;

    double  F[9], detF, fact, fact1, fact2, fact3, dvol, dvol0, Jac, volstrain, pbar;
    double  bb1, bb2, bb3, bb4, bb5, dUdJ, d2UdJ2, Rp;
    double  stre[6], cc[6][6], param[3], bforce[3], force[3], sig[3], Np[4];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double BULK = matDat[0];
    double eps  = 0.0; if(!SolnData->TRULY_INCOMPRESSIBLE) eps=1.0/BULK;
    double  mu  = matDat[1] ;
    double LAMBDA = BULK-2.0*mu/3.0;
    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    bforce[2]   = elmDat[8]*timeFunction[0].prop ;
    double  dt  = myTime.dt;
    double tPrev = myTime.cur - dt;

    double xNode[10], yNode[10], zNode[10], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal1.rows() < nsize)   Flocal1.resize(nsize);
    Flocal1.setZero();
    if(Flocal2.rows() < 4)   Flocal2.resize(4);
    Flocal2.setZero();
    if(Kpu.rows() < 4 || Kpu.cols() < nsize)   Kpu.resize(4, nsize);
    Kpu.setZero();
    if(Kup.rows() < nsize || Kup.cols() < 4)   Kup.resize(nsize, 4);
    Kup.setZero();
    if(Kpp.rows() < 4 || Kpp.cols() < 4)   Kpp.resize(4, 4);
    Kpp.setZero();

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

        xx = yy = zz = 0.0;
        for(ii=0;ii<nlbfU;ii++)
        {
          xx += xNode[ii]*N[ii];
          yy += yNode[ii]*N[ii];
          zz += zNode[ii]*N[ii];
        }


        F[0] = computeValuePrev(0, dN_dx) + 1.0;
        F[3] = computeValuePrev(0, dN_dy);
        F[6] = computeValuePrev(0, dN_dz);

        F[1] = computeValuePrev(1, dN_dx);
        F[4] = computeValuePrev(1, dN_dy) + 1.0;
        F[7] = computeValuePrev(1, dN_dz);

        F[2] = computeValuePrev(2, dN_dx);
        F[5] = computeValuePrev(2, dN_dy);
        F[8] = computeValuePrev(2, dN_dz) + 1.0;

        volstrain = F[0] + F[4] + F[8] - 3.0;

        detF = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[3]*(F[1]*F[8]-F[2]*F[7]) + F[6]*(F[1]*F[5]-F[2]*F[4]);

        if(detF < 0.0)
        {
          cerr << " Negative Jacobian in the element " << elenum << "'\t' at quadrature point " << gp << endl;
          return 1;
        }

        //matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0),  &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }
        elemVolCur  += dvol;

        // basis functions for pressure dof
        Np[0] = 1.0-param[0]-param[1]-param[2];    Np[1] = param[0];    Np[2] = param[1];    Np[3] = param[2];

        // evaluate pressure at the quadrature points
        pres = 0.0;
        for(ii=0; ii<4; ii++)
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);


        if(!finite)
        {
          d2UdJ2 = 1.0;

          Rp = volstrain-pres*eps;
        }

        if(SolnData->TRULY_INCOMPRESSIBLE)
        {
          stre[0] += pres;
          stre[1] += pres;
          stre[2] += pres;

          if(finite)
          {
            d2UdJ2 = 1.0;

            Rp = detF-1.0-pres*eps;
          }
        }
        else
        {
          fact = pres - (stre[0]+stre[1]+stre[2])/3.0;

          stre[0] += fact;
          stre[1] += fact;
          stre[2] += fact;

          if(finite)
          {
            if( (int) matDat[2] == 6 )
            {
              dUdJ   = LAMBDA*log(detF)/detF;
              d2UdJ2 = LAMBDA*(1.0-log(detF))/detF/detF;
            }
            else
            {
              dUdJ   = 0.5*BULK*(detF-1.0/detF);
              d2UdJ2 = 0.5*BULK*(1.0+1.0/detF/detF);
            }
            d2UdJ2 /= BULK;

            Rp = (dUdJ-pres)*eps;
          }
        }

        //printf("%12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \t %12.6f \n", stre[0], stre[1], stre[2], stre[3], stre[4], stre[5]);

        // Calculate Residual
        //===================

        fact = -1.0e7;
        double r1 = sqrt(xx*xx+yy*yy);
        double r2 = sqrt(xx*xx+yy*yy+zz*zz);
        double theta = atan(yy/xx);
        if( xx < 0.0 )
          theta += PI;
        if( CompareDoubles(xx,0.0) && CompareDoubles(yy,0.0) )
          theta = 0.0;
        double phi   = atan(r1/zz);
        if( zz < 0.0)
          phi += PI;

        //bforce[0]  =  fact*cos(theta)*sin(phi);
        //bforce[1]  =  fact*sin(theta)*sin(phi);
        //bforce[2]  =  fact*cos(phi);

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
          bb5 = dvol0*N[ii];

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          Flocal1(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3] - bb3*stre[5]) ;
          Flocal1(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1] - bb3*stre[4]) ;
          Flocal1(TIp2) += (bb5*force[2] - bb1*stre[5] - bb2*stre[4] - bb3*stre[2]) ;

          for(jj=0; jj<4; jj++)
          {
            Kup(TI,   jj)  += (bb1*Np[jj]);
            Kup(TIp1, jj)  += (bb2*Np[jj]);
            Kup(TIp2, jj)  += (bb3*Np[jj]);
          }
        }

        for(ii=0; ii<4; ii++)
        {
          bb4 = dvol*Np[ii]*d2UdJ2;
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

          bb5 *= eps;
          for(jj=0; jj<4; jj++)
          {
            Kpp(ii, jj)    -= bb5*Np[jj];
          }
        }
    }//gp

    //printMatrix(Kup);      printf("\n\n\n");
    //printMatrix(Kpu);      printf("\n\n\n");
    //printMatrix(Kpp);      printf("\n\n\n");
    //printVector(Flocal1);  printf("\n\n\n");
    //printVector(Flocal2);  printf("\n\n\n");

    return 0;
}




int MagnetoMech3DTet21::calcResidual(VectorXd& Flocal)
{
    int   err,  isw,  count,  count1, index, ii, jj, ll, gp, TI, TIp1, TIp2;

    double  F[9], detF, fact, dvol, dvol0, Jac, volstrain;
    double  bb1, bb2, bb3, bb4, bb5;
    double  stre[6], cc[6][6], param[3], bforce[3], force[3], Np[4];

    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double BULK = matDat[0] ;
    double eps  = 1.0/BULK;
    double rho0 = elmDat[5] ;
    double rho  = rho0 ;
    bforce[0]   = elmDat[6]*timeFunction[0].prop ;
    bforce[1]   = elmDat[7]*timeFunction[0].prop ;
    bforce[2]   = elmDat[8]*timeFunction[0].prop ;
    double  dt  = myTime.dt;

    double xNode[10], yNode[10], zNode[10], xx, yy, zz;
    for(ii=0;ii<npElem;ii++)
    {
      xNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][0];
      yNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][1];
      zNode[ii] = GeomData->NodePosOrig[nodeNums[ii]][2];
    }

    //  compute tr(F) and det(F) by averaging the shape function derivatives
    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    if(Flocal.rows() != nsize)
      Flocal.resize(nsize);
    Flocal.setZero();

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
        dvol  = dvol0;
        elemVolOrig += dvol0;

        F[0] = computeValuePrev(0, dN_dx) + 1.0;
        F[3] = computeValuePrev(0, dN_dy);
        F[6] = computeValuePrev(0, dN_dz);

        F[1] = computeValuePrev(1, dN_dx);
        F[4] = computeValuePrev(1, dN_dy) + 1.0;
        F[7] = computeValuePrev(1, dN_dz);

        F[2] = computeValuePrev(2, dN_dx);
        F[5] = computeValuePrev(2, dN_dy);
        F[8] = computeValuePrev(2, dN_dz) + 1.0;

        volstrain = F[0] + F[4] + F[8] - 3.0;

        detF = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[3]*(F[1]*F[8]-F[2]*F[7]) + F[6]*(F[1]*F[5]-F[2]*F[4]);

        if(detF < 0.0)
        {
          cerr << " Negative Jacobian in the element " << elenum << "'\t' at quadrature point " << gp << endl;
          return 1;
        }

        //matlib3d_(matDat, F, stre, cc[0], &(intVar1[ll]), &(intVar2[ll]), &dt, &matId, &nivGP, &finiteInt, &isw, &err, &count, NULL);
        count++;
        ll += nivGP;

        if(err !=0)    return 1;

        if(finite)
        {
          GeomData->computeBasisFunctions3D(CONFIG_DEFORMED, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0),  &dN_dz(0), Jac);
          dvol = gaussweights[gp]*Jac;
        }
        elemVolCur  += dvol;

        // basis functions for pressure dof
        Np[0] = 1.0-param[0]-param[1]-param[2];    Np[1] = param[0];    Np[2] = param[1];    Np[3] = param[2];

        // evaluate pressure at the quadrature points
        pres = 0.0;
        for(ii=0; ii<4; ii++)
          pres += (Np[ii]*SolnData->var2[nodeNums[ii]]);

        fact = pres - (stre[0]+stre[1]+stre[2])/3.0;

        stre[0] += fact;
        stre[1] += fact;
        stre[2] += fact;

        // Calculate Residual
        //===================

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
          bb5 = dvol0*N[ii];

          TI   = 3*ii;
          TIp1 = TI+1;
          TIp2 = TI+2;

          Flocal(TI)   += (bb5*force[0] - bb1*stre[0] - bb2*stre[3] - bb3*stre[5]) ;
          Flocal(TIp1) += (bb5*force[1] - bb1*stre[3] - bb2*stre[1] - bb3*stre[4]) ;
          Flocal(TIp2) += (bb5*force[2] - bb1*stre[5] - bb2*stre[4] - bb3*stre[2]) ;
        }
    }//gp

    return 0;
}



int MagnetoMech3DTet21::calcResidualPressure(VectorXd& Flocal)
{
    int gp, ii;
    VectorXd  N(nlbfU), dN_dx(nlbfU), dN_dy(nlbfU), dN_dz(nlbfU);

    double  F[9], detF, volstrain, param[3], dvol0, dvol;
    double  fact, Jac, r2d3 = 2.0/3.0, dUdJ, Np[4];
    double  BULK = matDat[0];
    double mu   = matDat[1];
    double LAMBDA = BULK-r2d3*mu;

    vector<double>  gausspoints1, gausspoints2, gausspoints3, gaussweights;
    getGaussPointsTetra(nGP, gausspoints1, gausspoints2, gausspoints3, gaussweights);

    // resize local matrices and initialise them to zero
    if(Flocal.rows() != 4)
      Flocal.resize(4);
    Flocal.setZero();

    elemVolOrig = 0.0;
    for(gp=0; gp<nGP; gp++)
    {
        param[0] = gausspoints1[gp];
        param[1] = gausspoints2[gp];
        param[2] = gausspoints3[gp];

        GeomData->computeBasisFunctions3D(CONFIG_ORIGINAL, ELEM_SHAPE, degree, param, nodeNums, &N(0), &dN_dx(0), &dN_dy(0), &dN_dz(0), Jac);

        dvol0 = gaussweights[gp]*Jac;
        elemVolOrig  += dvol0;

        if(finite)
        {
          F[0] = computeValue(0, dN_dx) + 1.0;
          F[3] = computeValue(0, dN_dy);
          F[6] = computeValue(0, dN_dz);

          F[1] = computeValue(1, dN_dx);
          F[4] = computeValue(1, dN_dy) + 1.0;
          F[7] = computeValue(1, dN_dz);

          F[2] = computeValue(2, dN_dx);
          F[5] = computeValue(2, dN_dy);
          F[8] = computeValue(2, dN_dz) + 1.0;

          detF = F[0]*(F[4]*F[8]-F[5]*F[7]) - F[3]*(F[1]*F[8]-F[2]*F[7]) + F[6]*(F[1]*F[5]-F[2]*F[4]);

          if( matDat[2] == 6 )
          {
            dUdJ   = LAMBDA*log(detF)/detF;
          }
          else
          {
            dUdJ   = 0.5*BULK*(detF-1.0/detF);
          }
          fact = dUdJ;
          //fact = BULK*(detF-1.0);
        }
        else
        {
          volstrain = computeValue(0, dN_dx) + computeValue(1, dN_dy) + computeValue(2, dN_dz);

          fact = BULK*volstrain;
        }

        // basis functions for pressure dof
        Np[0] = 1.0-param[0]-param[1]-param[2];    Np[1] = param[0];    Np[2] = param[1];    Np[3] = param[2];

        fact *= dvol0;
        for(ii=0; ii<4; ii++)
        {
          Flocal(ii) += (fact*Np[ii]);
        }
    }

    return 0;
}



int MagnetoMech3DTet21::calcLoadVector(VectorXd& Flocal)
{
  return 0;
}



int MagnetoMech3DTet21::calcInternalForces()
{
  return 0;
}


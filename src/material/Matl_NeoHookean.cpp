#include "Matl_NeoHookean.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"

using namespace std;


Matl_NeoHookean::Matl_NeoHookean()
{

}


Matl_NeoHookean::~Matl_NeoHookean()
{

}


double Matl_NeoHookean::computeValue(int sss,  MatrixXd&  F)
{
    double  Psi, r2d3 = 2.0/3.0;

    double  K      = 1.0/Kinv;                              // bulk modulus
    double  mu     = matData[0];                            // shear modulus

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jem2d3 = pow(J, -r2d3);

    /////////////////////////////////
    // Energy function value

    Psi = mu*(Jem2d3*b.trace()-3.0) + 0.5*K*(J-1.0)*(J-1.0);

    return Psi;
}



int Matl_NeoHookean::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = mu/2 (J^{-2/3} trb - 3) + U(J)
    // U(J) is selected based on the input value for 'Utype'

    int err =  computeStressAndTangent_NeoHooke(sss, MIXED_ELEMENT, Utype, Kinv, data_deviatoric, Fn, F, pres, stre, Cmat);

    return err;
}




int Matl_NeoHookean::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    int     Utype  = int(matData[0]); // volumetric energy function
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    double  mu     = matData[2];      // shear modulus

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jem2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jem2d3/J;
    double  Ib = b(0,0) + b(1,1) + b(2,2);

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);

    double  fact1 = mu*Jm5d3;
    double  fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre.setZero();

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    return 0;
}




int Matl_NeoHookean::getStressAndJacobianForInverseGrowth(int sss,  MatrixXd&  F, MatrixXd&  Fg, VectorXd&  stre, MatrixXd&  Jgp, MatrixXd&  Hgp)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r2d9 = 2.0/9.0;

    double  mu     = data_deviatoric[0];      // shear modulus
    double  BULK   = 1.0/Kinv;

    Matrix3d  FInv  = F.inverse();
    Matrix3d  FgInv = Fg.inverse();
    //printMatrix(FgInv);
    Matrix3d  Fe    = F*FgInv;
    //cout << " Fe " << endl;
    //printMatrix(Fe);
    //cout << Fe << endl;
    Matrix3d  FeInv = Fe.inverse();
    //printMatrix(FeInv);
    //cout << FeInv << endl;
    Matrix3d  Ce    = Fe.transpose()*Fe;

    double  Jg    = Fg.determinant();
    double  Je    = Fe.determinant();
    double  ICe   = Ce.trace();
    double  Jem2d3 = pow(Je, -r2d3), value;

    // compressible model
    // MatrixXd  PK1e = mu*(Fe - FeInv) + K*Je*(Je-1)*FeInv;

    // incompressible model
    Matrix3d  PK1e = mu*Jem2d3*(Fe - ICe*FeInv/3.0) + BULK*Je*(Je-1)*FeInv;

    Matrix3d  PK1  = Jg*PK1e*FgInv;

    stre[0] = PK1(0,0);    stre[3] = PK1(0,1);    stre[6] = PK1(0,2);
    stre[1] = PK1(1,0);    stre[4] = PK1(1,1);    stre[7] = PK1(1,2);
    stre[2] = PK1(2,0);    stre[5] = PK1(2,1);    stre[8] = PK1(2,2);


    int  map_tensor4[9][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  ii, jj, kk, i, J, iJ, K, L, KL, MN, M, N, I, p, Q, j, k, l, U, V, r, s;

    //cout << " 9999999999999 " << endl;

    //compute the fourthorder tensor
    // incompressible model
    double fact1 = mu*Jem2d3;
    double fact2 = BULK*Je*(Je-1.0);
    double fact3 = BULK*Je*(2*Je-1.0);

    MatrixXd  Cmat(9,9); Cmat.setZero();

    for(ii=0; ii<9; ii++)
    {
        i = map_tensor4[ii][0];
        j = map_tensor4[ii][1];

        for(jj=0; jj<9; jj++)
        {
          k = map_tensor4[jj][0];
          l = map_tensor4[jj][1];

          // incompressible model
          Cmat(ii,jj) +=      fact1*KronDelta(i,k)*KronDelta(j,l);
          Cmat(ii,jj) -= r2d3*fact1*FeInv(l,k)*(Fe(i,j)-r1d3*ICe*FeInv(j,i));
          Cmat(ii,jj) -= r2d3*fact1*Fe(k,l)*FeInv(j,i);
          Cmat(ii,jj) += r1d3*fact1*ICe*FeInv(j,k)*FeInv(l,i);

          Cmat(ii,jj) -= fact2*FeInv(j,k)*FeInv(l,i);
          Cmat(ii,jj) += fact3*FeInv(j,i)*FeInv(l,k);
        }
    }

    //cout << " 9999999999999 " << endl;

    int invmap[3][3] = { {0,3,6}, {1,4,7}, {2,5,8} };


    Jgp.setZero();

    int  tvec[3] = {0,4,8};
    for(iJ=0; iJ<9; iJ++)
    {
        i = map_tensor4[iJ][0];
        J = map_tensor4[iJ][1];

        for(kk=0; kk<3; kk++)
        {
          KL = tvec[kk];
          K = map_tensor4[KL][0];
          L = map_tensor4[KL][1];

          value = 0.0;
          for(I=0; I<3; I++)
          {
            value += PK1e(i,I)*FgInv(I,J)*FgInv(L,K);
            value -= PK1e(i,I)*FgInv(I,K)*FgInv(L,J);
            for(p=0; p<3; p++)
            {
              for(Q=0; Q<3; Q++)
              {
                value -= ( Cmat(invmap[i][I],invmap[p][Q])*Fe(p,K)*FgInv(Q,L)*FgInv(I,J) );
              }
            }
          }

          Jgp(iJ,kk) = Jg*value;
        }
    }

    double muJem2d3 = mu*Jem2d3;

    Hgp.setZero();

    for(r=0; r<3; r++)
    {
      KL = tvec[r];
      K = map_tensor4[KL][0];
      L = map_tensor4[KL][1];

      for(s=0; s<3; s++)
      {
        MN = tvec[s];
        M = map_tensor4[MN][0];
        N = map_tensor4[MN][1];

        value = 0.0;

        for(i=0; i<3; i++)
        {
          for(J=0; J<3; J++)
          {
            //
            // volumetric part
            fact1 = PK1(i,J)*BULK*Je*Jg*FInv(L,i);
            value -= fact1*(2*Je-1.0)*FgInv(K,J)*FgInv(N,M);
            value += fact1*  (Je-1.0)*FgInv(M,J)*FgInv(N,K);

            fact1 = PK1(i,J)*BULK*Je*Jg*FInv(N,i);
            value += fact1*  (Je-1.0)*FgInv(M,L)*FgInv(K,J);
            value -= fact1*(2*Je-1.0)*FgInv(K,L)*FgInv(M,J);

            fact1 = PK1(i,J)*BULK*Je*Jg*FInv(J,i);
            value += fact1*(4*Je-1.0)*FgInv(K,L)*FgInv(N,M);
            value -= fact1*(2*Je-1.0)*FgInv(M,L)*FgInv(N,K);
            //

            value += PK1(i,J)*( PK1(i,J)*FgInv(N,M)*FgInv(L,K) + PK1(i,M)*FgInv(N,K)*FgInv(L,J) + PK1(i,K)*FgInv(L,M)*FgInv(N,J) );
            value -= PK1(i,J)*( PK1(i,K)*FgInv(N,M)*FgInv(L,J) + PK1(i,J)*FgInv(L,M)*FgInv(N,K) + PK1(i,M)*FgInv(L,K)*FgInv(N,J) );

            for(I=0; I<3; I++)
            {
              //value += PK1(i,J)*( Jg*Cmat(invmap[i][I],invmap[M][N])*(FgInv(L,K)*FgInv(I,J) - FgInv(I,K)*FgInv(L,J)) );

              /*
              for(p=0; p<3; p++)
              {
                for(Q=0; Q<3; Q++)
                {
                  fact1 = PK1(i,J)*Jg*Cmat(invmap[i][I],invmap[p][Q]);

                  value += (fact1*Fe(p,M)*FgInv(N,K)*FgInv(Q,L)*FgInv(I,J) );
                  value += (fact1*Fe(p,K)*FgInv(Q,M)*FgInv(N,L)*FgInv(I,J) );
                  value += (fact1*Fe(p,K)*FgInv(Q,L)*FgInv(I,M)*FgInv(N,J) );
                  value -= (fact1*Fe(p,K)*FgInv(N,M)*FgInv(Q,L)*FgInv(I,J) );
                }
              }
              */

              // Part 1
              fact1 = -r2d3*muJem2d3*Jg*Fe(p,K)*FgInv(N,M);

              fact2  = KronDelta(i,p)*FgInv(L,I)*FgInv(I,J) - r2d3*(FInv(L,p)*Fe(i,I)*FgInv(I,J)+FInv(J,i)*Fe(p,I)*FgInv(I,L));
              fact2 += r1d3*ICe*(FInv(J,p)*FInv(L,i)+r2d3*FInv(L,p)*FInv(J,i));

              value += PK1(i,J)*fact1*fact2;

              // Part 2
              fact1  = -r2d3*muJem2d3*Jg;
              fact2  = -FgInv(L,M)*FgInv(I,J)*FgInv(N,K)*Fe(i,I) + FgInv(L,K)*FgInv(I,J)*FgInv(N,I)*Fe(i,M);
              fact2 +=  FgInv(L,I)*FgInv(J,i)*FgInv(N,I)*Ce(K,M) - FgInv(L,I)*FgInv(M,J)*FgInv(N,i)*Ce(K,I);

              value += PK1(i,J)*fact1*fact2;

              // Part 2
              fact1  = r2d3*muJem2d3*Jg*Ce(I,M)*FeInv(N,I);
              fact2  = FInv(L,i)*FgInv(K,J) + r2d3*FgInv(L,K)*FInv(J,i);

              value += PK1(i,J)*fact1*fact2;

            } //I

            // Part 2
            fact1 = -r1d3*muJem2d3*ICe*Jg;
            fact2   = FInv(L,i)*FgInv(M,J)*FgInv(N,K) + FInv(N,i)*FgInv(L,M)*FgInv(K,J);
            fact2  += r2d3*(FInv(J,i)*FgInv(L,M)*FgInv(N,K) + FInv(N,i)*FgInv(L,K)*FgInv(M,J));
              
            value += PK1(i,J)*fact1*fact2;
          }
        }

        Hgp(r,s) = value;
      }
    }


    return 0;
}




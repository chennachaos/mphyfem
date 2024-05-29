#include "Matl_SaintVenantKirchhoff.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"


using namespace std;


Matl_SaintVenantKirchhoff::Matl_SaintVenantKirchhoff()
{

}



Matl_SaintVenantKirchhoff::~Matl_SaintVenantKirchhoff()
{

}



int Matl_SaintVenantKirchhoff::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = mu/2 (J^{-2/3} trb - 3) + U(J)
    // U(J) is selected based on the input value for 'Utype'

    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r5d3 = 5.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  mu     = matData[0];      // shear modulus
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    int     Utype  = int(matData[2]); // volumetric energy function

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b(0,0) + b(1,1) + b(2,2);

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    stre.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    fact1 = mu*Jm5d3;
    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    /////////////////////////////////
    // material tangent tensor

        // without stress terms
        fact1 =  r1d3*mu*Jm5d3*Ib;
        fact2 = -r2d3*mu*Jm5d3;
        fact3 =  r2d3*fact1;

        // with stress terms
        //fact1 =  r1d3*mu*Jm5d3*Ib;
        //fact2 = -r2d3*mu*Jm5d3;
        //fact3 =  r2d3*r1d3*mu*Jm5d3*Ib;
        //fact4 =  mu*Jm5d3;

        fact7 = J*d2UdJ2 + dUdJ;
        fact8 = -dUdJ;
        fact9 = -dUdJ;

        for(i1=0; i1<9; i1++)
        {
            i = map_tensor4[i1][0];
            j = map_tensor4[i1][1];

            dij = KronDelta(i,j);

            for(j1=0; j1<9; j1++)
            {
              k = map_tensor4[j1][0];
              l = map_tensor4[j1][1];

              dik = KronDelta(i,k);
              dil = KronDelta(i,l);
              djk = KronDelta(j,k);
              djl = KronDelta(j,l);
              dkl = KronDelta(k,l);

              // Devioteric term (without stress terms)
              Cmat(i1,j1) += fact1*(dil*djk+dik*djl);
              Cmat(i1,j1) += fact2*(b(i,j)*dkl + dij*b(k,l));
              Cmat(i1,j1) += fact3*(dij*dkl);

              // Devioteric terms, with stress terms
              //Cmat(i1,j1) += fact1*(dil*djk);
              //Cmat(i1,j1) += fact2*(b(i,j)*dkl + dij*b(k,l));
              //Cmat(i1,j1) += fact3*(dij*dkl);
              //Cmat(i1,j1) += fact4*(b(j,l)*dik);

              if(MIXED_ELEMENT)
              {
                // pressure term
                Cmat(i1,j1) += pres*(dij*dkl - djk*dil - djl*dik);

                //Cmat(i1,j1) += pres*(dij*dkl - djk*dil);
              }
              else
              {
                // volumetric term (without stress terms)
                Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk)+fact9*(djl*dik);

                // volumetric term (with stress terms)
                //Cmat(i1,j1) += fact7*(dij*dkl)+fact8*(dil*djk);
              }
            }
        }

    return err;
}




int Matl_SaintVenantKirchhoff::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    double  mu     = matData[0];      // shear modulus
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    int     Utype  = int(matData[2]); // volumetric energy function

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
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






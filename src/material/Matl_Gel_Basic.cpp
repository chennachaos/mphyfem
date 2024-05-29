#include "Matl_Gel_Basic.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "util.h"


using namespace std;


Matl_Gel_Basic::Matl_Gel_Basic()
{

}


Matl_Gel_Basic::~Matl_Gel_Basic()
{

}



int Matl_Gel_Basic::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  chempotGrad, double& pres, VectorXd&  stre, VectorXd&  chemFlux, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact6, fact8, fact9;

    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K     = matData[0]; // bulk modulus
    double  G     = matData[1]; // shear modulus
    double  D     = matData[4];                            // solvent diffusivity
    //double  D     = 1.0e-9;
    double  Chi   = matData[5]; // Flory-Huggins interaction parameter

    double  R    = 8.3145;      // Gas constant (kg.m2/mol.s2.K)
    double  T     = 298.0;      // Reference temperature (K)


    VectorXd  b(9);
    Matrix3d  bb, FF, Finv;

    FF(0,0) = F(0);  FF(0,1) = F(3);  FF(0,2) = F(6);
    FF(1,0) = F(1);  FF(1,1) = F(4);  FF(1,2) = F(7);
    FF(2,0) = F(2);  FF(2,1) = F(5);  FF(2,2) = F(8);

    Finv = FF.inverse();

    computeLeftCGtensor(F, b);
    double  Ib = b(0) + b(4) + b(8);

    bb(0,0) = b(0);  bb(0,1) = b(3);  bb(0,2) = b(6);
    bb(1,0) = b(1);  bb(1,1) = b(4);  bb(1,2) = b(7);
    bb(2,0) = b(2);  bb(2,1) = b(5);  bb(2,2) = b(8);


    double  J     = determinant3b3(F);

    double  phi0 = 0.999;
    double  mu0  = -14382.975;                              // (Joules/mole)
    double  mu   = 0.0;
    double  Omega = 1.0e-4;
    double  Js = 1.0/phi0;
    double  Je = J/Js;
    double  CR = (Js-1.0)/Omega;
    double  md = D*CR/R/T/J;


    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    chemFlux.setZero();


        //Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + Js * U

        //Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + Js * K/2 * (J-1)^2
        //Psi = G/2 ( tr(C) - 3 - 2 ln(J) ) + Js * K/2 * (ln(J))^2

    //double  dUdJe = K*(Je-1.0), d2UJe2 = K;                 // U = K/2 * (J-1)^2
    double  dUdJe = K*log(Je)/Je, d2UJe2 = K*(1.0-log(Je))/Je/Je; // U = K/2 * (ln(J))^2


        /////////////////////////////////
        // Cauchy stress tensor

        // Mechanical part
        fact1 = G/J;
        fact2 = -fact1;

        if(MIXED_ELEMENT)
          fact2 += pres;
        else
          fact2 += (dUdJe - mu/Omega/Je);


        for(i=0; i<9; i++)
          stre[i] = fact1*b[i];

        stre[0] += fact2 ;
        stre[4] += fact2 ;
        stre[8] += fact2 ;

        // Mixing part

        /////////////////////////////////
        // Chemical flux vector

        chemFlux = - md*chempotGrad;

        /////////////////////////////////
        // material tangent tensor

        if(tangFlag)
        {
          // Mechanical part

          fact1 =  G/J;

          fact8 =  d2UJe2*Je + dUdJe - mu/Omega;
          fact9 =  -dUdJe + mu/Omega;

          for(i1=0; i1<9; i1++)
          {
            i = map_tensor4[i1][0];
            j = map_tensor4[i1][1];

            for(j1=0; j1<9; j1++)
            {
              k = map_tensor4[j1][0];
              l = map_tensor4[j1][1];

              Cmat(i1,j1) += fact1*(KronDelta(i,k)*bb(j,l)+KronDelta(i,l)*KronDelta(j,k));

              if(MIXED_ELEMENT)
              {
                // pressure term
                Cmat(i1,j1) += pres*(KronDelta(i,j)*KronDelta(k,l) - KronDelta(j,k)*KronDelta(i,l) - KronDelta(j,l)*KronDelta(i,k));
              }
              else
              {
                // volumetric term
                Cmat(i1,j1) += fact8*(KronDelta(i,j)*KronDelta(k,l));
                Cmat(i1,j1) += fact9*(KronDelta(i,l)*KronDelta(j,k));
              }
            }
          }

          /////////////////////////////////
          // chemical diffusivity tensor

          Amat.setZero();

          Amat(0,0) = -md;    Amat(1,1) = -md;    Amat(2,2) = -md;
        }


    return err;
}



int Matl_Gel_Basic::computeElectricComponents(VectorXd&  chempotGrad, VectorXd&  chemFlux, MatrixXd&  Amat)
{
    double  eps0  = matData[3]; // permittivity

    /////////////////////////////////
    // Electrical displacement vector

    chemFlux[0] = eps0*chempotGrad[0];
    chemFlux[1] = eps0*chempotGrad[1];
    chemFlux[2] = eps0*chempotGrad[2];


    /////////////////////////////////
    // electric permittivity tensor

    Amat.setZero();
    Amat(0,0) = eps0;  Amat(1,1) = eps0;  Amat(2,2) = eps0;

    return 0;
}





int Matl_Gel_Basic::computeElectricDisplacement(VectorXd&  chempotGrad, VectorXd&  chemFlux)
{
    double  eps0  = matData[3]; // permitivity

    chemFlux[0] = eps0*chempotGrad[0];
    chemFlux[1] = eps0*chempotGrad[1];
    chemFlux[2] = eps0*chempotGrad[2];

    return 0;
}



int Matl_Gel_Basic::computeMechanicalStress(VectorXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    double  K     = matData[0]; // bulk modulus
    double  G    = matData[1]; // shear modulus
    double  Im    = matData[2];
    double  eps0  = matData[3]; // permitivity
    double  Lamb  = K-2.0*G/3.0;

    VectorXd  b(9);

    computeLeftCGtensor(F, b);
    double  Ib = b(0) + b(4) + b(8);

    double  J     = determinant3b3(F);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  Ibbar = Jm2d3*Ib;
    double  gfact = 1.0/(1.0-(Ibbar-3.0)/Im);


    double  fact1 = G*Jm5d3*gfact;
    double  fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += K*(J-1.0);

    stre.setZero();
    for(int i=0; i<9; i++)
      stre[i] = fact1*b[i];

    stre[0] += fact2 ;
    stre[4] += fact2 ;
    stre[8] += fact2 ;

    return 0;
}



int Matl_Gel_Basic::computeMaxwellStress(VectorXd&  chempotGrad, VectorXd&  stre)
{
    double  eps0   = matData[3]; // permittivity

    stre.setZero();

    addElectricPartToCauchyStress(eps0, chempotGrad, stre);

    return 0;
}



int Matl_Gel_Basic::computeTotalStress(VectorXd&  F, VectorXd&  chempotGrad, double& pres, VectorXd&  stre)
{
    double  eps0   = matData[3]; // permittivity

    // Mechanical part
    computeMechanicalStress(F, pres, stre);

    // Electrical part
    addElectricPartToCauchyStress(eps0, chempotGrad, stre);

    return 0;
}



#include "Matl_ElecMech_Gent.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "util.h"


using namespace std;


Matl_ElecMech_Gent::Matl_ElecMech_Gent()
{

}


Matl_ElecMech_Gent::~Matl_ElecMech_Gent()
{

}



int Matl_ElecMech_Gent::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact6, fact8, fact9;

    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K     = matData[0]; // bulk modulus
    double  mu    = matData[1]; // shear modulus
    double  Im    = matData[2];
    double  eps0  = matData[3]; // permitivity
    double  Lamb  = K-2.0*mu/3.0;

    VectorXd  b(9);
    Matrix3d  bb,  bbdev;

    computeLeftCGtensor(F, b);
    double  Ib = b(0) + b(4) + b(8);

    bb(0,0) = b(0);  bb(0,1) = b(3);  bb(0,2) = b(6);
    bb(1,0) = b(1);  bb(1,1) = b(4);  bb(1,2) = b(7);
    bb(2,0) = b(2);  bb(2,1) = b(5);  bb(2,2) = b(8);

    fact = Ib*r1d3;

    bbdev = bb;
    bbdev(0,0) -= fact;
    bbdev(1,1) -= fact;
    bbdev(2,2) -= fact;

    double  J     = determinant3b3(F);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  Ibbar = Jm2d3*Ib;
    double  gfact = 1.0/(1.0-(Ibbar-3.0)/Im);

    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    elecDisp.setZero();


        //Psi = -mu/2 Im ln(1- (J^{-2/3} trb - 3)/Im) + K/2 (J-1)^2


        /////////////////////////////////
        // Cauchy stress tensor

        // Mechanical part
        fact1 = mu*Jm5d3*gfact;
        fact2 = -fact1*r1d3*Ib;

        if(MIXED_ELEMENT)
          fact2 += pres;
        else
          //fact2 += K*(J-1.0);
          fact2 += 0.1*K*(J*J*J*J-1.0/(J*J*J*J*J*J));


        for(i=0; i<9; i++)
          stre[i] = fact1*b[i];

        stre[0] += fact2 ;
        stre[4] += fact2 ;
        stre[8] += fact2 ;

        // Electrical part
        addElectricPartToCauchyStress(eps0, elecField, stre);

        /////////////////////////////////
        // Electrical displacement vector

        elecDisp[0] = eps0*elecField[0];
        elecDisp[1] = eps0*elecField[1];
        elecDisp[2] = eps0*elecField[2];

        /////////////////////////////////
        // material tangent tensor

        if(tangFlag)
        {
          // Mechanical part

          fact1 = -r2d3*mu*Jm5d3*gfact;
          fact2 =   2.0*mu*Jm7d3*gfact*gfact/Im;
          fact3 =  r1d3*mu*Jm5d3*gfact*Ib;
          fact4 = -r2d3*mu*Jm5d3*gfact;

          //fact8 = K*(2.0*J-1.0);
          //fact9 = -K*(J-1.0);

          fact8 =  0.5*K*(J*J*J*J + 1.0/(J*J*J*J*J*J));
          fact9 = -0.1*K*(J*J*J*J - 1.0/(J*J*J*J*J*J));

          for(i1=0; i1<9; i1++)
          {
            i = map_tensor4[i1][0];
            j = map_tensor4[i1][1];

            for(j1=0; j1<9; j1++)
            {
              k = map_tensor4[j1][0];
              l = map_tensor4[j1][1];

              // Devioteric terms, without stress terms
              Cmat(i1,j1) += fact1*(bbdev(i,j)*KronDelta(k,l));
              Cmat(i1,j1) += fact2*(bbdev(i,j)*bbdev(k,l));
              Cmat(i1,j1) += fact3*(KronDelta(i,k)*KronDelta(j,l)+KronDelta(i,l)*KronDelta(j,k));
              Cmat(i1,j1) += fact4*(bb(k,l)*KronDelta(i,j));

              if(MIXED_ELEMENT)
              {
                // pressure term
                Cmat(i1,j1) += pres*(KronDelta(i,j)*KronDelta(k,l) - KronDelta(j,k)*KronDelta(i,l) - KronDelta(j,l)*KronDelta(i,k));
              }
              else
              {
                // volumetric term
                Cmat(i1,j1) += fact8*(KronDelta(i,j)*KronDelta(k,l));
                Cmat(i1,j1) += fact9*(KronDelta(i,k)*KronDelta(j,l)+KronDelta(i,l)*KronDelta(j,k));
              }
            }
          }

          // Electrical part
          addElectricPartToMaterialTensor(eps0, elecField, Cmat);


          /////////////////////////////////
          // coupling tensor

          addToCouplingTensor(eps0, elecField, Bmat);

          /////////////////////////////////
          // electric permittivity tensor

          Amat(0,0) = eps0;  Amat(1,1) = eps0;  Amat(2,2) = eps0;
        }


    return err;
}



int Matl_ElecMech_Gent::computeElectricComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat)
{
    double  eps0  = matData[3]; // permittivity

    /////////////////////////////////
    // Electrical displacement vector

    elecDisp[0] = eps0*elecField[0];
    elecDisp[1] = eps0*elecField[1];
    elecDisp[2] = eps0*elecField[2];


    /////////////////////////////////
    // electric permittivity tensor

    Amat.setZero();
    Amat(0,0) = eps0;  Amat(1,1) = eps0;  Amat(2,2) = eps0;

    return 0;
}





int Matl_ElecMech_Gent::computeElectricDisplacement(VectorXd&  elecField, VectorXd&  elecDisp)
{
    double  eps0  = matData[3]; // permitivity

    elecDisp[0] = eps0*elecField[0];
    elecDisp[1] = eps0*elecField[1];
    elecDisp[2] = eps0*elecField[2];

    return 0;
}



int Matl_ElecMech_Gent::computeMechanicalStress(VectorXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    double  K     = matData[0]; // bulk modulus
    double  mu    = matData[1]; // shear modulus
    double  Im    = matData[2];
    double  eps0  = matData[3]; // permitivity
    double  Lamb  = K-2.0*mu/3.0;

    VectorXd  b(9);

    computeLeftCGtensor(F, b);
    double  Ib = b(0) + b(4) + b(8);

    double  J     = determinant3b3(F);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  Ibbar = Jm2d3*Ib;
    double  gfact = 1.0/(1.0-(Ibbar-3.0)/Im);


    double  fact1 = mu*Jm5d3*gfact;
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



int Matl_ElecMech_Gent::computeMaxwellStress(VectorXd&  elecField, VectorXd&  stre)
{
    double  eps0   = matData[3]; // permittivity

    stre.setZero();

    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}



int Matl_ElecMech_Gent::computeTotalStress(VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre)
{
    double  eps0   = matData[3]; // permittivity

    // Mechanical part
    computeMechanicalStress(F, pres, stre);

    // Electrical part
    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}



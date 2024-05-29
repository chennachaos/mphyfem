#include "Matl_ElecMech_ArrudaBoyce.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "util.h"


using namespace std;


Matl_ElecMech_ArrudaBoyce::Matl_ElecMech_ArrudaBoyce()
{

}


Matl_ElecMech_ArrudaBoyce::~Matl_ElecMech_ArrudaBoyce()
{

}



int Matl_ElecMech_ArrudaBoyce::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp)
{
    //Psi = c1(I - 3) + c2(I^2 - 3^2) + c3(I^3 - 3^3) + c4(I^4 - 3^4) + c5(I^5 - 3^5) + K/2 (J-1)^2

    double  fact, fact1, fact2, fact3, fact4, fact5, fact8, fact9;

    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = matData[0]; // bulk modulus
    double  mu     = matData[1]; // shear modulus
    double  N      = matData[2]; // number of links per polymer chain 
    double  eps0   = matData[3]; // permittivity
    double  c1     = mu/2.0;
    double  c2     = mu/20.0/N;
    double  c3     = mu*11.0/1050.0/N/N;
    double  c4     = mu*19.0/7000.0/N/N/N;
    double  c5     = mu*519.0/673750.0/N/N/N/N;

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
    double  Jm5d3 = pow(J, -r5d3);
    double  Jm7d3 = pow(J, -r7d3);
    double  Ibbar = Jm2d3*Ib;

    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    elecDisp.setZero();

    double a1 = c1 + 2.0*c2*Ibbar + 3.0*c3*Ibbar*Ibbar + 4.0*c4*Ibbar*Ibbar*Ibbar + 5.0*c5*Ibbar*Ibbar*Ibbar*Ibbar;
    double a2 = 2.0*c2 + 6.0*c3*Ibbar + 12.0*c4*Ibbar*Ibbar + 20.0*c5*Ibbar*Ibbar*Ibbar;

    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    fact1 = 2.0*a1*Jm5d3;
    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += K*(J-1.0);

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

        fact1 =  4.0*Jm7d3*a2;
        fact2 = -4.0*r1d3*a1*Jm5d3;
        fact3 =      r2d3*a1*Jm5d3*Ib;
        fact4 =  4.0*r1d3*r1d3*a1*Jm5d3*Ib;

        fact8 = -K*(J-1.0);
        fact9 =  K*(2.0*J-1.0);

        for(i1=0; i1<9; i1++)
        {
            i = map_tensor4[i1][0];
            j = map_tensor4[i1][1];

            for(j1=0; j1<9; j1++)
            {
              k = map_tensor4[j1][0];
              l = map_tensor4[j1][1];

              // Devioteric terms, without stress terms
              Cmat(i1,j1) += fact1*(bbdev(i,j)*bbdev(k,l));
              Cmat(i1,j1) += fact2*(bb(i,j)*KronDelta(k,l)+bb(k,l)*KronDelta(i,j));
              Cmat(i1,j1) += fact3*(KronDelta(i,k)*KronDelta(j,l)+KronDelta(i,l)*KronDelta(j,k));
              Cmat(i1,j1) += fact4*(KronDelta(i,j)*KronDelta(k,l));

              if(MIXED_ELEMENT)
              {
                // pressure term
                Cmat(i1,j1) += pres*(KronDelta(i,j)*KronDelta(k,l) - KronDelta(j,k)*KronDelta(i,l) - KronDelta(j,l)*KronDelta(i,k));
              }
              else
              {
                // volumetric term
                Cmat(i1,j1) += fact8*(KronDelta(i,k)*KronDelta(j,l)+KronDelta(i,l)*KronDelta(j,k));
                Cmat(i1,j1) += fact9*(KronDelta(i,j)*KronDelta(k,l));
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




int Matl_ElecMech_ArrudaBoyce::computeElectricComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat)
{
    double  eps0  = matData[3]; // permitivity

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



int Matl_ElecMech_ArrudaBoyce::computeElectricDisplacement(VectorXd&  elecField, VectorXd&  elecDisp)
{
    double  eps0  = matData[3]; // permitivity

    elecDisp[0] = eps0*elecField[0];
    elecDisp[1] = eps0*elecField[1];
    elecDisp[2] = eps0*elecField[2];

    return 0;
}



int Matl_ElecMech_ArrudaBoyce::computeMechanicalStress(VectorXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    int  i1, j1, i, j, k, l, err=0;

    double  K      = matData[0]; // bulk modulus
    double  mu     = matData[1]; // shear modulus
    double  N      = matData[2]; // number of links per polymer chain 
    double  eps0   = matData[3]; // permittivity
    double  c1     = mu/2.0;
    double  c2     = mu/20.0/N;
    double  c3     = mu*11.0/1050.0/N/N;
    double  c4     = mu*19.0/7000.0/N/N/N;
    double  c5     = mu*519.0/673750.0/N/N/N/N;

    VectorXd  b(9);
    Matrix3d  bb,  bbdev;

    computeLeftCGtensor(F, b);
    double  Ib = b(0) + b(4) + b(8);

    double  J     = determinant3b3(F);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = pow(J, -r5d3);
    double  Jm7d3 = pow(J, -r7d3);
    double  Ibbar = Jm2d3*Ib;

    double a1 = c1 + 2.0*c2*Ibbar + 3.0*c3*Ibbar*Ibbar + 4.0*c4*Ibbar*Ibbar*Ibbar + 5.0*c5*Ibbar*Ibbar*Ibbar*Ibbar;
    double a2 = 2.0*c2 + 6.0*c3*Ibbar + 12.0*c4*Ibbar*Ibbar + 20.0*c5*Ibbar*Ibbar*Ibbar;

    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    double  fact1 = 2.0*a1*Jm5d3;
    double  fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += K*(J-1.0);

    stre.setZero();

    for(i=0; i<9; i++)
      stre[i] = fact1*b[i];

    stre[0] += fact2 ;
    stre[4] += fact2 ;
    stre[8] += fact2 ;

    return 0;
}



int Matl_ElecMech_ArrudaBoyce::computeMaxwellStress(VectorXd&  elecField, VectorXd&  stre)
{
    double  eps0   = matData[3]; // permittivity

    stre.setZero();

    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}



int Matl_ElecMech_ArrudaBoyce::computeTotalStress(VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre)
{
    double  eps0   = matData[3]; // permittivity

    // Mechanical part
    computeMechanicalStress(F, pres, stre);

    // Electrical part
    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}




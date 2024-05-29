#include "Matl_ElecMech_MooneyRivlin.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "util.h"


using namespace std;


Matl_ElecMech_MooneyRivlin::Matl_ElecMech_MooneyRivlin()
{

}


Matl_ElecMech_MooneyRivlin::~Matl_ElecMech_MooneyRivlin()
{

}



int Matl_ElecMech_MooneyRivlin::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact8, fact9;

    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = matData[0]; // bulk modulus
    double  mu     = matData[1]; // shear modulus
    double  eps0   = matData[2]; // permittivity

    VectorXd  b(9);
    Matrix3d  bb;

    computeLeftCGtensor(F, b);

    bb(0,0) = b(0);  bb(0,1) = b(3);  bb(0,2) = b(6);
    bb(1,0) = b(1);  bb(1,1) = b(4);  bb(1,2) = b(7);
    bb(2,0) = b(2);  bb(2,1) = b(5);  bb(2,2) = b(8);

    double  J     = determinant3b3(F);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b(0) + b(4) + b(8);

    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    elecDisp.setZero();


        //Psi = mu/2 (J^{-2/3} trb - 3) + K/2 (J-1)^2


        /////////////////////////////////
        // Cauchy stress tensor

        // Mechanical part
        fact1 = mu*Jm5d3;

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

          fact1 = mu*Jm5d3*r1d3;

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

              // Devioteric term (without stress terms)
              Cmat(i1,j1) += fact1*Ib*(KronDelta(i,l)*KronDelta(j,k)+KronDelta(i,k)*KronDelta(j,l));
              Cmat(i1,j1) += fact1*(-2.0*bb(i,j)*KronDelta(k,l) - 2.0*KronDelta(i,j)*bb(k,l));
              Cmat(i1,j1) += fact1*r2d3*Ib*(KronDelta(i,j)*KronDelta(k,l));

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



/*
int Matl_ElecMech_MooneyRivlin::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp)
{
    double  fact, fact1, fact2, fact3, fact4, fact5;

    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    int  NHtype = max(1,(int) matData[0] );
    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K     = matData[1]; // bulk modulus
    double  mu    = matData[2]; // shear modulus
    double  Lamb  = K-2.0*mu/3.0;
    double  eps0  = matData[3]; // permitivity
    double  eFeF = elecField[0]*elecField[0]+elecField[1]*elecField[1]+elecField[2]*elecField[2];
    double  eFeFd2 = 0.5*eFeF;

    VectorXd  b(9);
    Matrix3d  bb;

    computeLeftCGtensor(F, b);

    bb(0,0) = b(0);  bb(0,1) = b(3);  bb(0,2) = b(6);
    bb(1,0) = b(1);  bb(1,1) = b(4);  bb(1,2) = b(7);
    bb(2,0) = b(2);  bb(2,1) = b(5);  bb(2,2) = b(8);

    double  J     = determinant3b3(F);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b(0) + b(4) + b(8);

    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    elecDisp.setZero();

    switch(NHtype)
    {
      case 11:

        //Psi = mu/2 (trb-3) - mu lnJ + Lamb/2 (lnJ)^2
        //Psi = mu/2 (trb-3) - mu lnJ + K/2 (J-1)^2


        /////////////////////////////////
        // Cauchy stress tensor

        // Mechanical part
        fact1 = mu/J;
        //fact2 = Lamb*log(J)/J - mu/J;
        fact2 = K*(J-1.0) - mu/J;

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
          fact1 = mu/J;
          fact2 = mu/J;  // - Lamb*log(J)/J;
          //fact3 = Lamb/J;

          fact3 = -K*(J-1.0);
          fact5 =  K*(2.0*J-1.0);

          for(i1=0; i1<9; i1++)
          {
            i = map_tensor4[i1][0];
            j = map_tensor4[i1][1];

            for(j1=0; j1<9; j1++)
            {
              k = map_tensor4[j1][0];
              l = map_tensor4[j1][1];

              //  stress terms included
              //Cmat(i1,j1) += fact1*bb(j,l)*KronDelta(i,k);
              //Cmat(i1,j1) += fact2*(KronDelta(j,k)*KronDelta(l,i));
              //Cmat(i1,j1) += fact3*(KronDelta(j,i)*KronDelta(l,k));

              //  stress terms excluded
              //Cmat(i1,j1) += fact3*(KronDelta(j,i)*KronDelta(l,k));
              Cmat(i1,j1) += fact2*(KronDelta(j,k)*KronDelta(l,i));
              Cmat(i1,j1) += fact2*(KronDelta(j,l)*KronDelta(i,k));

              Cmat(i1,j1) += fact3*(KronDelta(i,k)*KronDelta(j,l)+KronDelta(i,l)*KronDelta(j,k));
              Cmat(i1,j1) += fact5*(KronDelta(i,j)*KronDelta(k,l));
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

        break;

      case 12:

        //Psi = mu/2 (J^{-2/3} trb - 3) + Lamb/2 (lnJ)^2
        //Psi = mu/2 (J^{-2/3} trb - 3) + K/2 (J-1)^2


        /////////////////////////////////
        // Cauchy stress tensor

        // Mechanical part
        fact1 = mu*Jm5d3;
        //fact2 = Lamb*log(J)/J - mu*Jm5d3*r1d3*Ib;
        fact2 = -fact1*r1d3*Ib + pres;

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

          fact2 =  r2d3*r1d3*mu*Jm5d3*Ib;
          fact  =  r1d3*mu*Jm5d3*Ib;
          fact1 =  fact * 2.0;
          fact3 = -r2d3*mu*Jm5d3;

          //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n", fact, fact1, fact2, fact3);
          //cout <<  "pres = " <<  pres <<  endl;

          fact1 = mu*Jm5d3*r1d3;

          for(i1=0; i1<9; i1++)
          {
            i = map_tensor4[i1][0];
            j = map_tensor4[i1][1];

            for(j1=0; j1<9; j1++)
            {
              k = map_tensor4[j1][0];
              l = map_tensor4[j1][1];

              // Devioteric terms, with stress terms
              //Cmat(i1,j1) += fact1*Ib*(KronDelta(i,l)*KronDelta(j,k));
              //Cmat(i1,j1) += fact1*(3.0*bb(j,l)*KronDelta(i,k) - 2.0*bb(i,j)*KronDelta(k,l) - 2.0*KronDelta(i,j)*bb(k,l));
              //Cmat(i1,j1) += fact1*r2d3*Ib*(KronDelta(i,j)*KronDelta(k,l));

              // Devioteric terms, without stress terms
              Cmat(i1,j1) += fact1*Ib*(KronDelta(i,l)*KronDelta(j,k)+KronDelta(i,k)*KronDelta(j,l));
              Cmat(i1,j1) += fact1*(-2.0*bb(i,j)*KronDelta(k,l) - 2.0*KronDelta(i,j)*bb(k,l));
              Cmat(i1,j1) += fact1*r2d3*Ib*(KronDelta(i,j)*KronDelta(k,l));

              // pressure term
              Cmat(i1,j1) += pres*(KronDelta(i,j)*KronDelta(k,l) - KronDelta(j,k)*KronDelta(i,l) - KronDelta(j,l)*KronDelta(i,k));
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

        break;

      default:

          cerr << " This NeoHooke Material is not avaialble yet " << endl;
          err = -2;

        break;
    }

    return err;
}
*/


int Matl_ElecMech_MooneyRivlin::computeElectricComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat)
{
    double  eps0  = matData[2]; // permitivity

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



int Matl_ElecMech_MooneyRivlin::computeElectricDisplacement(VectorXd&  elecField, VectorXd&  elecDisp)
{
    double  eps0  = matData[2]; // permitivity

    elecDisp[0] = eps0*elecField[0];
    elecDisp[1] = eps0*elecField[1];
    elecDisp[2] = eps0*elecField[2];

    return 0;
}



int Matl_ElecMech_MooneyRivlin::computeMechanicalStress(VectorXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    double  K      = matData[0]; // bulk modulus
    double  mu     = matData[1]; // shear modulus

    VectorXd  b(9);

    computeLeftCGtensor(F, b);

    double  J     = determinant3b3(F);
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b(0) + b(4) + b(8);

    double  fact1 = mu*Jm5d3;
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



int Matl_ElecMech_MooneyRivlin::computeMaxwellStress(VectorXd&  elecField, VectorXd&  stre)
{
    double  eps0   = matData[2]; // permittivity

    stre.setZero();

    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}



int Matl_ElecMech_MooneyRivlin::computeTotalStress(VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre)
{
    double  eps0   = matData[2]; // permittivity

    // Mechanical part
    computeMechanicalStress(F, pres, stre);

    // Electrical part
    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}




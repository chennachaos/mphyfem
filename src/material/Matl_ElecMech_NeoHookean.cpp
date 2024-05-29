#include "Matl_ElecMech_NeoHookean.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"


using namespace std;


Matl_ElecMech_NeoHookean::Matl_ElecMech_NeoHookean()
{

}


Matl_ElecMech_NeoHookean::~Matl_ElecMech_NeoHookean()
{

}



int Matl_ElecMech_NeoHookean::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn,  MatrixXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, MatrixXd& intVarPrev, MatrixXd& intVar, int gp, double dt)
{
    //Psi = mu/2 (J^{-2/3} trb - 3) + K/2 (J-1)^2

    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;

    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  K      = matData[0]; // bulk modulus
    double  mu     = matData[1]; // shear modulus
    double  eps0   = matData[2]; // permittivity
    int     Utype  = int(matData[2]); // volumetric energy function
    double  Lambda = K - 2.0*mu/3.0;

    Matrix3d  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Ib = b.trace();

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);
    double  dij, dik, dil, djk, djl, dkl;

    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    elecDisp.setZero();



    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    fact1 = mu*Jm5d3;

    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 +=  volumetricfunctions_getFirstDerivative(8, K, J);


    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

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

        // without stress terms
        fact1 =  r1d3*mu*Jm5d3*Ib;
        fact2 = -r2d3*mu*Jm5d3;
        fact3 =  r2d3*fact1;

        //printf(" %14.12f \t %14.12f \t %14.12f \t %14.12f \n", fact, fact1, fact2, fact3);
        //cout <<  "pres = " <<  pres <<  endl;

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
int Matl_ElecMech_NeoHookean::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp)
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


int Matl_ElecMech_NeoHookean::computeElectricComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat)
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



int Matl_ElecMech_NeoHookean::computeElectricDisplacement(VectorXd&  elecField, VectorXd&  elecDisp)
{
    double  eps0  = matData[2]; // permitivity

    elecDisp[0] = eps0*elecField[0];
    elecDisp[1] = eps0*elecField[1];
    elecDisp[2] = eps0*elecField[2];

    return 0;
}



int Matl_ElecMech_NeoHookean::computeMechanicalStress(VectorXd&  F, double&  pres, VectorXd&  stre)
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



int Matl_ElecMech_NeoHookean::computeMaxwellStress(VectorXd&  elecField, VectorXd&  stre)
{
    double  eps0   = matData[2]; // permittivity

    stre.setZero();

    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}



int Matl_ElecMech_NeoHookean::computeTotalStress(VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre)
{
    double  eps0   = matData[2]; // permittivity

    // Mechanical part
    computeMechanicalStress(F, pres, stre);

    // Electrical part
    addElectricPartToCauchyStress(eps0, elecField, stre);

    return 0;
}




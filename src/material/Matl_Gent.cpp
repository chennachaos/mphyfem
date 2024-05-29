#include "Matl_Gent.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"


using namespace std;


Matl_Gent::Matl_Gent()
{

}


Matl_Gent::~Matl_Gent()
{

}



int Matl_Gent::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = -mu/2 Im ln(1- (J^{-2/3} trb - 3)/Im) + U(J)
    // U(J) is selected based on the input value for 'Utype'

    int err =  computeStressAndTangent_Gent(sss, MIXED_ELEMENT, Utype, Kinv, data_deviatoric, Fn, F, pres, stre, Cmat);

    return err;
}




int Matl_Gent::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r5d3 = 5.0/3.0, r7d3 = 7.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  mu     = matData[0];      // shear modulus
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    double  Im     = matData[2];
    int     Utype  = int(matData[3]); // volumetric energy function

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = Jm2d3/J;
    double  Jm7d3 = pow(J, -r7d3);
    double  Ib = b(0,0) + b(1,1) + b(2,2);
    double  Ibbar = Jm2d3*Ib;
    double  gfact = 1.0/(1.0-(Ibbar-3.0)/Im);

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);
    double  d2UdJ2 = volumetricfunctions_getSecondDerivative(Utype, K, J);

    stre.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    fact1 = mu*Jm5d3*gfact;
    fact2 = -fact1*r1d3*Ib;

    if(MIXED_ELEMENT)
      fact2 += pres;
    else
      fact2 += dUdJ;

    stre[0] = fact1*b(0,0);    stre[3] = fact1*b(0,1);    stre[6] = fact1*b(0,2);
    stre[1] = fact1*b(1,0);    stre[4] = fact1*b(1,1);    stre[7] = fact1*b(1,2);
    stre[2] = fact1*b(2,0);    stre[5] = fact1*b(2,1);    stre[8] = fact1*b(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;

    return 0;
}






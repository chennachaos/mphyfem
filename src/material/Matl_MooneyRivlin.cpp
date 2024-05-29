#include "Matl_MooneyRivlin.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"


using namespace std;


Matl_MooneyRivlin::Matl_MooneyRivlin()
{

}


Matl_MooneyRivlin::~Matl_MooneyRivlin()
{

}



int Matl_MooneyRivlin::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = \Psi_Dev + U(J)
    // U(J) is selected based on the input value for 'Utype'

    int err =  computeStressAndTangent_MooneyRivlin(sss, MIXED_ELEMENT, Utype, Kinv, data_deviatoric, Fn, F, pres, stre, Cmat);

    return err;
}




int Matl_MooneyRivlin::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
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






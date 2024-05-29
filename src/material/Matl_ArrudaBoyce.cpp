#include "Matl_ArrudaBoyce.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"


using namespace std;


Matl_ArrudaBoyce::Matl_ArrudaBoyce()
{

}


Matl_ArrudaBoyce::~Matl_ArrudaBoyce()
{

}



int Matl_ArrudaBoyce::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = c1(I - 3) + c2(I^2 - 3^2) + c3(I^3 - 3^3) + c4(I^4 - 3^4) + c5(I^5 - 3^5) + U(J()

    int err =  computeStressAndTangent_ArrudaBoyce(sss, MIXED_ELEMENT, Utype, Kinv, matData, Fn, F, pres, stre, Cmat);

    return err;
}





int Matl_ArrudaBoyce::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    int  i1, j1, i, j, k, l, err=0;

    int     Utype  = int(matData[3]);                       // volumetric energy function

    double  mu     = matData[0];      // shear modulus
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus
    double  N      = matData[2]; // number of links per polymer chain 
    double  c1     = mu/2.0;
    double  c2     = mu/20.0/N;
    double  c3     = mu*11.0/1050.0/N/N;
    double  c4     = mu*19.0/7000.0/N/N/N;
    double  c5     = mu*519.0/673750.0/N/N/N/N;

    // Left Cauchy-Green tensor;
    MatrixXd  b = F*F.transpose();

    double  J     = F.determinant();
    double  Jm2d3 = pow(J, -r2d3);
    double  Jm5d3 = pow(J, -r5d3);
    double  Jm7d3 = pow(J, -r7d3);
    double  Ib    = b(0,0) + b(1,1) + b(2,2);
    double  Ibbar = Jm2d3*Ib;

    double a1 = c1 + 2.0*c2*Ibbar + 3.0*c3*Ibbar*Ibbar + 4.0*c4*Ibbar*Ibbar*Ibbar + 5.0*c5*Ibbar*Ibbar*Ibbar*Ibbar;
    double a2 = 2.0*c2 + 6.0*c3*Ibbar + 12.0*c4*Ibbar*Ibbar + 20.0*c5*Ibbar*Ibbar*Ibbar;

    double  dUdJ   = volumetricfunctions_getFirstDerivative(Utype, K, J);

    /////////////////////////////////
    // Cauchy stress tensor

    // Mechanical part
    double  fact1 = 2.0*a1*Jm5d3;
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






#include "Matl_Longevin8chain_Viscoelastic.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"


using namespace std;


Matl_Longevin8chain_Viscoelastic::Matl_Longevin8chain_Viscoelastic()
{

}


Matl_Longevin8chain_Viscoelastic::~Matl_Longevin8chain_Viscoelastic()
{

}





int Matl_Longevin8chain_Viscoelastic::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = mu/2 (J^{-2/3} trb - 3) + U(J)
    // U(J) is selected based on the input value for 'Utype'

    // Compute the stress and tangent tensors for the elastic part
    int err =  computeStressAndTangent_Longevin8chain(sss, MIXED_ELEMENT, Utype, Kinv, data_deviatoric, Fn, F, pres, stre, Cmat);


    //printVector(stre);
    // Compute the stress and tangent tensors for the viscoelastic part
    SetTimeParametersFluid(tis, spectralRadius, dt, td);

    err = computeStressAndTangent_Viscoelasticity_Model1(data_viscoelastic, gp, F, td, ivar, dt, stre, Cmat);
    //printVector(stre);


    return err;
}




int Matl_Longevin8chain_Viscoelastic::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
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






#include "Matl_MagnMech_Gent_Viscoelastic.h"
#include "headersEigen.h"
#include "utilitiesmaterial.h"
#include "utilitiesmaterialHyperelastic.h"
#include "util.h"
#include "TimeFunction.h"
#include "MyTime.h"

extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime                myTime;
extern bool debug;

using namespace std;


Matl_MagnMech_Gent_Viscoelastic::Matl_MagnMech_Gent_Viscoelastic()
{

}


Matl_MagnMech_Gent_Viscoelastic::~Matl_MagnMech_Gent_Viscoelastic()
{

}



int Matl_MagnMech_Gent_Viscoelastic::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = mu/2 (J^{-2/3} trb - 3) + U(J)
    // U(J) is selected based on the input value for 'Utype'

    int err =  computeStressAndTangent_Gent(sss, MIXED_ELEMENT, Utype, Kinv, data_deviatoric, Fn, F, pres, stre, Cmat);

    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l;

    double  mu0   = data_magnfield[0];      // permeability
    double  J     = F.determinant();

    /////////////////////////////////
    // Cauchy stress tensor

    MatrixXd  streMagn(3,3);

    //printVector(ApplMagnfield);
    //printVector(ResiMagnfield);

    //printMatrix(streMagn);

    VectorXd  ResiMagnfieldCur = F*(ResiMagnfield/J);
    streMagn = (-1.0/mu0) * ((timeFunctions[0]->getValue()*ApplMagnfield)* (ResiMagnfieldCur.transpose()) );


    stre(0) += streMagn(0,0);    stre(3) += streMagn(0,1);    stre(6) += streMagn(0,2);
    stre(1) += streMagn(1,0);    stre(4) += streMagn(1,1);    stre(7) += streMagn(1,2);
    stre(2) += streMagn(2,0);    stre(5) += streMagn(2,1);    stre(8) += streMagn(2,2);

    /////////////////////////////////
    // material tangent tensor
    // this part due to the residual-applied field is zero


    /////////////////////////////////
    // Viscoelastic part

    SetTimeParametersFluid(tis, spectralRadius, dt, td);

    err = computeStressAndTangent_Viscoelasticity_Model1(data_viscoelastic, gp, F, td, ivar, dt, stre, Cmat);

    return err;
}






int Matl_MagnMech_Gent_Viscoelastic::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre, VectorXd&  magnDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt)
{
    //Psi = -mu/2 Im ln(1- (J^{-2/3} trb - 3)/Im) + K/2 (J-1)^2

    int err =  computeStressAndTangent_Gent(sss, MIXED_ELEMENT, Utype, Kinv, matData, Fn, F, pres, stre, Cmat);


    double  mu0   = data_magnfield[0];      // permeability
    double  J     = F.determinant();

    Bmat.setZero();
    Amat.setZero();
    magnDisp.setZero();

    // Magnetic part
    addElectricPartToCauchyStress(mu0, magnField, stre);

    /////////////////////////////////
    // Magnetic displacement vector

    magnDisp[0] = mu0*magnField[0];
    magnDisp[1] = mu0*magnField[1];
    magnDisp[2] = mu0*magnField[2];

    /////////////////////////////////
    // material tangent tensor

    // Magnetic part
    addElectricPartToMaterialTensor(mu0, magnField, Cmat);


    /////////////////////////////////
    // coupling tensor

    addToCouplingTensor(mu0, magnField, Bmat);

    /////////////////////////////////
    // magnetic permeability tensor

    Amat(0,0) = mu0;  Amat(1,1) = mu0;  Amat(2,2) = mu0;

    return err;
}



int Matl_MagnMech_Gent_Viscoelastic::computeMagneticComponents(VectorXd&  magnField, VectorXd&  magnDisp, MatrixXd&  Amat)
{
    double  mu0   = data_magnfield[0];      // permeability

    /////////////////////////////////
    // Electrical displacement vector

    magnDisp[0] = mu0*magnField[0];
    magnDisp[1] = mu0*magnField[1];
    magnDisp[2] = mu0*magnField[2];


    /////////////////////////////////
    // electric permittivity tensor

    Amat.setZero();
    Amat(0,0) = mu0;  Amat(1,1) = mu0;  Amat(2,2) = mu0;

    return 0;
}





int Matl_MagnMech_Gent_Viscoelastic::computeMagneticDisplacement(VectorXd&  magnField, VectorXd&  magnDisp)
{
    double  mu0   = data_magnfield[0];      // permeability

    magnDisp[0] = mu0*magnField[0];
    magnDisp[1] = mu0*magnField[1];
    magnDisp[2] = mu0*magnField[2];

    return 0;
}




int Matl_MagnMech_Gent_Viscoelastic::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    double  r1d3 = 1.0/3.0;
    double  r2d3 = 2.0/3.0;
    double  r5d3 = 5.0/3.0;
    double  r7d3 = 7.0/3.0;

    double  K     = 1.0/Kinv;                               // bulk modulus
    double  mu    = matData[0];                             // shear modulus
    double  Im    = matData[1];
    double  Lamb  = K-2.0*mu/3.0;

    double  mu0   = data_magnfield[0];                      // permeability

    VectorXd  b(9), FF(9);
    FF(0) = F(0,0);    FF(3) = F(0,1);    FF(6) = F(0,2);
    FF(1) = F(1,0);    FF(4) = F(1,1);    FF(7) = F(1,2);
    FF(2) = F(2,0);    FF(5) = F(2,1);    FF(8) = F(2,2);
    computeLeftCGtensor(FF, b);

    double  Ib = b(0) + b(4) + b(8);

    double  J     = F.determinant();
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



int Matl_MagnMech_Gent_Viscoelastic::computeMaxwellStress(VectorXd&  magnField, VectorXd&  stre)
{
    double  mu0   = matData[3]; // permittivity

    stre.setZero();

    addElectricPartToCauchyStress(mu0, magnField, stre);

    return 0;
}



int Matl_MagnMech_Gent_Viscoelastic::computeTotalStress(MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre)
{
    double  mu0   = matData[3]; // permittivity

    // Mechanical part
    computeMechanicalStress(F, pres, stre);

    // Electrical part
    addElectricPartToCauchyStress(mu0, magnField, stre);

    return 0;
}



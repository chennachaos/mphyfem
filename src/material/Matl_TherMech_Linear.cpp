#include "Matl_TherMech_Linear.h"
#include "headersEigen.h"
#include "util.h"
#include "utilitiesmaterial.h"


using namespace std;


Matl_TherMech_Linear::Matl_TherMech_Linear()
{

}


Matl_TherMech_Linear::~Matl_TherMech_Linear()
{

}



int Matl_TherMech_Linear::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  tempGrad, double& temp, double& pres, VectorXd&  stre, VectorXd&  heatflux, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Kmat, double* iv1, double* iv2, double dt, int& nivGP, int gp)
{
    /*
     xx -> 11   yx -> 21  zx -> 31
     xy -> 12   yy -> 22  zy -> 32
     xz -> 13   yz -> 23  zz -> 33
    */

    int isw = 1,  err = 0;

    // compute stress and tangent matrix

    int ii, jj;
    double  fact1, fact2;

    double  T0       = matData[0];                           // T0
    double  Tinf     = matData[1];                           // Tinf
    double  rho0     = matData[2];                           // density
    double  BULK     = matData[3];                           // bulk modulus
    double  mu       = matData[4];                           // shear modulus
    double  alpha    = matData[5];                           // thermal expansion coefficient
    double  cV       = matData[6];                           // specific heat per unit volume at constant strain
    double  condCoef = matData[7];                           // thermal conductivity
    double  convCoef = matData[8];                           // convection coefficient

    double  kappa = alpha*3.0*BULK;

    double  Lamb  = BULK-2.0*mu/3.0;

    VectorXd  eps(9);

    Cmat.setZero();
    Bmat.setZero();
    Kmat.setZero();
    stre.setZero();
    heatflux.setZero();


        // small-strain tensor
        defgrad2eps(F, eps);

        // ()*dij*dkl
        Cmat(0,0) += Lamb;  Cmat(0,4) += Lamb;  Cmat(0,8) += Lamb;
        Cmat(4,0) += Lamb;  Cmat(4,4) += Lamb;  Cmat(4,8) += Lamb;
        Cmat(8,0) += Lamb;  Cmat(8,4) += Lamb;  Cmat(8,8) += Lamb;

        //printMatrix(Cmat);

        // ()*dil*djk
        Cmat(0,0) += mu;  Cmat(1,3) += mu;  Cmat(2,6) += mu;
        Cmat(3,1) += mu;  Cmat(4,4) += mu;  Cmat(5,7) += mu;
        Cmat(6,2) += mu;  Cmat(7,5) += mu;  Cmat(8,8) += mu;

        //printMatrix(Cmat);

        // ()*dik*dlj
        for(ii=0; ii<9; ii++)
          Cmat(ii,ii) += mu;



        Kmat(0,0) =  condCoef;    Kmat(1,1) =  condCoef;    Kmat(2,2) =  condCoef;

        double  fact = -kappa*(temp-T0);

        stre     = Cmat*eps;
        stre(0) += fact;
        stre(4) += fact;
        stre(8) += fact;

        heatflux = -Kmat*tempGrad;

    //printVector(eps);
    //printVector(stre);
    //printVector(elecDisp);

    return err;
}

#include "Matl_LinearElastic.h"
#include "headersEigen.h"
#include "util.h"
#include "utilitiesmaterial.h"


using namespace std;


Matl_LinearElastic::Matl_LinearElastic()
{

}


Matl_LinearElastic::~Matl_LinearElastic()
{

}



int Matl_LinearElastic::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    /*
     xx -> 11   yx -> 21  zx -> 31
     xy -> 12   yy -> 22  zy -> 32
     xz -> 13   yz -> 23  zz -> 33
    */

    int err = 0;

    double  YoungsModulus = matData[0], nu = matData[1];

    double  K     = YoungsModulus/(3.0*(1.0-2.0*nu)); // bulk modulus
    double  mu    = YoungsModulus/(2.0*(1+nu)); // shear modulus

    Cmat.setZero();
    stre.setZero();


    MatrixXd  eps(3,3), epsDev(3,3);

    get_smallstraintensor(F, eps);
    get_deviatoricpart(eps, epsDev);

    double  fact1 = 2.0*mu, fact2;

    if(MIXED_ELEMENT)
      fact2 = pres;
    else
      fact2 = K*eps.trace();

    stre[0] = fact1*epsDev(0,0);    stre[3] = fact1*epsDev(0,1);    stre[6] = fact1*epsDev(0,2);
    stre[1] = fact1*epsDev(1,0);    stre[4] = fact1*epsDev(1,1);    stre[7] = fact1*epsDev(1,2);
    stre[2] = fact1*epsDev(2,0);    stre[5] = fact1*epsDev(2,1);    stre[8] = fact1*epsDev(2,2);
    stre[0] += fact2 ;              stre[4] += fact2 ;              stre[8] += fact2 ;


    if(MIXED_ELEMENT)
        fact2 = -2.0*mu/3.0;
    else
        fact2 = K-2.0*mu/3.0;

    // ()*dij*dkl
    Cmat(0,0) += fact2;  Cmat(0,4) += fact2;  Cmat(0,8) += fact2;
    Cmat(4,0) += fact2;  Cmat(4,4) += fact2;  Cmat(4,8) += fact2;
    Cmat(8,0) += fact2;  Cmat(8,4) += fact2;  Cmat(8,8) += fact2;

    // ()*dil*djk
    Cmat(0,0) += mu;  Cmat(1,3) += mu;  Cmat(2,6) += mu;
    Cmat(3,1) += mu;  Cmat(4,4) += mu;  Cmat(5,7) += mu;
    Cmat(6,2) += mu;  Cmat(7,5) += mu;  Cmat(8,8) += mu;

    // ()*dik*dlj
    for(int ii=0; ii<9; ii++)
        Cmat(ii,ii) += mu;


    return err;
}




int  Matl_LinearElastic::computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
{
    int err = 0;

    double  YoungsModulus = matData[0], nu = matData[1];

    double  K     = YoungsModulus/(3.0*(1.0-2.0*nu)); // bulk modulus
    double  mu    = YoungsModulus/(2.0*(1+nu)); // shear modulus

    stre.setZero();


    MatrixXd  eps(3,3), epsDev(3,3);

    get_smallstraintensor(F, eps);
    get_deviatoricpart(eps, epsDev);

    double  fact1 = 2.0*mu, fact2;

    if(MIXED_ELEMENT)
      fact2 = pres;
    else
      fact2 = K*eps.trace();

    stre[0] = fact1*epsDev(0,0);    stre[3] = fact1*epsDev(0,1);    stre[6] = fact1*epsDev(0,2);
    stre[1] = fact1*epsDev(1,0);    stre[4] = fact1*epsDev(1,1);    stre[7] = fact1*epsDev(1,2);
    stre[2] = fact1*epsDev(2,0);    stre[5] = fact1*epsDev(2,1);    stre[8] = fact1*epsDev(2,2);
    stre[0] += fact2 ;              stre[4] += fact2 ;              stre[8] += fact2 ;

    return err;
}

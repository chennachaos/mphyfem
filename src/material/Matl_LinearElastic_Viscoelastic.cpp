#include "Matl_LinearElastic_Viscoelastic.h"
#include "headersEigen.h"
#include "util.h"
#include "utilitiesmaterial.h"


using namespace std;


Matl_LinearElastic_Viscoelastic::Matl_LinearElastic_Viscoelastic()
{

}


Matl_LinearElastic_Viscoelastic::~Matl_LinearElastic_Viscoelastic()
{

}



int Matl_LinearElastic_Viscoelastic::computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
{
    /*
     xx -> 11   yx -> 21  zx -> 31
     xy -> 12   yy -> 22  zy -> 32
     xz -> 13   yz -> 23  zz -> 33
    */


    int isw = 1,  err = 0;

    // check for the validity of input data and assumptions

    int ii, jj;
    double  K      = 1.0/matData[1];                        // bulk modulus
    double  muInf  = matData[2]; // shear modulus

    double  mu0    = matData[3]; // mu0

    int  numSeries = (int) matData[4];  // number of terms in the Prony series

    double  muVec[numSeries], lambdaVec[numSeries];

    for(ii=0; ii<numSeries; ii++)
    {
        jj = 2*ii;
        muVec[ii]     = matData[5+jj];
        lambdaVec[ii] = matData[5+jj+1];
    }

    double  mu1 = matData[5];
    double  lambda1 = matData[6];


    Cmat.setZero();
    stre.setZero();


    MatrixXd  eps(3,3), epsn(3,3), epsDev(3,3), epsnDev(3,3);

    get_smallstraintensor(F, eps);
    get_deviatoricpart(eps, epsDev);

    get_smallstraintensor(Fn, epsn);
    get_deviatoricpart(epsn, epsnDev);


    VectorXd  qVecPrev =  ivar.varPrev.block<9,1>(0,gp), qVec(9);
    MatrixXd  qVecPrevMat(3,3);
    vector2matrix(qVecPrev, qVecPrevMat);

    MatrixXd  qVecMat = exp(-dt/lambda1)*qVecPrevMat + (lambda1/dt)*(1.0-exp(-dt/lambda1)) * (epsDev-epsnDev);
    matrix2vector(qVecMat, qVec);


    // update the internal variables vector
    ivar.var.block<9,1>(0,gp) = qVec ;


    MatrixXd  streDev = 2.0*muInf*(mu0*epsDev+mu1*qVecMat);

    double  fact2;
    if(MIXED_ELEMENT)
      fact2 = pres;
    else
      fact2 = K*eps.trace();

    stre[0] = streDev(0,0);    stre[3] = streDev(0,1);    stre[6] = streDev(0,2);
    stre[1] = streDev(1,0);    stre[4] = streDev(1,1);    stre[7] = streDev(1,2);
    stre[2] = streDev(2,0);    stre[5] = streDev(2,1);    stre[8] = streDev(2,2);
    stre[0] += fact2 ;         stre[4] += fact2 ;         stre[8] += fact2 ;


    double  mu = muInf*(mu0 + (mu1*lambda1/dt)*(1.0-exp(-dt/lambda1)) );

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
    for(ii=0; ii<9; ii++)
        Cmat(ii,ii) += mu;

    return err;
}


int  Matl_LinearElastic_Viscoelastic::computeMechanicalStress(MatrixXd&  F, double&  volstrainnew, VectorXd&  stre)
{
    double  fact, fact1, fact2, fact3, fact4, fact5, fact7, fact8, fact9;
    double  r1d3 = 1.0/3.0,  r2d3 = 2.0/3.0,  r5d3 = 5.0/3.0;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };
    int  i1, j1, i, j, k, l, err=0;

    double  mu     = matData[0];      // shear modulus
    double  K      = 1.0/matData[1];  // matData[1] is the inverse of bulk modulus

    // gradient of displacement tensor
    MatrixXd  grad = F;
    grad(0,0) -= 1.0;    grad(1,1) -= 1.0;    grad(2,2) -= 1.0;
    // small-strain tensor
    MatrixXd  eps = 0.5*(grad+grad.transpose());

    double  volstrain = eps(0,0) + eps(1,1) + eps(2,2);

    stre.setZero();

    /////////////////////////////////
    // Cauchy stress tensor

    fact1 = 2.0*mu;
    fact2 = (-2.0*mu/3.0)*volstrain;

    if(MIXED_ELEMENT)
      fact2 += K*volstrainnew;
    else
      fact2 += K*volstrain;

    stre[0] = fact1*eps(0,0);    stre[3] = fact1*eps(0,1);    stre[6] = fact1*eps(0,2);
    stre[1] = fact1*eps(1,0);    stre[4] = fact1*eps(1,1);    stre[7] = fact1*eps(1,2);
    stre[2] = fact1*eps(2,0);    stre[5] = fact1*eps(2,1);    stre[8] = fact1*eps(2,2);
    stre[0] += fact2 ;           stre[4] += fact2 ;           stre[8] += fact2 ;

    return err;
}

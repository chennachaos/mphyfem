#include "Matl_ElecMech_Piezoelectric.h"
#include "headersEigen.h"
#include "util.h"
#include "utilitiesmaterial.h"


using namespace std;


Matl_ElecMech_Piezoelectric::Matl_ElecMech_Piezoelectric()
{

}


Matl_ElecMech_Piezoelectric::~Matl_ElecMech_Piezoelectric()
{

}



int Matl_ElecMech_Piezoelectric::computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGP, int gp)
{
    /*
     xx -> 11   yx -> 21  zx -> 31
     xy -> 12   yy -> 22  zy -> 32
     xz -> 13   yz -> 23  zz -> 33
    */

    int isw = 1,  err = 0;

    // check for the validity of input data and assumptions

    /*
    if(finite)
    {
      cerr << " finite eps assumption is not valid in 'small eps elasticity' " << endl;
      cerr << " \n \n Check the input data \n \n " << endl;
      cerr << " \n \n The program has been aborted. \n \n " << endl;

      return -2;
    }
    */

    // compute stress and tangent matrix

    int ii, jj;
    double  fact1, fact2;
    double  K     = matData[1]; // bulk modulus
    double  mu    = matData[2]; // shear modulus
    double  eps0  = matData[3]; // permitivity
    double  R2mu  = 2.0*mu;
    double  Lamb  = K-2.0*mu/3.0;

    VectorXd  eps(9);

    Cmat.setZero();
    Bmat.setZero();
    Amat.setZero();
    stre.setZero();
    elecDisp.setZero();


    if(sss == 1) // plane stress
    {
      /*
      double  E  = 9.0*K*G/(3.0*K+G);
      double  nu = 0.5*(3.0*K-2.0*G)/(3.0*K+G);

      fact1 = E/(1.0-nu*nu);

      // small-eps tensor
      eps[0] = F[0][0] - 1.0;
      eps[1] = F[1][1] - 1.0;
      eps[3] = 0.5*(F[0][1] + F[1][0]);

      // Cauchy stress tensor
      stre[0] = fact1*(   eps[0] + nu*eps[1]);
      stre[1] = fact1*(nu*eps[0] +    eps[1]);
      stre[2] = 0.0;
      stre[3] = R2G*eps[3];
      stre[4] = 0.0;
      stre[5] = 0.0;

      // this is to be done only for the plane stress case
      eps[2]  = -nu*(eps[0]+eps[1])/(1.0-nu); //or
      //eps[2]  = -nu*(stre[0]+stre[1])/E;
      F[2][2] = 1.0+eps[2];

      cc[0][0] =    fact1;      cc[0][1] = fact1*nu;
      cc[1][0] = nu*fact1;      cc[1][1] = fact1;
      cc[3][3] = G;
      */
    }
    // plane-eps and axisymmetric
    else
    {
        // small-eps tensor
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

        //       xx                yx               zx                xy                yy                  zy
        Bmat(0,0) =  0.0;  Bmat(0,1) = 17.0;  Bmat(0,2) = 0.0;  Bmat(0,3) = 17.0;  Bmat(0,4) =  0.0;  Bmat(0,5) =  0.0;
        Bmat(1,0) = -6.5;  Bmat(1,1) = 17.0;  Bmat(1,2) = 0.0;  Bmat(1,3) = 17.0;  Bmat(1,4) = 23.3;  Bmat(1,5) =  0.0;
        //Bmat *=  1.0e-6;
        Bmat *=  1.0e-2;

        //Bmat.setZero();

        Amat(0,0) =  1.5;      Amat(0,1) =  0.0;     Amat(0,2) =  0.0;
        Amat(1,0) =  0.0;      Amat(1,1) =  1.5;     Amat(1,2) =  0.0;
        Amat(2,0) =  0.0;      Amat(2,1) =  0.0;     Amat(2,2) =  1.5;

        /*
        Cmat(0,0) = 12.6;  Cmat(0,1) =  5.3;  Cmat(0,2) = 0.0;  Cmat(0,3) = 0.0;
        Cmat(1,0) =  5.3;  Cmat(1,1) = 12.6;  Cmat(1,2) = 0.0;  Cmat(1,3) = 0.0;
        Cmat(2,0) =  0.0;  Cmat(2,1) =  0.0;  Cmat(2,2) = 0.0;  Cmat(2,3) = 0.0;
        Cmat(3,0) =  0.0;  Cmat(3,1) =  0.0;  Cmat(3,2) = 0.0;  Cmat(3,3) = 3.53;
        //Cmat *=  10.0e4;

        Bmat(0,0) =  0.0;  Bmat(0,1) =  0.0;  Bmat(0,2) = 0.0;  Bmat(0,3) = 17.0;
        Bmat(1,0) = -6.5;  Bmat(1,1) = 23.3;  Bmat(1,2) = 0.0;  Bmat(1,3) =  0.0;
        Bmat *=  10.0e-6;

        //Bmat.setZero();

        Amat(0,0) =  1.5e-6;   Amat(0,1) =  0.0;
        Amat(1,0) =  0.0;      Amat(1,1) =  1.5e-6;
        */

        stre     = Cmat*eps - Bmat.transpose()*elecField;
        elecDisp = Bmat*eps + Amat*elecField;
    }

    //printVector(eps);
    //printVector(stre);
    //printVector(elecDisp);

    return err;
}

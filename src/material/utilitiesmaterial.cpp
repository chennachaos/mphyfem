#include "utilitiesmaterial.h"
#include <iostream>
#include <string>
#include <unordered_map>
#include <stdexcept>

using namespace std;


void   matrix2vector(const MatrixXd& mat,  VectorXd& vec)
{
    vec[0] = mat(0,0);    vec[3] = mat(0,1);    vec[6] = mat(0,2);
    vec[1] = mat(1,0);    vec[4] = mat(1,1);    vec[7] = mat(1,2);
    vec[2] = mat(2,0);    vec[5] = mat(2,1);    vec[8] = mat(2,2);

    return;
}


void   vector2matrix(const VectorXd& vec, MatrixXd& mat)
{
    mat(0,0) = vec[0];    mat(0,1) = vec[3];    mat(0,2) = vec[6];
    mat(1,0) = vec[1];    mat(1,1) = vec[4];    mat(1,2) = vec[7];
    mat(2,0) = vec[2];    mat(2,1) = vec[5];    mat(2,2) = vec[8];

    return;
}


void  get_smallstraintensor(const MatrixXd& F, MatrixXd& eps)
{
    // gradient of displacement tensor
    MatrixXd  grad = F;
    grad(0,0) -= 1.0;    grad(1,1) -= 1.0;    grad(2,2) -= 1.0;

    // small-strain tensor
    eps = 0.5*(grad+grad.transpose());

    return;
}


void  get_deviatoricpart(const MatrixXd& F, MatrixXd& eps)
{
    double  val = F.trace()/3.0;

    eps = F;
    eps(0,0) -= val;    eps(1,1) -= val;    eps(2,2) -= val;

    return;
}


void  defgrad2eps(VectorXd&  F,  VectorXd&  eps)
{
    // gradient of displacement
    eps = F;
    // epsXX
    eps[0] -= 1.0;
    // epsYY
    eps[4] -= 1.0;
    // epsZZ
    eps[8] -= 1.0;

    // epsXY
    eps[1] = 0.5*(eps[1]+eps[3]);
    eps[3] = eps[1];

    // epsYZ
    eps[5] = 0.5*(eps[5]+eps[7]);
    eps[7] = eps[5];

    // epsZX
    eps[2] = 0.5*(eps[2]+eps[6]);
    eps[6] = eps[2];

    return;
}


//  b = F*Ft
void  computeLeftCGtensor(VectorXd& F, VectorXd& b)
{
    b(0) = F(0)*F(0) + F(3)*F(3) + F(6)*F(6);
    b(1) = F(1)*F(0) + F(4)*F(3) + F(7)*F(6);
    b(2) = F(2)*F(0) + F(5)*F(3) + F(8)*F(6);

    b(3) = F(0)*F(1) + F(3)*F(4) + F(6)*F(7);
    b(4) = F(1)*F(1) + F(4)*F(4) + F(7)*F(7);
    b(5) = F(2)*F(1) + F(5)*F(4) + F(8)*F(7);

    b(6) = F(0)*F(2) + F(3)*F(5) + F(6)*F(8);
    b(7) = F(1)*F(2) + F(4)*F(5) + F(7)*F(8);
    b(8) = F(2)*F(2) + F(5)*F(5) + F(8)*F(8);

    return;
}



//  C = Ft*F
void  computeRightCGtensor(VectorXd& F, VectorXd& C)
{
    C(0) = F(0)*F(0) + F(1)*F(1) + F(2)*F(2);
    C(1) = F(3)*F(0) + F(4)*F(1) + F(5)*F(2);
    C(2) = F(6)*F(0) + F(7)*F(1) + F(8)*F(2);

    C(3) = F(0)*F(3) + F(1)*F(4) + F(2)*F(5);
    C(4) = F(3)*F(3) + F(4)*F(4) + F(5)*F(5);
    C(5) = F(6)*F(3) + F(7)*F(4) + F(8)*F(5);

    C(6) = F(0)*F(6) + F(1)*F(7) + F(2)*F(8);
    C(7) = F(3)*F(6) + F(4)*F(7) + F(5)*F(8);
    C(8) = F(6)*F(6) + F(7)*F(7) + F(8)*F(8);

    return;
}



// ADJUST F FOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
// Subroutine for displacement formulation
void  adjust_deformationgradient_2D(int sss, bool finite, MatrixXd&  F)
{
    double  detF = F.determinant();

    if(sss == 1)  // plane stress
    {
        if(finite)
          F(2,2) = 1.0/sqrt(detF);
        else
          F(2,2) = 3.0 - F(0,0) - F(1,1);
    }
    else if(sss == 2)    // plane strain
    {
        F(2,2) = 1.0;
    }

    return;
}



// ADJUST F FOR 2D PROBLEMS BASED ON THE ASSUMPTIONS OF PLANE STRESS/PLANE STRAIN/AXISYMMETRIC
// Subroutine for F-bar formulation
void  adjust_deformationgradient_2D_Fbar(int sss, bool finite, MatrixXd&  F, double detFc)
{
    double  detF = F.determinant();

    if(sss == 1)  // plane stress
    {
        if(finite)
          F(2,2) = 1.0/sqrt(detF);
        else
          F(2,2) = 3.0 - F(0,0) - F(1,1);
    }
    else if(sss == 2)    // plane strain
    {
        F *= sqrt(detFc/detF) ;

        F(2,2) = 1.0;
    }

    return;
}


void getAmtx2D(bool axsy, double cc[4][4], double stre[4], double aa[5][5])
{
  // this routine computes "a" matrix from D matrix
  // and stress vector
  // stre[0] ->s11, stre[1] ->s22, stre[2] ->s33, stre[3] ->s12

    int ii, jj, nn=5, n1, n2, mapping[5] = {0, 3, 3, 1, 2};

    // create A matrix
    for(ii=0; ii<nn; ii++)
    {
      n1 = mapping[ii];
      for(jj=0; jj<nn; jj++)
      {
        n2 = mapping[jj];

        aa[ii][jj] = cc[n1][n2];
      }
    }

    // add stress terms
    aa[0][0] += stre[0]; aa[0][2] += stre[3];
    aa[1][1] += stre[0]; aa[1][3] += stre[3];
    aa[2][0] += stre[3]; aa[2][2] += stre[1];
    aa[3][1] += stre[3]; aa[3][3] += stre[1];

    if(axsy) // axisymmetric case
    {
      aa[4][4] += stre[2];
    }

    return;
}



void getQmtx2D(bool axsy, double aa[5][5], double stre[4], double q[5])
{
  // this routine computes "a" matrix from D matrix
  // and stress vector
  // stre[0] ->s11, stre[1] ->s22, stre[2] ->s33, stre[3] ->s12

    int ii, jj, nn, n1, n2, mapping[3] = {0, 3, 4};
    double  vec1[5], vec2[5], fact1, fact2;

    if(axsy) // axisymmetric case
    {
      nn = 3;

      fact1 = 1.0/3.0;
      fact2 = 2.0/3.0;
    }
    else // 2D case
    {
      nn = 2;

      fact1 = 0.5;
      fact2 = 0.5;
    }

    // compute contribution from a:(IxI)
    for(ii=0; ii<5; ii++)
    {
      vec1[ii] = 0.0;
      for(jj=0; jj<nn; jj++)
      {
        vec1[ii] += aa[ii][mapping[jj]];
      }
    }

    // compute contribution from (sigmaxI)
    vec2[0] = stre[0];
    vec2[1] = stre[3];
    vec2[2] = stre[3];
    vec2[3] = stre[1];
    vec2[4] = stre[2];

    // compute q
    for(ii=0; ii<5; ii++)
      q[ii] = fact1*vec1[ii] - fact2*vec2[ii];

    return;
}




void  addElectricPartToCauchyStress(double  eps0, VectorXd&  elecField,  VectorXd& stre)
{
    double  eFeFd2 = 0.5*(elecField[0]*elecField[0]+elecField[1]*elecField[1]+elecField[2]*elecField[2]);

    stre[0] += eps0*(elecField[0]*elecField[0] - eFeFd2);
    stre[1] += eps0*(elecField[1]*elecField[0]);
    stre[2] += eps0*(elecField[2]*elecField[0]);

    stre[3] += eps0*(elecField[0]*elecField[1]);
    stre[4] += eps0*(elecField[1]*elecField[1] - eFeFd2);
    stre[5] += eps0*(elecField[2]*elecField[1]);

    stre[6] += eps0*(elecField[0]*elecField[2]);
    stre[7] += eps0*(elecField[1]*elecField[2]);
    stre[8] += eps0*(elecField[2]*elecField[2] - eFeFd2);

    return;
}




void  addToCouplingTensor(double  e0, VectorXd&  elecField,  MatrixXd& Bmat)
{
    Bmat(0,0) -=  e0*elecField(0);  Bmat(0,1) +=  e0*elecField(1);  Bmat(0,2) +=  e0*elecField(2);
    Bmat(1,0) -=  e0*elecField(1);  Bmat(1,1) -=  e0*elecField(0);  Bmat(1,2) +=  0.0;
    Bmat(2,0) -=  e0*elecField(2);  Bmat(2,1) +=  0.0;              Bmat(2,2) -=  e0*elecField(0);
    Bmat(3,0) -=  e0*elecField(1);  Bmat(3,1) -=  e0*elecField(0);  Bmat(3,2) +=  0.0;
    Bmat(4,0) +=  e0*elecField(0);  Bmat(4,1) -=  e0*elecField(1);  Bmat(4,2) +=  e0*elecField(2);
    Bmat(5,0) +=  0.0;              Bmat(5,1) -=  e0*elecField(2);  Bmat(5,2) -=  e0*elecField(1);
    Bmat(6,0) -=  e0*elecField(2);  Bmat(6,1) +=  0.0;              Bmat(6,2) -=  e0*elecField(0);
    Bmat(7,0) +=  0.0;              Bmat(7,1) -=  e0*elecField(2);  Bmat(7,2) -=  e0*elecField(1);
    Bmat(8,0) +=  e0*elecField(0);  Bmat(8,1) +=  e0*elecField(1);  Bmat(8,2) -=  e0*elecField(2);

    return;
}



// for the formulation using gradient tensor
void addElectricPartToMaterialTensor(double  eps0, VectorXd&  ev,  MatrixXd&  Cmat)
{
    double  eFeF = ev[0]*ev[0] + ev[1]*ev[1] + ev[2]*ev[2];
    double  fact1, fact2, fact3,  eFeFd2 = 0.5*eFeF;
    int  i1, j1, i, j, k, l;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };

    for(i1=0; i1<9; i1++)
    {
      i = map_tensor4[i1][0];
      j = map_tensor4[i1][1];

      for(j1=0; j1<9; j1++)
      {
        k = map_tensor4[j1][0];
        l = map_tensor4[j1][1];

        // without stress terms
        //Cmat(i1,j1) += eps0*ev(i)*( ev(j)*KronDelta(k,l) - ev(k)*KronDelta(j,l) - ev(l)*KronDelta(j,k) );
        //Cmat(i1,j1) += eps0*eFeFd2*( KronDelta(i,k)*KronDelta(j,l) + KronDelta(i,l)*KronDelta(j,k) - KronDelta(i,j)*KronDelta(k,l) );
        //Cmat(i1,j1) += eps0*( ev[k]*ev[l]*KronDelta(i,j) - ev[k]*ev[j]*KronDelta(i,l) - ev(l)*ev(j)*KronDelta(i,k) );

        // with stress terms
        Cmat(i1,j1) += eps0*ev(i)*( ev(j)*KronDelta(k,l) - ev(k)*KronDelta(j,l) - ev(l)*KronDelta(j,k) );
        Cmat(i1,j1) += eps0*ev(k)*( ev(l)*KronDelta(i,j) - ev(j)*KronDelta(i,l) );
        Cmat(i1,j1) += eps0*eFeFd2*( KronDelta(j,k)*KronDelta(l,i) - KronDelta(j,i)*KronDelta(l,k) );
      }
    }

    return;
}


// for the formulation using gradient tensor
void addMagneticPartToMaterialTensor(double  mu0, VectorXd&  ev,  MatrixXd&  Cmat)
{
    double  eFeF = ev[0]*ev[0] + ev[1]*ev[1] + ev[2]*ev[2];
    double  fact1, fact2, fact3,  eFeFd2 = 0.5*eFeF;
    int  i1, j1, i, j, k, l;

    int  map_tensor4[][2] = { {0,0}, {1,0}, {2,0}, {0,1}, {1,1}, {2,1}, {0,2}, {1,2}, {2,2} };

    for(i1=0; i1<9; i1++)
    {
      i = map_tensor4[i1][0];
      j = map_tensor4[i1][1];

      for(j1=0; j1<9; j1++)
      {
        k = map_tensor4[j1][0];
        l = map_tensor4[j1][1];

        // without stress terms
        //Cmat(i1,j1) += mu0*ev(i)*( ev(j)*KronDelta(k,l) - ev(k)*KronDelta(j,l) - ev(l)*KronDelta(j,k) );
        //Cmat(i1,j1) += mu0*eFeFd2*( KronDelta(i,k)*KronDelta(j,l) + KronDelta(i,l)*KronDelta(j,k) - KronDelta(i,j)*KronDelta(k,l) );
        //Cmat(i1,j1) += mu0*( ev[k]*ev[l]*KronDelta(i,j) - ev[k]*ev[j]*KronDelta(i,l) - ev(l)*ev(j)*KronDelta(i,k) );

        // with stress terms
        Cmat(i1,j1) += mu0*ev(i)*( ev(j)*KronDelta(k,l) - ev(k)*KronDelta(j,l) - ev(l)*KronDelta(j,k) );
        Cmat(i1,j1) += mu0*ev(k)*( ev(l)*KronDelta(i,j) - ev(j)*KronDelta(i,l) );
        Cmat(i1,j1) += mu0*eFeFd2*( KronDelta(j,k)*KronDelta(l,i) - KronDelta(j,i)*KronDelta(l,k) );
      }
    }

    return;
}





/*
// for the formulation using strain tensor
inline  void addElectricPartToMaterialTensor(double  eps0, VectorXd&  elecField,  MatrixXd  Cmat)
{
    double  eFeF = elecField[0]*elecField[0]+elecField[1]*elecField[1]+elecField[2]*elecField[2];
    double  fact1, fact2, fact3;

    // -0.5*eps0*es*es*dij*dkl
    fact1 = -0.5*eps0*eFeF;
    Cmat(0,0) += fact1;  Cmat(0,4) += fact1;  Cmat(0,8) += fact1;
    Cmat(4,0) += fact1;  Cmat(4,4) += fact1;  Cmat(4,8) += fact1;
    Cmat(8,0) += fact1;  Cmat(8,4) += fact1;  Cmat(8,8) += fact1;

    // eps0*es*es*dik*djl
    fact1 = eps0*eFeF;
    for(int ii=0; ii<9; ii++)
      Cmat(ii,ii) += fact1;

    // eps0*ei*ej*dkl
    // eps0*ek*el*dij

    fact1 = eps0*elecField[0]*elecField[0];
    fact2 = eps0*elecField[0]*elecField[1];
    fact3 = eps0*elecField[0]*elecField[2];
    // eps0*ei*ej*dkl
    Cmat(0,0) += fact1;  Cmat(0,4) += fact1;  Cmat(0,8) += fact1;
    Cmat(1,0) += fact2;  Cmat(1,4) += fact2;  Cmat(1,8) += fact2;
    Cmat(2,0) += fact3;  Cmat(2,4) += fact3;  Cmat(2,8) += fact3;
    // eps0*ek*el*dij
    Cmat(0,0) += fact1;  Cmat(4,0) += fact1;  Cmat(8,0) += fact1;
    Cmat(0,1) += fact2;  Cmat(4,1) += fact2;  Cmat(8,1) += fact2;
    Cmat(0,2) += fact3;  Cmat(4,2) += fact3;  Cmat(8,2) += fact3;

    fact1 = eps0*elecField[1]*elecField[0];
    fact2 = eps0*elecField[1]*elecField[1];
    fact3 = eps0*elecField[1]*elecField[2];
    // eps0*ei*ej*dkl
    Cmat(3,0) += fact1;  Cmat(3,4) += fact1;  Cmat(3,8) += fact1;
    Cmat(4,0) += fact2;  Cmat(4,4) += fact2;  Cmat(4,8) += fact2;
    Cmat(5,0) += fact3;  Cmat(5,4) += fact3;  Cmat(5,8) += fact3;
    // eps0*ek*el*dij
    Cmat(0,3) += fact1;  Cmat(4,3) += fact1;  Cmat(8,3) += fact1;
    Cmat(0,4) += fact2;  Cmat(4,4) += fact2;  Cmat(8,4) += fact2;
    Cmat(0,5) += fact3;  Cmat(4,5) += fact3;  Cmat(8,5) += fact3;

    fact1 = eps0*elecField[2]*elecField[0];
    fact2 = eps0*elecField[2]*elecField[1];
    fact3 = eps0*elecField[2]*elecField[2];
    // eps0*ei*ej*dkl
    Cmat(6,0) += fact1;  Cmat(6,4) += fact1;  Cmat(6,8) += fact1;
    Cmat(7,0) += fact2;  Cmat(7,4) += fact2;  Cmat(7,8) += fact2;
    Cmat(8,0) += fact3;  Cmat(8,4) += fact3;  Cmat(8,8) += fact3;
    // eps0*ek*el*dij
    Cmat(0,6) += fact1;  Cmat(4,6) += fact1;  Cmat(8,6) += fact1;
    Cmat(0,7) += fact2;  Cmat(4,7) += fact2;  Cmat(8,7) += fact2;
    Cmat(0,8) += fact3;  Cmat(4,8) += fact3;  Cmat(8,8) += fact3;



    // -2*eps0*ei*el*djk
    // -2*eps0*ek*ej*dil

    fact1 = -2.0*eps0*elecField[0]*elecField[0];
    fact2 = -2.0*eps0*elecField[0]*elecField[1];
    fact3 = -2.0*eps0*elecField[0]*elecField[2];

    // -2*eps0*ei*el*djk
    Cmat(0,0) += fact1;  Cmat(3,1) += fact1;  Cmat(6,2) += fact1;
    Cmat(1,0) += fact2;  Cmat(4,1) += fact2;  Cmat(7,2) += fact2;
    Cmat(2,0) += fact3;  Cmat(5,1) += fact3;  Cmat(8,2) += fact3;
    // -2*eps0*ek*ej*dil
    Cmat(0,0) += fact1;  Cmat(1,3) += fact1;  Cmat(2,6) += fact1;
    Cmat(0,1) += fact2;  Cmat(1,4) += fact2;  Cmat(2,7) += fact2;
    Cmat(0,2) += fact3;  Cmat(1,5) += fact3;  Cmat(2,8) += fact3;


    fact1 = -2.0*eps0*elecField[1]*elecField[0];
    fact2 = -2.0*eps0*elecField[1]*elecField[1];
    fact3 = -2.0*eps0*elecField[1]*elecField[2];

    // -2*eps0*ei*el*djk
    Cmat(0,3) += fact1;  Cmat(3,4) += fact1;  Cmat(6,5) += fact1;
    Cmat(1,3) += fact2;  Cmat(4,4) += fact2;  Cmat(7,5) += fact2;
    Cmat(2,3) += fact3;  Cmat(5,4) += fact3;  Cmat(8,5) += fact3;
    // -2*eps0*ek*ej*dil
    Cmat(3,0) += fact1;  Cmat(4,3) += fact1;  Cmat(5,6) += fact1;
    Cmat(3,1) += fact2;  Cmat(4,4) += fact2;  Cmat(5,7) += fact2;
    Cmat(3,2) += fact3;  Cmat(4,5) += fact3;  Cmat(5,8) += fact3;

    fact1 = -2.0*eps0*elecField[2]*elecField[0];
    fact2 = -2.0*eps0*elecField[2]*elecField[1];
    fact3 = -2.0*eps0*elecField[2]*elecField[2];

    // -2*eps0*ei*el*djk
    Cmat(0,6) += fact1;  Cmat(3,7) += fact1;  Cmat(6,8) += fact1;
    Cmat(1,6) += fact2;  Cmat(4,7) += fact2;  Cmat(7,8) += fact2;
    Cmat(2,6) += fact3;  Cmat(5,7) += fact3;  Cmat(8,8) += fact3;
    // -2*eps0*ek*ej*dil
    Cmat(6,0) += fact1;  Cmat(7,3) += fact1;  Cmat(8,6) += fact1;
    Cmat(6,1) += fact2;  Cmat(7,4) += fact2;  Cmat(8,7) += fact2;
    Cmat(6,2) += fact3;  Cmat(7,5) += fact3;  Cmat(8,8) += fact3;


    return;
}
*/







#ifndef incl_material_utilities_h
#define incl_material_utilities_h

#include <iostream>
#include <cmath>
#include "headersEigen.h"


inline  double  KronDelta(int a, int b)
{
  return ( (a==b) ? 1.0 : 0.0 );
}

inline  double  determinant3b3(VectorXd& F)
{
    return ( F(0)*(F(4)*F(8)-F(5)*F(7)) - F(3)*(F(1)*F(8)-F(7) *F(2)) + F(6)*(F(1)*F(5)-F(2)*F(4)) );
}

void   matrix2vector(const MatrixXd& mat,  VectorXd& vec);


void   vector2matrix(const VectorXd& vec, MatrixXd& mat);


void  adjust_deformationgradient_2D(int sss, bool finite, MatrixXd&  F);


void  adjust_deformationgradient_2D_Fbar(int sss, bool finite, MatrixXd&  F, double detFc);


void  get_smallstraintensor(const MatrixXd& F, MatrixXd& eps);


void  get_deviatoricpart(const MatrixXd& F, MatrixXd& eps);


void  defgrad2eps(VectorXd&  F,  VectorXd&  eps);


//  b = F*Ft
void  computeLeftCGtensor(VectorXd& F, VectorXd& b);


//  C = Ft*F
void  computeRightCGtensor(VectorXd& F, VectorXd& C);


void getAmtx2D(bool axsy, double cc[4][4], double stre[4], double aa[5][5]);


void getQmtx2D(bool axsy, double aa[5][5], double stre[4], double q[5]);


void  addElectricPartToCauchyStress(double  eps0, VectorXd&  elecField,  VectorXd& stre);


void  addToCouplingTensor(double  e0, VectorXd&  elecField,  MatrixXd& Bmat);


// for the formulation using gradient tensor
void addElectricPartToMaterialTensor(double  eps0, VectorXd&  ev,  MatrixXd&  Cmat);


// for the formulation using gradient tensor
void addMagneticPartToMaterialTensor(double  mu0, VectorXd&  ev,  MatrixXd&  Cmat);



int  getMaterialID(std::string& matkey);


int  getElementID_Standard(std::string& matkey);


int  getElementID_MagnetoMech(std::string& matkey);


int  getElementID_ElectroMech(std::string& matkey);


int  getElementID_Growth(std::string& matkey);






#endif

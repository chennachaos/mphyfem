#ifndef incl_Matl_TherMech_Linear_h
#define incl_Matl_TherMech_Linear_h

#include "MaterialBase.h"
#include "headersEigen.h"

class Matl_TherMech_Linear : public MaterialBase
{
  public:

    //member functions

    Matl_TherMech_Linear();

    virtual ~Matl_TherMech_Linear();

    virtual int computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& temp, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp);
};

#endif


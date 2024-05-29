#ifndef incl_Matl_ElecMech_Piezoelectric_h
#define incl_Matl_ElecMech_Piezoelectric_h

#include "MaterialBase.h"
#include "headersEigen.h"

class Matl_ElecMech_Piezoelectric : public MaterialBase
{
  public:

    //member functions

    Matl_ElecMech_Piezoelectric();

    virtual ~Matl_ElecMech_Piezoelectric();

    virtual int computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp);
};

#endif


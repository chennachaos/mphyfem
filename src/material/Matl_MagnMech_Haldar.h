#ifndef incl_Matl_MagnMech_Haldar_h
#define incl_Matl_MagnMech_Haldar_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_MagnMech_Haldar : public MaterialBase
{
  public:

    //member functions

    Matl_MagnMech_Haldar();

    virtual ~Matl_MagnMech_Haldar();

    virtual int getMaterialTypeNameNumber()
    {  return 2005; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre, VectorXd&  magnDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


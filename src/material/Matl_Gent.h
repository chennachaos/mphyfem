#ifndef incl_Matl_Gent_h
#define incl_Matl_Gent_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_Gent : public MaterialBase
{
  public:

    //member functions

    Matl_Gent();

    virtual ~Matl_Gent();

    virtual int getMaterialTypeNameNumber()
    {  return 5; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


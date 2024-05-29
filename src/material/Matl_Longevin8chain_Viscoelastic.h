#ifndef incl_Matl_Longevin8chain_Viscoelastic_h
#define incl_Matl_Longevin8chain_Viscoelastic_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_Longevin8chain_Viscoelastic : public MaterialBase
{
  public:

    //member functions

    Matl_Longevin8chain_Viscoelastic();

    virtual ~Matl_Longevin8chain_Viscoelastic();

    virtual int getMaterialTypeNameNumber()
    {  return 106; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


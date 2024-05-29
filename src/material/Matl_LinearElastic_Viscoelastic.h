#ifndef incl_Matl_LinearElastic_Viscoelastic_h
#define incl_Matl_LinearElastic_Viscoelastic_h

#include "MaterialBase.h"
#include "headersEigen.h"

class Matl_LinearElastic_Viscoelastic : public MaterialBase
{
  public:

    //member functions

    Matl_LinearElastic_Viscoelastic();

    virtual ~Matl_LinearElastic_Viscoelastic();

    virtual int getMaterialTypeNameNumber()
    {  return 101; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


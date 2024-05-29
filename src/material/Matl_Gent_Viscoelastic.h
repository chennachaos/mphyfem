#ifndef incl_Matl_Gent_Viscoelastic_h
#define incl_Matl_Gent_Viscoelastic_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_Gent_Viscoelastic : public MaterialBase
{
  public:

    //member functions

    Matl_Gent_Viscoelastic();

    virtual ~Matl_Gent_Viscoelastic();

    virtual int getMaterialTypeNameNumber()
    {  return 105; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


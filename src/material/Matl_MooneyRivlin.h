#ifndef incl_Matl_MooneyRivlin_h
#define incl_Matl_MooneyRivlin_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_MooneyRivlin : public MaterialBase
{
  public:

    //member functions

    Matl_MooneyRivlin();

    virtual ~Matl_MooneyRivlin();

    virtual int getMaterialTypeNameNumber()
    {  return 4; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


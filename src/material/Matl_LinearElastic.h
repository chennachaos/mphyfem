#ifndef incl_Matl_LinearElastic_h
#define incl_Matl_LinearElastic_h

#include "MaterialBase.h"
#include "headersEigen.h"

class Matl_LinearElastic : public MaterialBase
{
  public:

    //member functions

    Matl_LinearElastic();

    virtual ~Matl_LinearElastic();

    virtual int getMaterialTypeNameNumber()
    {  return 1; }

    virtual bool isFiniteStrain()
    {
        return false;
    }

    virtual  double  getKinv()
    {
        return (3.0*(1.0-2.0*matData[1]))/matData[0];
    }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


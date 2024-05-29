#ifndef incl_Matl_NeoHookean_h
#define incl_Matl_NeoHookean_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_NeoHookean : public MaterialBase
{
  public:

    //member functions

    Matl_NeoHookean();

    virtual ~Matl_NeoHookean();

    virtual int getMaterialTypeNameNumber()
    {  return 3; }

    virtual double computeValue(int sss,  MatrixXd&  F);

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
    
    virtual int getStressAndJacobianForInverseGrowth(int sss,  MatrixXd&  F, MatrixXd&  Fg, VectorXd&  stre, MatrixXd&  Jgp, MatrixXd&  Hgp);
};

#endif


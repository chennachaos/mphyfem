#ifndef incl_Matl_ArrudaBoyce_h
#define incl_Matl_ArrudaBoyce_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_ArrudaBoyce : public MaterialBase
{
  public:

    //member functions

    Matl_ArrudaBoyce();

    virtual ~Matl_ArrudaBoyce();

    virtual int getMaterialTypeNameNumber()
    {  return 6; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);

};

#endif


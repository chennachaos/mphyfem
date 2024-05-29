#ifndef incl_Matl_SaintVenantKirchhoff_h
#define incl_Matl_SaintVenantKirchhoff_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_SaintVenantKirchhoff : public MaterialBase
{
  public:

    //member functions

    Matl_SaintVenantKirchhoff();

    virtual ~Matl_SaintVenantKirchhoff();

    virtual int getMaterialTypeNameNumber()
    {  return 3; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);
};

#endif


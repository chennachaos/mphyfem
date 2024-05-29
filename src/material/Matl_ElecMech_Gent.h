#ifndef incl_Matl_ElecMech_Gent_h
#define incl_Matl_ElecMech_Gent_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_ElecMech_Gent : public MaterialBase
{
  public:

    //member functions

    Matl_ElecMech_Gent();

    virtual ~Matl_ElecMech_Gent();

    virtual double  getPermittivity()
    {return matData[3];}

    virtual int computeStressAndTangent(bool tangFlag, int sss,  VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, double* iv1, double* iv2, double dt, int& nivGp, int gp);

    virtual int computeElectricComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat);

    virtual int computeElectricDisplacement(VectorXd&  elecField, VectorXd&  elecDisp);

    virtual int computeMechanicalStress(VectorXd&  F, double&  pres, VectorXd&  stre);

    virtual int computeMaxwellStress(VectorXd&  elecField, VectorXd&  elecDisp);

    virtual int computeTotalStress(VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre);
};

#endif


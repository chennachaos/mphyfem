#ifndef incl_Matl_ElecMech_NeoHookean_h
#define incl_Matl_ElecMech_NeoHookean_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_ElecMech_NeoHookean : public MaterialBase
{
  public:

    //member functions

    Matl_ElecMech_NeoHookean();

    virtual ~Matl_ElecMech_NeoHookean();

    virtual double  getPermittivity()
    {return matData[2];}

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn,  MatrixXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, MatrixXd& intVarPrev, MatrixXd& intVar, int gp, double dt);

    virtual int computeElectricComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat);

    virtual int computeElectricDisplacement(VectorXd&  elecField, VectorXd&  elecDisp);

    virtual int computeMechanicalStress(VectorXd&  F, double&  pres, VectorXd&  stre);

    virtual int computeMaxwellStress(VectorXd&  elecField, VectorXd&  elecDisp);

    virtual int computeTotalStress(VectorXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre);

};

#endif


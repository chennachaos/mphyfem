#ifndef incl_Matl_MagnMech_NeoHookean_h
#define incl_Matl_MagnMech_NeoHookean_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_MagnMech_NeoHookean : public MaterialBase
{
  public:

    //member functions

    Matl_MagnMech_NeoHookean();

    virtual ~Matl_MagnMech_NeoHookean();

    virtual int getMaterialTypeNameNumber()
    {  return 2001; }

    virtual double  getPermittivity()
    {return matData[2];}

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre, VectorXd&  magnDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMagneticComponents(VectorXd&  magnField, VectorXd&  magnDisp, MatrixXd&  Amat);

    virtual int computeMagneticDisplacement(VectorXd&  magnField, VectorXd&  magnDisp);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);

    virtual int computeMaxwellStress(VectorXd&  magnField, VectorXd&  magnDisp);

    virtual int computeTotalStress(MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre);

};

#endif


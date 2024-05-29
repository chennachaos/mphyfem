#ifndef incl_Matl_MagnMech_NeoHookean_Viscoelastic_h
#define incl_Matl_MagnMech_NeoHookean_Viscoelastic_h

#include "MaterialBase.h"
#include "headersEigen.h"


class Matl_MagnMech_NeoHookean_Viscoelastic : public MaterialBase
{
  public:

    //member functions

    Matl_MagnMech_NeoHookean_Viscoelastic();

    virtual ~Matl_MagnMech_NeoHookean_Viscoelastic();

    virtual int getMaterialTypeNameNumber()
    {  return 2006; }

    virtual double  getPermittivity()
    {return data_magnfield[0];}

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt);

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre, VectorXd&  magnDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt);

    virtual int computeMagneticComponents(VectorXd&  magnField, VectorXd&  magnDisp, MatrixXd&  Amat);

    virtual int computeMagneticDisplacement(VectorXd&  magnField, VectorXd&  magnDisp);

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre);

    virtual int computeMaxwellStress(VectorXd&  magnField, VectorXd&  magnDisp);

    virtual int computeTotalStress(MatrixXd&  F, VectorXd&  magnField, double& pres, VectorXd&  stre);

};

#endif


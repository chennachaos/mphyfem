
#ifndef incl_MagnetoMech3DTet221_h
#define incl_MagnetoMech3DTet221_h


#include "MagnetoMech3DMixedBase221.h"


class  MagnetoMech3DTet221 : public MagnetoMech3DMixedBase221
{
  private:
    int  AlgoType;

  public:

    MagnetoMech3DTet221();

    virtual ~MagnetoMech3DTet221();

    virtual int getElmTypeNameNum()
    {  return 2054; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcStiffnessAndResidualElecField(MatrixXd& Kff, VectorXd& FlocalF, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int  calcInternalForces();

    virtual int  calcLoadVector(VectorXd& Flocal);

    virtual double computeVolume(bool init = false);
};







#endif



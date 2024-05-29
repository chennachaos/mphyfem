
#ifndef incl_MagnetoMech3DWedge21_h
#define incl_MagnetoMech3DWedge21_h


#include "MagnetoMech3DMixedBase21.h"


class  MagnetoMech3DWedge21 : public MagnetoMech3DMixedBase21
{
  private:
    int  AlgoType;

  public:

    MagnetoMech3DWedge21();

    virtual ~MagnetoMech3DWedge21();

    virtual int getElmTypeNameNum()
    {  return 2052; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcResidualPressure(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int calcInternalForces();

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual double computeVolume(bool init);
};







#endif



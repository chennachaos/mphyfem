
#ifndef incl_MagnetoMech2DQua21_h
#define incl_MagnetoMech2DQua21_h


#include "MagnetoMech2DMixedBase21.h"


class  MagnetoMech2DQua21 : public MagnetoMech2DMixedBase21
{
  private:
    int  AlgoType;

  public:

    MagnetoMech2DQua21();

    virtual ~MagnetoMech2DQua21();

    virtual int getElmTypeNameNum()
    {  return 2002; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcResidualPressure(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int  calcInternalForces();

    virtual int  calcLoadVector(VectorXd& Flocal);

    virtual int  solveForPressure();

    virtual double computeVolume(bool init = false);

};







#endif



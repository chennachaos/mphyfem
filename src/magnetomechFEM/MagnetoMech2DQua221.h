
#ifndef incl_MagnetoMech2DQua221_h
#define incl_MagnetoMech2DQua221_h


#include "MagnetoMech2DMixedBase221.h"


class  MagnetoMech2DQua221 : public MagnetoMech2DMixedBase221
{
  private:
    int  AlgoType;

  public:

    MagnetoMech2DQua221();

    virtual ~MagnetoMech2DQua221();

    virtual int getElmTypeNameNum()
    {  return 2004; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcStiffnessAndResidualElecField(MatrixXd& Kff, VectorXd& FlocalF, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcResidualPressure(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int  calcInternalForces();

    virtual int  calcLoadVector(VectorXd& Flocal);

    virtual int  solveForPressure();

    virtual int  solveForPressureTIC(VectorXd& matK1, double beta, double dt, VectorXd& Flocal);

    virtual double computeVolume(bool init = false);

};







#endif



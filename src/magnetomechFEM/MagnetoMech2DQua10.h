
#ifndef incl_MagnetoMech2DQua10_h
#define incl_MagnetoMech2DQua10_h


#include "ElementBase.h"


class  MagnetoMech2DQua10 : public ElementBase
{
  private:
    int  AlgoType;

  public:

    MagnetoMech2DQua10();

    virtual ~MagnetoMech2DQua10();

    virtual int getElmTypeNameNum()
    {  return 2010; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidual(MatrixXd& Kuu, VectorXd& FlocalU, bool firstIter=false);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcResidualPressure(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int  calcInternalForces();

    virtual int  calcLoadVector(VectorXd& Flocal);

    virtual int  solveForPressure();

    virtual double computeVolume(bool init = false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);

    virtual int calcError(int index);
};







#endif



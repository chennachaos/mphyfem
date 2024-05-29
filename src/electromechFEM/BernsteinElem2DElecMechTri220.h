
#ifndef incl_BernsteinElem2DElecMechTri220_h
#define incl_BernsteinElem2DElecMechTri220_h


#include "ElementBase.h"


class  BernsteinElem2DElecMechTri220 : public ElementBase
{
  private:
    int  AlgoType;

  public:

    BernsteinElem2DElecMechTri220();

    virtual ~BernsteinElem2DElecMechTri220();

    virtual int getElmTypeNameNum()
    {  return 1005; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcResidualPressure(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int  toComputeInfSupCondition(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);

    virtual int  calcInternalForces();

    virtual int  calcLoadVector(VectorXd& Flocal);

    virtual int  solveForPressure();

    virtual int  solveForPressureTIC(VectorXd& matK1, double beta, double dt, VectorXd& Flocal);

    virtual double computeVolume(bool init = false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif



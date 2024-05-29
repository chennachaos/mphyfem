
#ifndef incl_GrowthElem2DQua21_h
#define incl_GrowthElem2DQua21_h


#include "GrowthElem2DMixedBase.h"


class  GrowthElem2DQua21 : public GrowthElem2DMixedBase
{
  private:
    int  AlgoType;

  public:

    GrowthElem2DQua21();

    virtual ~GrowthElem2DQua21();

    virtual int getElmTypeNameNum()
    {  return 6003; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcResidualPressure(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int  calcInternalForces();

    virtual int calcForceVectorGrowthModel(VectorXd& Flocal1, VectorXd& Flocal2);

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



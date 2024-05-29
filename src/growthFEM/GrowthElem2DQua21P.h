
#ifndef incl_GrowthElem2DQua21P_h
#define incl_GrowthElem2DQua21P_h


#include "ElementBase.h"


class  GrowthElem2DQua21P : public ElementBase
{
  private:
    int  AlgoType;
    MatrixXd  KppInv;

  public:

    GrowthElem2DQua21P();

    virtual ~GrowthElem2DQua21P();

    virtual int getElmTypeNameNum()
    {  return 6005; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidual(MatrixXd& Kuu, VectorXd& FlocalU, bool firstIter=false);

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



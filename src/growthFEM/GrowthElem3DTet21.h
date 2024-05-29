
#ifndef incl_GrowthElem3DTet21_h
#define incl_GrowthElem3DTet21_h


#include "GrowthElem3DMixedBase.h"


class  GrowthElem3DTet21 : public GrowthElem3DMixedBase
{
  private:
    int  AlgoType;

  public:

    GrowthElem3DTet21();

    virtual ~GrowthElem3DTet21();

    virtual int getElmTypeNameNum()
    {  return 6055; }

    virtual  void prepareElemData();

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcResidualPressure(VectorXd& Flocal);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int calcInternalForces();

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual double computeVolume(bool init);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif



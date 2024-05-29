
#ifndef incl_GrowthElem2DQua1_h
#define incl_GrowthElem2DQua1_h


#include "ElementBase.h"


class  GrowthElem2DQua1 : public ElementBase
{
  public:

    GrowthElem2DQua1();

    virtual ~GrowthElem2DQua1();

    virtual int getElmTypeNameNum()
    {  return 6001; }

    virtual void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual double calcCriticalTimeStep(bool flag);

    virtual int calcInternalForces();

    virtual int calcForceVectorGrowthModel(VectorXd& Flocal1, VectorXd& Flocal2);

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual double computeVolume(bool init = false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);

    virtual int calcError(int index);
};







#endif



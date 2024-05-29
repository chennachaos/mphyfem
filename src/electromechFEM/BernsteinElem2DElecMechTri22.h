
#ifndef incl_BernsteinElem2DElecMechTri22_h
#define incl_BernsteinElem2DElecMechTri22_h


#include "ElementBase.h"


class  BernsteinElem2DElecMechTri22 : public ElementBase
{
  public:

    BernsteinElem2DElecMechTri22();

    virtual ~BernsteinElem2DElecMechTri22();

    virtual int getElmTypeNameNum()
    {  return 1002; }

    virtual void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int calcInternalForces();

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual double computeVolume(bool init = false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif



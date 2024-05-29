
#ifndef incl_BernsteinElem2DElecMechQua11F_h
#define incl_BernsteinElem2DElecMechQua11F_h


#include "ElementBase.h"


class  BernsteinElem2DElecMechQua11F : public ElementBase
{
  public:

    BernsteinElem2DElecMechQua11F();

    virtual ~BernsteinElem2DElecMechQua11F();

    virtual int getElmTypeNameNum()
    {  return 1009; }

    virtual void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter);

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



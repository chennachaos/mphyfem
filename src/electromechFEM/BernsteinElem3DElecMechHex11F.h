
#ifndef incl_BernsteinElem3DElecMechHex11F_h
#define incl_BernsteinElem3DElecMechHex11F_h


#include "ElementBase.h"


class  BernsteinElem3DElecMechHex11F : public ElementBase
{
  public:

    BernsteinElem3DElecMechHex11F();

    virtual ~BernsteinElem3DElecMechHex11F();

    virtual int getElmTypeNameNum()
    {  return 1049;     }

    virtual void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

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



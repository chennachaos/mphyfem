
#ifndef incl_BernsteinElem2DElecMechTri22F_h
#define incl_BernsteinElem2DElecMechTri22F_h


#include "ElementBase.h"


class  BernsteinElem2DElecMechTri22F : public ElementBase
{
  public:

    BernsteinElem2DElecMechTri22F();

    virtual ~BernsteinElem2DElecMechTri22F();

    virtual int getElmTypeNameNum()
    {  return 1003; }

    virtual  void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter=false);
    int calcStiffnessAndResidualSS(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter=false);
    int calcStiffnessAndResidualFS1(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter=false);
    int calcStiffnessAndResidualFS2(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);
    virtual int  calcResidualSS(VectorXd& Flocal);
    virtual int  calcResidualFS1(VectorXd& Flocal);
    virtual int  calcResidualFS2(VectorXd& Flocal);

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual double  calcCriticalTimeStep(bool flag);

    virtual int  calcInternalForces();

    virtual int  calcLoadVector(VectorXd& Flocal);

    virtual double computeVolume(bool init = false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif



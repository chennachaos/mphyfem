
#ifndef incl_BernsteinElem3DElecMechHex221P_h
#define incl_BernsteinElem3DElecMechHex221P_h


#include "ElementBase.h"


class  BernsteinElem3DElecMechHex221P : public ElementBase
{
  private:
    int  AlgoType;
    MatrixXd  KppInv;
    VectorXd  presDOF;

  public:

    BernsteinElem3DElecMechHex221P();

    virtual ~BernsteinElem3DElecMechHex221P();

    virtual int getElmTypeNameNum()
    {  return 1053; }

    virtual  void prepareElemData();

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, VectorXd& FlocalU, VectorXd& FlocalF, bool firstIter=false);

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual int  calcStiffnessAndResidualElecField(MatrixXd& Kff, VectorXd& FlocalF, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

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
};







#endif



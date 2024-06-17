
#ifndef incl_GrowthElem3DHex1_h
#define incl_GrowthElem3DHex1_h


#include "ElementBase.h"


class  GrowthElem3DHex1 : public ElementBase
{
  public:

    GrowthElem3DHex1();

    virtual ~GrowthElem3DHex1();

    virtual int getElmTypeNameNum()
    {  return 6051; }

    virtual void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcResidual(VectorXd& Flocal);

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




#ifndef incl_BernsteinElem2DEdge3Node_h
#define incl_BernsteinElem2DEdge3Node_h


#include "ElementBase.h"


class  BernsteinElem2DEdge3Node : public ElementBase
{
  public:

    BernsteinElem2DEdge3Node();

    virtual ~BernsteinElem2DEdge3Node();

    virtual int getElmTypeNameNum()
    {  return 1082; }

    virtual void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

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



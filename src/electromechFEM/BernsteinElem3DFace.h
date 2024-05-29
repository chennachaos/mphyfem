
#ifndef incl_BernsteinElem3DFace_h
#define incl_BernsteinElem3DFace_h


#include "ElementBase.h"


class  BernsteinElem3DFace : public ElementBase
{
  public:

    BernsteinElem3DFace();

    virtual ~BernsteinElem3DFace();

    virtual int getElmTypeNameNum()
    {  return 2083; }

    virtual void prepareElemData();

    virtual int calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal, bool firstIter=false);

    virtual int  calcResidual(VectorXd& Flocal);

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual double computeVolume(bool init = false);
};







#endif



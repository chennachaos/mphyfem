
#ifndef incl_GrowthElem3DMixedBase_h
#define incl_GrowthElem3DMixedBase_h


#include "ElementBase.h"


class  GrowthElem3DMixedBase : public ElementBase
{
  public:

    GrowthElem3DMixedBase();

    virtual ~GrowthElem3DMixedBase();

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);
};







#endif



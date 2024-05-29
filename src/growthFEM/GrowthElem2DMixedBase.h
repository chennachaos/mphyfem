
#ifndef incl_GrowthElem2DMixedBase_h
#define incl_GrowthElem2DMixedBase_h


#include "ElementBase.h"


class  GrowthElem2DMixedBase : public ElementBase
{
  public:

    GrowthElem2DMixedBase();

    virtual ~GrowthElem2DMixedBase();

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);
};







#endif



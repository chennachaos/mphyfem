
#ifndef incl_MagnetoMech2DMixedBase21_h
#define incl_MagnetoMech2DMixedBase21_h


#include "ElementBase.h"


class  MagnetoMech2DMixedBase21 : public ElementBase
{
  public:

    MagnetoMech2DMixedBase21();

    virtual ~MagnetoMech2DMixedBase21();

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);

    virtual int calcError(int index);
};







#endif



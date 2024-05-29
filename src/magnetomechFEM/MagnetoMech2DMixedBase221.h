
#ifndef incl_MagnetoMech2DMixedBase221_h
#define incl_MagnetoMech2DMixedBase221_h


#include "ElementBase.h"


class  MagnetoMech2DMixedBase221 : public ElementBase
{
  public:

    MagnetoMech2DMixedBase221();

    virtual ~MagnetoMech2DMixedBase221();

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalF, VectorXd& FlocalP, bool firstIter=false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif




#ifndef incl_MagnetoMech3DMixedBase21_h
#define incl_MagnetoMech3DMixedBase21_h


#include "ElementBase.h"


class  MagnetoMech3DMixedBase21 : public ElementBase
{
  public:

    MagnetoMech3DMixedBase21();

    virtual ~MagnetoMech3DMixedBase21();

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2, bool firstIter=false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif



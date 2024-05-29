
#ifndef incl_Elem_Magnmech_3D_SM_221_h
#define incl_Elem_Magnmech_3D_SM_221_h


#include "ElementBase.h"


class  Elem_Magnmech_3D_SM_221 : public ElementBase
{
  public:

    Elem_Magnmech_3D_SM_221();

    virtual ~Elem_Magnmech_3D_SM_221();

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalF, VectorXd& FlocalP, bool firstIter=false);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif



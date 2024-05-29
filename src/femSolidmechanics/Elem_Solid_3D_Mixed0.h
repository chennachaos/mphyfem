
#ifndef incl_Elem_Solid_3D_Mixed0_h
#define incl_Elem_Solid_3D_Mixed0_h


#include "ElementBase.h"


class  Elem_Solid_3D_Mixed0 : public ElementBase
{
  public:

    Elem_Solid_3D_Mixed0();

    virtual ~Elem_Solid_3D_Mixed0();

    virtual int getElmTypeNameNum()
    {
      return ELEM_SOLID_3D_MIXED0;
    }

    virtual void prepareElemData();

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal);

    virtual int  solveForPressure();

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);

};







#endif



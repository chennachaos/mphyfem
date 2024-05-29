
#ifndef incl_Elem_Solid_3D_Disp_h
#define incl_Elem_Solid_3D_Disp_h


#include "ElementBase.h"


class  Elem_Solid_3D_Disp : public ElementBase
{
  public:

    Elem_Solid_3D_Disp();

    virtual ~Elem_Solid_3D_Disp();

    virtual int getElmTypeNameNum()
    {
      return ELEM_SOLID_3D_DISP;
    }

    virtual void prepareElemData();

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidual(MatrixXd& Kuu, VectorXd& Flocal);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);

};







#endif



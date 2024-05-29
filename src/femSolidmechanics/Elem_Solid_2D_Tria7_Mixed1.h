
#ifndef incl_Elem_Solid_2D_Tria7_Mixed1_h
#define incl_Elem_Solid_2D_Tria7_Mixed1_h


#include "ElementBase.h"


class  Elem_Solid_2D_Tria7_Mixed1 : public ElementBase
{
  public:

    Elem_Solid_2D_Tria7_Mixed1();

    virtual ~Elem_Solid_2D_Tria7_Mixed1();

    virtual int getElmTypeNameNum()
    {
      return ELEM_SOLID_2D_TRIA7_MIXED1;
    }

    virtual void prepareElemData();

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalP);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);

};







#endif



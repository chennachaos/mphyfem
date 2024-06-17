
#ifndef incl_Elem_Growth_3D_Mixed21_h
#define incl_Elem_Growth_3D_Mixed21_h


#include "ElementBase.h"


class  Elem_Growth_3D_Mixed21 : public ElementBase
{
  public:

    Elem_Growth_3D_Mixed21();

    virtual ~Elem_Growth_3D_Mixed21();

    virtual int getElmTypeNameNum()
    {
      return ELEM_GROWTH_3D_MIXED21;
    }

    virtual void prepareElemData();

    virtual int calcMassMatrix(MatrixXd& Mlocal, bool MassLumping);

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2);

    virtual void elementContourplot(int, int, int);

    virtual void projectToNodes(bool, int, int, int);

    virtual void projectStrain(bool, int, int, int, double*);

    virtual void projectStress(bool, int, int, int, double*);

    virtual void projectInternalVariable(bool, int, int, int, double*);
};







#endif



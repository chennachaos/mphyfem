
#ifndef incl_Elem_Magnmech_3D_HM_21_h
#define incl_Elem_Magnmech_3D_HM_21_h


#include "ElementBase.h"


class  Elem_Magnmech_3D_HM_21 : public ElementBase
{
  public:

    Elem_Magnmech_3D_HM_21();

    virtual ~Elem_Magnmech_3D_HM_21();

    virtual int getElmTypeNameNum()
    {
      return ELEM_MAGNMECH_3D_HM_21;
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



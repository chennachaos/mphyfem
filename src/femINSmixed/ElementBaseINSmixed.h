#ifndef incl_ElementBaseINSmixed_h
#define incl_ElementBaseINSmixed_h

#include "ElementBase.h"
#include "headersEigen.h"
#include <vector>
#include "AABB.h"
#include "FluidMaterialBase.h"


using std::vector;
using std::cout;
using namespace myGeom;



class ElementBaseINSmixed : public ElementBase
{
  public:

    //member variables
    double  elemVol, charlen, dtCrit;

    vector<double>  elemVolGP;

    vector<VectorXd>  Nv, dNvdx, dNvdy, dNvdz, Np, dNpdx, dNpdy, dNpdz;

    AABB  bbox;

    FluidMaterialBase  *FluidMatlData;

    //member functions

    ElementBaseINSmixed();

    virtual ~ElementBaseINSmixed();

/*
    virtual void calcBoundingBox()
    { cout << "   'calcBoundingBox' is not defined for this element!\n\n"; return; }

    virtual bool isPointInside(myPoint& target_point)
    { cout << "   'isPointInside' is not defined for this element!\n\n"; return 0; }

    virtual void findLocalCoordinates(vector<myPoint>& nodeCoords, myPoint& coords_global, myPoint& coords_local)
    { cout << "   'findLocalCoordinates' is not defined for this element!\n\n"; return; }

    virtual int calcLoadVector(VectorXd& Flocal)
    { cout << "   'calcLoadVector' is not defined for this element!\n\n"; return 0; }

    virtual double  calcCriticalTimeStep(double* elemData, double* timeData, VectorXd&  veloVec)
    { cout << "   'calcCriticalTimeStep' is not defined for this element!\n\n"; return 0.0; }

    virtual  int calcError(int index)
    { cout << "  'calcError' is not available for this element!\n\n"; return -1; }

    virtual  void  prepareElemData()
    { cout << "   'prepareElemData' is not defined for this element!\n\n"; return; }

    virtual void prepareElemData(vector<myPoint>& nodeCoords)
    { cout << "   'prepareElemData' is not defined for this element!\n\n"; return; }

    virtual int MassMatrices(vector<myPoint>& nodeCoords, double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2)
    { cout << "   'MassMatrices' is not defined for this element!\n\n"; return -1; }

    virtual double  ResidualIncNavStokesAlgo1(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& veloDotVec, VectorXd& veloDotVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  Flocal1, VectorXd&  Flocal2, double timeCur)
    { cout << "   'ResidualIncNavStokesAlgo1' is not defined for this element!\n\n"; return -1; }

    virtual double  ResidualIncNavStokesSemiImpl(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& veloDotVec, VectorXd& veloDotVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  Flocal1, VectorXd&  Flocal2, double timeCur)
    { cout << "   'ResidualIncNavStokesSemiImpl' is not defined for this element!\n\n"; return -1; }

    virtual int  ResidualIncNavStokesAlgo2(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloCur, VectorXd& acceCur, VectorXd& presCur, VectorXd&  Flocal2)
    { cout << "   'ResidualIncNavStokesAlgo2' is not defined for this element!\n\n"; return -1; }

    virtual double CalculateError(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloDotPrev, VectorXd& presPrev, double timeCur, int index)
    { cout << "   'CalculateError' is not defined for this element!\n\n"; return 0; }

    virtual int  StiffnessAndResidualFullyImplicit(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& velo, VectorXd& veloPrev, VectorXd& veloCur, VectorXd& veloDotCur, VectorXd& presCur, MatrixXd& Kuu, MatrixXd& Kup, VectorXd& Fu, VectorXd& Fp, double dt, double timeCur)
    { cout << "   'StiffnessAndResidual' is not defined for this element!\n\n"; return -1; }

    virtual int  StiffnessAndResidualFullyImplicit(vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& dispCur, VectorXd& veloCur, VectorXd& acceCur, VectorXd& dispPrev, VectorXd& veloPrev, VectorXd& accePrev, MatrixXd& Klocal, VectorXd& Flocal, double dt, double timeCur)
    { cout << "   'StiffnessAndResidual' is not defined for this element!\n\n"; return -1; }

    virtual int  StiffnessForSemiImpl(double* elemData, double* timeData, MatrixXd& Kup)
    { cout << "   'StiffnessForSemiImpl' is not defined for this element!\n\n"; return -1; }

    virtual int toComputeInfSupCondition(vector<myPoint>& nodeCoords, double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp)
    { cout << "   'toComputeInfSupCondition' is not defined for this element!\n\n"; return -1; }

    virtual int CalculateForces(int side, vector<myPoint>& nodeCoords, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& presVec, VectorXd&  Flocal1)
    { cout << "   'CalculateForces' is not defined for this element!\n\n"; return -1; }
*/
};



#endif //incl_Lagrange_Element_h


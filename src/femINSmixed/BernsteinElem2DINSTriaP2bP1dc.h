
#ifndef incl_BernsteinElem2DINSTriaP2bP1dc_h
#define incl_BernsteinElem2DINSTriaP2bP1dc_h

#include "ElementBaseINSmixed.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "SolutionData.h"


class  BernsteinElem2DINSTriaP2bP1dc : public ElementBaseINSmixed
{
  public:

    BernsteinElem2DINSTriaP2bP1dc();

    virtual ~BernsteinElem2DINSTriaP2bP1dc();

    void prepareElemData(vector<myPoint>& nodeCoords);

    //virtual bool isPointInside(myPoint& target_point);

    virtual void findLocalCoordinates(vector<myPoint>& node_coods, myPoint& coords_global, myPoint& coords_local);

    virtual double calcCriticalTimeStep(double* elemData, double* timeData, VectorXd&  veloVec);

    virtual int calcLoadVector(VectorXd& Flocal);

    virtual int MassMatrices(vector<myPoint>& node_coods, double* elemData, VectorXd&  Mlocal1, VectorXd&  Mlocal2);

    virtual double  ResidualIncNavStokesAlgo1(vector<myPoint>& node_coods, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& veloDotVec, VectorXd& veloDotVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  Flocal1, VectorXd&  Flocal2, double timeCur);

    virtual double  ResidualIncNavStokesSemiImpl(vector<myPoint>& node_coods, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& veloVecPrev, VectorXd& veloDotVec, VectorXd& veloDotVecPrev, VectorXd& presVec, VectorXd& presVecPrev, VectorXd&  Flocal1, VectorXd&  Flocal2, double timeCur);

    virtual int  ResidualIncNavStokesAlgo2(vector<myPoint>& node_coods, double* elemData, double* timeData, VectorXd& veloCur, VectorXd& acceCur, VectorXd& presCur, VectorXd&  Flocal2);

    virtual double CalculateError(vector<myPoint>& node_coords, double* elemData, double* timeData, VectorXd& veloPrev, VectorXd& veloDotPrev, VectorXd& presPrev, double timeCur, int index);

    virtual int  StiffnessAndResidualFullyImplicit(vector<myPoint>& node_coords, double* elemData, double* timeData, VectorXd& velo, VectorXd& veloPrev, VectorXd& veloCur, VectorXd& veloDotCur, VectorXd& presCur, MatrixXd& Kuu, MatrixXd& Kup, VectorXd& Fu, VectorXd& Fp);

    virtual int  StiffnessForSemiImpl(double* elemData, double* timeData, MatrixXd& Kup);

    virtual int toComputeInfSupCondition(vector<myPoint>& node_coords, double* elemData, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);

    virtual int CalculateForces(int side, vector<myPoint>& node_coords, double* elemData, double* timeData, VectorXd& veloVec, VectorXd& presVec, VectorXd&  Flocal1);
};







#endif



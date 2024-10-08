#ifndef incl_ImmersedIntegrationElement_h
#define incl_ImmersedIntegrationElement_h


#include "util.h"
#include "ElementBase.h"

using namespace std;
using namespace Eigen;

class SolutionData;
class GeomDataLagrange;
class GeomDataHBSplines;



class ImmersedIntegrationElement
{
    private:

        static  int  itemcount;

        int ndim, nnode, id, nGP;

        bool  IS_ACTIVE;

    public:

        vector<int>  nodeNums, forAssyVec, elemNums;
        vector<double>  gausspoints, gaussweights;

        myPoint coords, param;

        VectorXd  Flocal, Flocal2;

        MatrixXd  Khorz, Kvert;

        vector<ElementBase*>  elems;

        SolutionData  *SolnData;
        GeomDataLagrange  *GeomDataLag;
        //GeomDataHBSplines *GeomDataHBS;

        ImmersedIntegrationElement();

        ~ImmersedIntegrationElement();

        int getID()
        {  return id; }

        static int  getCount()
        { return itemcount; }

        void setDimension(int dd)
        {  ndim = dd; return;  }

        int getDimension()
        {  return ndim;  }

        bool isActive()
        {  return (IS_ACTIVE == true);  }

        void deactivate()
        {  IS_ACTIVE = false;  return ;  }

        void activate()
        {  IS_ACTIVE = true;  return ;  }

        void  printSelf();

        void  findElements();

        void  initialiseDOFvalues();

        void  resetMatrixAndVector();
        
        void  computePointAtGP(int ind, myPoint& pt);

        void  integrateForceAndMoment(VectorXd& vectemp, myPoint& centLoc);

        void  integrateForceFlexible(int ind1, int ind2, VectorXd& vectemp);

        void  computeKhorzKvertRigid(MatrixXd& tempMat, MatrixXd& tempMat2);

        void  computeKhorzKvertFlexible(int ind1, int ind2, MatrixXd& tempMat, MatrixXd& tempMat2);

        void  calcStiffnessAndResidual(int ind1=0, int ind2=0, double inp1=0.0, double inp2=0.0);

        void  assembleElementMatrix(int, SparseMatrixXd&);

        void  assembleElementVector(int ind, bool flag, double* rhs);

        void  assembleMatrixAndVector(int, SparseMatrixXd&, double*);

        void  prepareElemData();

        void  getDataFromGlobalSolutionVector(bool flag);

        void  reset();

        void  computeVelocity(const VectorXd& NN, myPoint&  velSpec);

        void  computeVelocityCur(const VectorXd& NN, myPoint&  velSpec);

        void  computeAcceleration(const VectorXd& NN, myPoint&  velSpec);

        void  computeAccelerationCur(const VectorXd& NN, myPoint&  velSpec);
};






#endif



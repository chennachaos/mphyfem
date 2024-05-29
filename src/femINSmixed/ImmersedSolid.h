#ifndef incl_ImmersedSolid_h
#define incl_ImmersedSolid_h

#include "SolutionData.h"
#include "GeomDataLagrange.h"
#include "AABB.h"
#include "headersVTK.h"
#include "SolverEigen.h"


using namespace myGeom;

class  ImmersedIntegrationElement;


enum  {BC_ENFORCE_TYPE_LAGRANGE=0, BC_ENFORCE_TYPE_PENALTY};


class ImmersedSolid
{
    public:

        static  int  count;

        int  id, ndim, fileCount;
        int  nNode, nElem, npElem, totalDOF;
        int  BC_ENFORCE_TYPE;

        double  tol, PENALTY, NitscheFact;

        bool  isNitsche;

        vector<myPoint>  nodePosOrig, nodePosCur, specVelocity;

        vector<vector<int> >  OutputData;
        vector<vector<double> > rigidBodyMotionLimits;
        //vector<TimeFunctionCore> PrescMotionTimeFuncs;

        VectorXd  totalForce;
        myPoint  centroid;

        AABB  bbox;

        vector<ImmersedIntegrationElement*>  ImmIntgElems;

        //vector<myPoly*>  ImmersedFaces;

        vtkSmartPointer<vtkUnstructuredGrid>  uGrid;

        vtkSmartPointer<vtkPolyData>  polyDataVTK;

        //////////////////////////////////////////////////////
        // member functions
        //////////////////////////////////////////////////////

        ImmersedSolid();

        virtual ~ImmersedSolid();

        int getID()
        {  return id; }

        void  setDimension(int dd)
        { ndim = dd; }

        void  setNumberOfNodes(int nn)
        {  nNode = nn; }

        int  getNumberOfNodes()
        { return nNode; } 

        void  setNumberOfElements(int nn)
        {  nElem = nn; }

        int  getNumberOfElements()
        {  return nElem; }

        int getTotalDOF()
        { return  totalDOF; }

        void  setPenaltyParameter(double tt)
        { PENALTY = tt; }

        double  getPenaltyParameter()
        { return PENALTY; }

        void  setNitscheFact(double tt)
        { NitscheFact = tt; }

        double  getNitscheFact()
        { return NitscheFact; }

        void  setNitscheFlag(bool tt)
        { isNitsche = tt; }

        bool  getNitscheFlag()
        { return isNitsche; }

       
        void  readInput(string& fname);
        
        void  readDimension(ifstream& infile, string& line);

        void  readNodes(ifstream& infile, string& line);

        void  readImmersedBoundaryElements(ifstream& infile, string& line);

        void  readPrescribedMotion(ifstream& infile, string& line);

        void  setBoundaryConditionType(int tt)
        {
          if(tt)
            BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_LAGRANGE;
          else
            BC_ENFORCE_TYPE = BC_ENFORCE_TYPE_PENALTY;
        }

        bool isBoundaryConditionTypeLagrange()
        {
          return (BC_ENFORCE_TYPE == BC_ENFORCE_TYPE_LAGRANGE);
        }

        void  adjustBoundaryPoints(double* minVal, double* maxVal);

        void  setDataForOutput(vector<vector<int> >& vectemp);

        void  printSelf();

        void  reset();

        void  computePointAtGP(int elnum, double gpcoord, myPoint& pt);

        void  computeCentroid(int index);

        myPoint&  getCentroid(int index)
        {
          computeCentroid(index);

          return  centroid;
        }

        void  computeAABB(int index);

        virtual  int  within(myPoint& pt);

        virtual  int  doAABBintersect(AABB&  bb2);

        virtual  int  doIntersect2D(AABB&  bb2, bool  flag, vector<int>&  vecTemp, vector<myPoint>&  ptOut);

        virtual  int  doIntersect2Dfor3D(int sideTemp, double coord3, AABB& bbTemp, bool flag, vector<int>&  vecTemp,  vector<myPoint>& ptOut);

        virtual  int  doIntersect3D(AABB&  bb2, bool  flag, vector<int>&  vecTemp, vector<myPoint>&  ptOut);

        double  distanceFromPoint(myPoint&  pt);

        double  distanceFromPoint(double xx=0.0, double yy=0.0, double zz=0.0);

        void  computeTotalForce();

        virtual void  writeOutput()
        { cout << "   'writeOutput()' is not defined for this Solid!\n\n"; return; }

        void  postProcess(string infilename, int fileCount);
        //{ cout << "   'postProcess()' is not defined for this Solid!\n\n"; return; }

        int  updatePointPositions(double timeNow);
        
        int  updatePointPositions(int ndof, VectorXd& disp);

        int  updatePointVelocities(int ndof, VectorXd& velo);

        //{ cout << "   'updatePointPositions()' is not defined for this Solid!\n\n"; return 0; }

        virtual  void  writeResult(ofstream&)
        { cout << " writeResult() ... is not defined for the class... " << endl; return; }

        virtual  void  readResult(ifstream&)
        { cout << " readResult() ... is not defined for the class... " << endl; return; }
};






#endif

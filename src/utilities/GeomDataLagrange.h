
#ifndef incl_GeomDataLagrange_h
#define incl_GeomDataLagrange_h

#include "util.h"
#include <vector>
#include "headersBasic.h"
#include "headersEigen.h"

using std::vector;


class  GeomDataLagrange
{
  private:
    int  ndim, nNode, bodyForceTimeFunction;
    FemModelBehaviour  modelbehaviour;
    double  thickness;
    myPoint  bodyForce;

  public:

    vector<myPoint>  NodePosOrig, NodePosNew, NodePosCur, NodePosPrev, NodePosPrevCur;
    vector<myPoint>  specValNew, specValCur, acceNew, acceCur;

    vector<int>  node_map_get_old, node_map_get_new;
    vector<int>  assy4r;

    GeomDataLagrange();

    ~GeomDataLagrange();


    void  setDimension(int dd)
    {      ndim = dd;    }

    int  getNumNode()
    {  return      nNode;    }

    void  setModelBehaviour(FemModelBehaviour sss)
    {
      modelbehaviour = sss;
      return;
    }

    FemModelBehaviour  getModelBehaviour()
    {
      return modelbehaviour;
    }

    void  setThicnkess(double thk)
    {
      thickness = thk;
      return;
    }

    double  getThickness()
    {
      return thickness;
    }

    void  setBodyForce(myPoint& bforce)
    {
      bodyForce = bforce;
      return;
    }

    myPoint  getBodyForce()
    {
      return bodyForce;
    }

    void  setBodyForceTimeFunction(int bftf)
    {
      bodyForceTimeFunction = bftf;
      return;
    }

    int  getBodyForceTimeFunction()
    {
      return bodyForceTimeFunction;
    }

    void  setNodalPositions(vector<myPoint>&  datatemp);

    void  updateNodesAfterPartition();

    void  addBubbleNodes(vector<vector<int> >& elemConn);

    void  computeBasisFunctions1D(double uu, double *N, double *dN_dx);

    void  computeBasisFunctions2D(int deg, double* param, double *N);

    void  computeBasisFunctions2D(bool flag, int eltype, double* param, vector<int>& nodeNum, double *N, double *dN_dx, double *dN_dy, double& Jac);

    void  computeDeformationGradient2D(bool flag, vector<int>& nodeNum, double* dN_dx, double* dN_dy, double* F, double& detF);

    void  computeBasisFunctions3D(int deg, double* param, double *N);

    void  computeBasisFunctions3D(bool flag, int eltype, double* param, vector<int>& nodeNum, double *N, double *dN_dx, double *dN_dy, double* dN_dz, double& Jac);

    void  computeDeformationGradient3D(bool flag, vector<int>& nodeNum, double* dN_dx, double* dN_dy, double* dN_dz, double* F, double& detF);

};








#endif

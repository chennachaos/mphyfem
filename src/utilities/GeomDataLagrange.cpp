
#include "GeomDataLagrange.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "BasisFunctionsLagrange.h"
#include "BasisFunctionsBernstein.h"
#include "QuadratureUtil.h"
#include "ElementBase.h"

extern MyTime myTime;
extern vector<unique_ptr<TimeFunction> > timeFunctions;




GeomDataLagrange::GeomDataLagrange()
{
  bodyForceTimeFunction = -1;

  bodyForce.setZero();
}



GeomDataLagrange::~GeomDataLagrange()
{
}




void GeomDataLagrange::setNodalPositions(vector<myPoint>&  datatemp)
{
    int  ii=0, jj=0, ind=0;
    double  val=0.0;
  
    nNode = datatemp.size();

    NodePosOrig.resize(nNode);
    NodePosCur.resize(nNode);
    NodePosNew.resize(nNode);

    specValCur.resize(nNode);
    specValNew.resize(nNode);

    for(ii=0;ii<nNode;ii++)
    {
      for(jj=0;jj<ndim;jj++)
      {
        //cout << ii << '\t' << jj << '\t' << datatemp[ii][jj] << endl;
        val = datatemp[ii][jj];

        NodePosOrig[ii][jj] = val;
        NodePosCur[ii][jj]  = val;
        NodePosNew[ii][jj]  = val;

        specValCur[ii][jj]  = val;
        specValNew[ii][jj]  = val;
      }
    }

    if(ndim == 2)
    {
      jj=2;
      for(ii=0;ii<nNode;ii++)
      {
        NodePosOrig[ii][jj] = 0.0;
        NodePosCur[ii][jj]  = 0.0;
        NodePosNew[ii][jj]  = 0.0;

        specValCur[ii][jj]  = 0.0;
        specValNew[ii][jj]  = 0.0;
      }
    }

    return;
}



void  GeomDataLagrange::addBubbleNodes(vector<vector<int> >& elemConn)
{
    int nElem = elemConn.size();
    int ee, ii, jj, npElem, ind = nNode+nElem;
    vector<int>  nodeNums;
    myPoint pt;

    NodePosOrig.resize(ind);
    NodePosCur.resize(ind);
    NodePosNew.resize(ind);

    specValCur.resize(ind);
    specValNew.resize(ind);

    for(int ee=0; ee<nElem; ee++)
    {
      npElem = elemConn[ee].size()-5;
      nodeNums.resize(npElem);
      for(ii=0; ii<npElem; ii++)
        nodeNums[ii] = elemConn[ee][5+ii]-1;

      pt.setZero();
      for(jj=0; jj<npElem; jj++)
      {
        pt += NodePosOrig[ nodeNums[jj] ];
      }
      pt /= npElem;

      elemConn[ee].push_back(nNode+ee+1);

      ind = nNode+ee;
      NodePosOrig[ind] = pt;
      NodePosCur[ind]  = pt;
      NodePosNew[ind]  = pt;
    }

    nNode = NodePosCur.size();

    return;
}




void GeomDataLagrange::updateNodesAfterPartition()
{
    myPoint  pt;
    for(int ii=0;ii<nNode;ii++)
    {
        pt = NodePosCur[node_map_get_old[ii]];

        NodePosOrig[ii] = pt;
    }

    for(int ii=0;ii<nNode;ii++)
    {
        pt = NodePosOrig[ii];

        NodePosCur[ii]  = pt;
        NodePosNew[ii]  = pt;
    }

    return;
}




void GeomDataLagrange::computeBasisFunctions1D(double uu, double *N, double *dN_dx)
{
    int  ROWS = 2, ii=0;
/*
    double** ders = new double*[ROWS], du_dx;

    for(ii=0;ii<ROWS;ii++)
       ders[ii] = new double[degree[0]+1];

    //HB_DersBasisFuns(degree[0], start, incr, uu, 1, ders);

    //  dx_du = J;
    //  du_dx = 1.0/dx_du

    //double  du_dx = 1.0/Jacobian[0];
    
    for(ii=0; ii<=degree[0]; ii++)
    {
        N[ii]     =  ders[0][ii];
        dN_dx[ii] =  ders[1][ii] * du_dx;
    }

    for(ii=0;ii<ROWS;ii++)
      delete [] ders[ii];

    delete [] ders;
*/
    return;
}




void GeomDataLagrange::computeBasisFunctions2D(int deg, double* param, double *N)
{
    int  ii=0, jj=0, count = deg + 1;

    vector<double>  N1(count), N2(count);

    Lagrange_BasisFuns1D(deg, param[0], &N1[0]);
    Lagrange_BasisFuns1D(deg, param[1], &N2[0]);

    count = 0;
    for(jj=0; jj<=deg; jj++)
    {
      for(ii=0; ii<=deg; ii++)
        N[count++]   =  N2[jj] *  N1[ii];
    }

    return;
}




void  GeomDataLagrange::computeBasisFunctions2D(bool flag, int eltype, double* param, vector<int>& nodeNums, double *N, double *dN_dx, double *dN_dy, double& Jac)
{
    int  npElem = nodeNums.size();

    double  dN_du1[npElem], dN_du2[npElem];

    if(eltype == ELEM_SHAPE_TRIA) // triangular elements
    {
      LagrangeBasisFunsTria(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else if(eltype == ELEM_SHAPE_TRIA_BERNSTEIN) // triangular elements
    {
      BernsteinBasisFunsTria(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else  if(eltype == ELEM_SHAPE_QUAD) // Lagrange quad elements
    {
      LagrangeBasisFunsQuad(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else  if(eltype == ELEM_SHAPE_QUAD_BERNSTEIN) // Bernstein quad elements
    {
      BernsteinBasisFunsQuad(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else
    {
      throw runtime_error("GeomDataLagrange::computeBasisFunctions2D ... Invalid element type");
    }


    // Gradient of mapping from parameter space to physical space
    double  xx, yy;
    MatrixXd  B(2,2), Binv(2,2);
    B.setZero();

    if(flag)
    {
       for(int ii=0; ii<npElem; ii++)
       {
          xx = NodePosCur[nodeNums[ii]][0];
          yy = NodePosCur[nodeNums[ii]][1];

          B(0,0) +=  (xx * dN_du1[ii]) ;
          B(1,0) +=  (xx * dN_du2[ii]) ;
          B(0,1) +=  (yy * dN_du1[ii]) ;
          B(1,1) +=  (yy * dN_du2[ii]) ;
       }
    }
    else
    {
       for(int ii=0; ii<npElem; ii++)
       {
          xx = NodePosOrig[nodeNums[ii]][0];
          yy = NodePosOrig[nodeNums[ii]][1];

          B(0,0) +=  (xx * dN_du1[ii]) ;
          B(1,0) +=  (xx * dN_du2[ii]) ;
          B(0,1) +=  (yy * dN_du1[ii]) ;
          B(1,1) +=  (yy * dN_du2[ii]) ;
       }
    }

    Jac  = B.determinant();
    Binv = B.inverse();

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(int ii=0; ii<npElem; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv(0,0) + dN_du2[ii] * Binv(0,1);
      dN_dy[ii] = dN_du1[ii] * Binv(1,0) + dN_du2[ii] * Binv(1,1);
    }

    return;
}



void  GeomDataLagrange::computeDeformationGradient2D(bool flag, vector<int>& nodeNums, double* dN_dx, double* dN_dy, double* F, double& detF)
{
    /////////////////////////////////////////////
    //                                         //
    //    F(1,1) = F[0]   //   F(1,2) = F[1]   //
    //                                         //
    //    F(2,1) = F[2]   //   F(2,2) = F[3]   //
    //                                         //
    /////////////////////////////////////////////

    //   confflag = 1 --> coordiantes from current configuration
    //                    shape function derivatives w.r.t reference configuration

    //   confflag = 2 ==> coordiantes from reference configuration
    //                    shape function derivatives w.r.t current configuration

    double  xx=0.0, yy=0.0;

    F[0] = F[1] = F[2] = F[3] = 0.0;

    for(int ii=0; ii<nodeNums.size(); ii++)
    {
      //xx = NodePosCur[node_map_new_to_old[nodeNums[ii]]][0];
      //yy = NodePosCur[node_map_new_to_old[nodeNums[ii]]][1];

      xx = NodePosCur[nodeNums[ii]][0];
      yy = NodePosCur[nodeNums[ii]][1];

      F[0] += xx * dN_dx[ii];
      F[2] += xx * dN_dy[ii];
      F[1] += yy * dN_dx[ii];
      F[3] += yy * dN_dy[ii];
    }

    detF = F[0]*F[3] - F[1]*F[2];

    return;
}





void GeomDataLagrange::computeBasisFunctions3D(int deg, double* param, double *N)
{
    int  ii=0, jj=0, kk=0, count = deg + 1;

    vector<double>  N1(count), N2(count), N3(count);

    Lagrange_BasisFuns1D(deg, param[0], &N1[0]);
    Lagrange_BasisFuns1D(deg, param[1], &N2[0]);
    Lagrange_BasisFuns1D(deg, param[2], &N3[0]);

    count = 0;
    for(kk=0; kk<=deg; kk++)
    {
      for(jj=0; jj<=deg; jj++)
      {
        for(ii=0; ii<=deg; ii++)
          N[count++]   =  N3[kk] * N2[jj] * N1[ii];
      }
    }

    return;
}





void  GeomDataLagrange::computeBasisFunctions3D(bool flag, int eltype, double* param, vector<int>& nodeNums, double *N, double *dN_dx, double *dN_dy, double *dN_dz, double& Jac)
{
    int  npElem = nodeNums.size();

    double  dN_du1[npElem], dN_du2[npElem], dN_du3[npElem];

    if(eltype == ELEM_SHAPE_TETRA) // Lagrange tetrahedral elements
    {
      LagrangeBasisFunsTetra(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(eltype == ELEM_SHAPE_TETRA_BERNSTEIN) // Bernstein tetrahedral elements
    {
      BernsteinBasisFunsTetra(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(eltype == ELEM_SHAPE_PYRAMID) // pyramid elements
    {
      LagrangeBasisFunsPyramid(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(eltype == ELEM_SHAPE_WEDGE_BERNSTEIN) // wedge elements
    {
      BernsteinBasisFunsWedge(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(eltype == ELEM_SHAPE_HEXA) // hex elements
    {
      LagrangeBasisFunsHexa(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(eltype == ELEM_SHAPE_HEXA_BERNSTEIN) // hex elements
    {
      BernsteinBasisFunsHexa(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else
    {
      throw runtime_error("GeomDataLagrange::computeBasisFunctions3D ... Invalid element type");
    }


    double  xx=0.0, yy=0.0, zz=0.0;
    MatrixXd  B(3,3), Binv(3,3);

    // Gradient of mapping from parameter space to physical space
    B.setZero();

    if(flag)
    {
       for(int ii=0; ii<npElem; ii++)
       {
          //xx = NodePosCur[node_map_new_to_old[nodeNums[ii]]][0];
          //yy = NodePosCur[node_map_new_to_old[nodeNums[ii]]][1];
          //zz = NodePosCur[node_map_new_to_old[nodeNums[ii]]][2];

          xx = NodePosCur[nodeNums[ii]][0];
          yy = NodePosCur[nodeNums[ii]][1];
          zz = NodePosCur[nodeNums[ii]][2];

          B(0,0) +=  (xx * dN_du1[ii]) ;
          B(1,0) +=  (xx * dN_du2[ii]) ;
          B(2,0) +=  (xx * dN_du3[ii]) ;

          B(0,1) +=  (yy * dN_du1[ii]) ;
          B(1,1) +=  (yy * dN_du2[ii]) ;
          B(2,1) +=  (yy * dN_du3[ii]) ;

          B(0,2) +=  (zz * dN_du1[ii]) ;
          B(1,2) +=  (zz * dN_du2[ii]) ;
          B(2,2) +=  (zz * dN_du3[ii]) ;
      }
    }
    else
    {
       for(int ii=0; ii<npElem; ii++)
      {
          //xx = NodePosOrig[node_map_new_to_old[nodeNums[ii]]][0];
          //yy = NodePosOrig[node_map_new_to_old[nodeNums[ii]]][1];
          //zz = NodePosOrig[node_map_new_to_old[nodeNums[ii]]][2];

          xx = NodePosOrig[nodeNums[ii]][0];
          yy = NodePosOrig[nodeNums[ii]][1];
          zz = NodePosOrig[nodeNums[ii]][2];

          B(0,0) +=  (xx * dN_du1[ii]) ;
          B(1,0) +=  (xx * dN_du2[ii]) ;
          B(2,0) +=  (xx * dN_du3[ii]) ;

          B(0,1) +=  (yy * dN_du1[ii]) ;
          B(1,1) +=  (yy * dN_du2[ii]) ;
          B(2,1) +=  (yy * dN_du3[ii]) ;

          B(0,2) +=  (zz * dN_du1[ii]) ;
          B(1,2) +=  (zz * dN_du2[ii]) ;
          B(2,2) +=  (zz * dN_du3[ii]) ;
      }
    }

    //printMatrix(B);

    Jac  = B.determinant();
    Binv = B.inverse();

    //printMatrix(Binv);

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(int ii=0; ii<npElem; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv(0,0) + dN_du2[ii] * Binv(0,1) + dN_du3[ii] * Binv(0,2);
      dN_dy[ii] = dN_du1[ii] * Binv(1,0) + dN_du2[ii] * Binv(1,1) + dN_du3[ii] * Binv(1,2);
      dN_dz[ii] = dN_du1[ii] * Binv(2,0) + dN_du2[ii] * Binv(2,1) + dN_du3[ii] * Binv(2,2);
    }

    return;
}




void  GeomDataLagrange::computeDeformationGradient3D(bool flag, vector<int>& nodeNums, double* dN_dx, double* dN_dy, double* dN_dz, double* F, double& detF)
{
    double  xx=0.0, yy=0.0, zz=0.0;

    F[0] = F[1] = F[2] = 0.0;
    F[3] = F[4] = F[5] = 0.0;
    F[6] = F[7] = F[8] = 0.0;

    for(int ii=0; ii<nodeNums.size(); ii++)
    {
        xx = NodePosCur[node_map_get_old[nodeNums[ii]]][0];
        yy = NodePosCur[node_map_get_old[nodeNums[ii]]][1];
        zz = NodePosCur[node_map_get_old[nodeNums[ii]]][2];

        //cout << xx << '\t' << yy << endl;

        F[0] += xx * dN_dx[ii];
        F[3] += xx * dN_dy[ii];
        F[6] += xx * dN_dz[ii];

        F[1] += yy * dN_dx[ii];
        F[4] += yy * dN_dy[ii];
        F[7] += yy * dN_dz[ii];

        F[2] += zz * dN_dx[ii];
        F[5] += zz * dN_dy[ii];
        F[8] += zz * dN_dz[ii];
    }

    detF = F[0]*(F[4]*F[8] - F[5]*F[7]) - F[3]*(F[1]*F[8] - F[2]*F[7]) + F[6]*(F[1]*F[5] - F[2]*F[4]);

    return;
}





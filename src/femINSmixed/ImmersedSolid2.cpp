#include "ImmersedSolid.h"
#include "SolutionData.h"
#include "GeomDataLagrange.h"
#include "ImmersedIntegrationElement.h"
#include "util.h"
#include <boost/algorithm/string.hpp>


using namespace myGeom;



void ImmersedSolid::postProcess(string dirname, int fileCount)
{
    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();

    //vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    //vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    //vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK     =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    //vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK     =  vtkSmartPointer<vtkBiQuadraticQuad>::New();

    vtkSmartPointer<vtkFloatArray>           vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2;

    vtkIdType pt[10];

    if(ndim == 2)
    {
      for(ii=0; ii<nNode; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(nodePosCur[ii][0], nodePosCur[ii][1], 0.0);
      }

      for(ee=0; ee<nElem; ee++)
      {
        if( npElem == 1 )
        {
          vertexVTK->GetPointIds()->SetId(0, ImmIntgElems[ee]->nodeNums[0] );

          uGridVTK->InsertNextCell(vertexVTK->GetCellType(), vertexVTK->GetPointIds());
        }
        if( npElem == 2 )
        {
          for(ii=0; ii<npElem; ii++)
            lineVTK->GetPointIds()->SetId(ii, ImmIntgElems[ee]->nodeNums[ii] );

          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
        }
      }
    }
    else // (ndim == 3)
    {
    }

    uGridVTK->SetPoints(pointsVTK);

    //Write the file.

    char VTKfilename[200];

    sprintf(VTKfilename,"%s%s%02d%s%06d%s", dirname.c_str(), "-IB", id, "-", fileCount,".vtu");

    writerUGridVTK->SetFileName(VTKfilename);

    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    return;
}




void ImmersedSolid::adjustBoundaryPoints(double* minVal, double* maxVal)
{
/*
    //cout << " ImmersedSolid::adjustBoundaryPoints .... needs to be modified " << endl;
//#ifdef CUTCELL3D_ALGO_VTK
    int ii=0, jj=0;
    double  tol1=1.0e-8, tol2=1.0e-4;

    myPoint  pt1, pt2;

    for(ii=0; ii<GeomData.NodePosOrig.size(); ii++)
    {
      pt1 = GeomData.NodePosOrig[ii];

      for(jj=0; jj<DIM; jj++)
      {
        if( abs(pt1[jj]-minVal[jj]) < tol1 )
        {
          GeomData.NodePosOrig[ii][jj] -= tol2;
          GeomData.NodePosCur[ii][jj]  -= tol2;
        }

        if( abs(pt1[jj]-maxVal[jj]) < tol1 )
        {
          GeomData.NodePosOrig[ii][jj] += tol2;
          GeomData.NodePosCur[ii][jj]  += tol2;
        }
      }
    }
//#endif
*/
    return;
}





int ImmersedSolid::within(myPoint& ptTemp)
{
/*
  // check if the point is outside/inside the boundingbox of the polygon
  // if the point is outside the boundingbox then
  // the point is outside the polygon
  // otherwise, the point may be inside the polygon and needs further checking

  if( !(bbox.within(ptTemp)) )
    return 0;

  int i=0, c=0;

  if(DIM == 2)
  {
    ray1.updateOrigin(ptTemp);

    vector<myPoint>  ptOut;

    c=0;
    for(i=0; i<ImmersedFaces.size(); i++)
    {
      if( ImmersedFaces[i]->IntersectWithRay(ray1, ptOut ) )
        c = !c;
    }
  }
  else
  {
#ifdef CUTCELL3D_ALGO_VTK
    c =  selectEnclosedPoints->IsInsideSurface( ptTemp[0], ptTemp[1], ptTemp[2] );
#else
    c = checkBoundedSideCGAL((*point_inside_tester)(CGAL_Point(ptTemp[0], ptTemp[1], ptTemp[2])));
#endif
  }
  return c;
*/
  return 0;
}


int  ImmersedSolid::doAABBintersect(AABB&  bb2)
{
  if( bbox.doIntersect(bb2) )
    return 1;
  else
    return 0;
}



int  ImmersedSolid::doIntersect2D(AABB&  bbTemp, bool flag, vector<int>& vecCorners, vector<myPoint>&  ptOut)
{
    return 0;
}




int  ImmersedSolid::doIntersect2Dfor3D(int sideTemp, double coord3, AABB&  bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>&  ptOut)
{
    return 0;
}




int  ImmersedSolid::doIntersect3D(AABB&  bbTemp, bool flag, vector<int>& vecTemp, vector<myPoint>&  ptOut)
{
    return 0;
}



double  ImmersedSolid::distanceFromPoint(myPoint&  pt)
{
  return 0.0;
}


double  ImmersedSolid::distanceFromPoint(double xx, double yy, double zz)
{
  return 0.0;
}




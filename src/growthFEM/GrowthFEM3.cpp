
#include "GrowthFEM.h"
#include "MyTime.h"
#include "Files.h"
#include "MyString.h"
#include "ElementBase.h"

//#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>


extern MyTime myTime;
extern Files files;



void GrowthFEM::plotGeom()
{
    int  dd, ii, jj, kk, ll, nlocal, index, ind1, ind2, e, ee, count, gcount, ind;
    double xx, yy, zz;
    //                        {0,1,2,3,4,5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    int  hexa27nodemap[27] =  {0,1,2,3,4,5,6,7,8,11,16, 9,17,10,18,19,12,15,13,14,24,22,20,21,23,25,26};
    //                        {0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15,16,17};
    int  wedge18nodemap[18] = {0,1,2,3,4,5,6,8,12,7,13,14, 9,11,10,15,17,16};
    vtkIdType pt[8];

    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();
    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();

    vtkSmartPointer<vtkPolygon>              polygonVTK   =  vtkSmartPointer<vtkPolygon>::New();

    vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkHexahedron>           hexaVTK      =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK     =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK     =  vtkSmartPointer<vtkBiQuadraticQuad>::New();
    vtkSmartPointer<vtkQuadraticTetra>       tetra10VTK   =  vtkSmartPointer<vtkQuadraticTetra>::New();
    vtkSmartPointer<vtkBiQuadraticQuadraticWedge>  wedge18VTK = vtkSmartPointer<vtkBiQuadraticQuadraticWedge>::New();
    vtkSmartPointer<vtkTriQuadraticHexahedron>  hexa27VTK    =  vtkSmartPointer<vtkTriQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           cellDataVTKmatltype  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKelemtype  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKBappl  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKBresi  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    count = GeomData.NodePosOrig.size();

    if(ndim == 2)
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        if( midnodeData[ii][0] ) // if midnode
        {
          xx  = 0.50*GeomData.NodePosOrig[ii][0];
          xx += 0.25*GeomData.NodePosOrig[midnodeData[ii][1]][0];
          xx += 0.25*GeomData.NodePosOrig[midnodeData[ii][2]][0];

          yy  = 0.50*GeomData.NodePosOrig[ii][1];
          yy += 0.25*GeomData.NodePosOrig[midnodeData[ii][1]][1];
          yy += 0.25*GeomData.NodePosOrig[midnodeData[ii][2]][1];
        }
        else
        {
          xx = GeomData.NodePosOrig[ii][0];
          yy = GeomData.NodePosOrig[ii][1];
        }

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);
      }

      if( elems[0]->nodeNums.size() == 2 ) // line
      {
        for(ee=0;ee<nElem;ee++)
        {
          lineVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          lineVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);

          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
        }
      }
      else if( elems[0]->nodeNums.size() == 3 )             // tria
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
      }
      else if( elems[0]->nodeNums.size() == 6 ) // tria, 6 node
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<6;ll++)
            tria6VTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
      }
      else if( elems[0]->nodeNums.size() == 4 ) // quad, 4 node
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<4;ll++)
            quadVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
      }
      else if( elems[0]->nodeNums.size() == 9 ) // quad, 9 node
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<9;ll++)
            quad9VTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
      }
    }
    else
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        if( midnodeData[ii][0] ) // if midnode
        {
          xx  = 0.50*GeomData.NodePosOrig[ii][0];
          xx += 0.25*GeomData.NodePosOrig[midnodeData[ii][1]][0];
          xx += 0.25*GeomData.NodePosOrig[midnodeData[ii][2]][0];

          yy  = 0.50*GeomData.NodePosOrig[ii][1];
          yy += 0.25*GeomData.NodePosOrig[midnodeData[ii][1]][1];
          yy += 0.25*GeomData.NodePosOrig[midnodeData[ii][2]][1];

          zz  = 0.50*GeomData.NodePosOrig[ii][2];
          zz += 0.25*GeomData.NodePosOrig[midnodeData[ii][1]][2];
          zz += 0.25*GeomData.NodePosOrig[midnodeData[ii][2]][2];
        }
        else
        {
          xx = GeomData.NodePosOrig[ii][0];
          yy = GeomData.NodePosOrig[ii][1];
          zz = GeomData.NodePosOrig[ii][2];
        }

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);
      }

      if(elems[0]->nodeNums.size() == 4) // tet
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<4;ll++)
            tetraVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 8) // hexa
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<8;ll++)
            hexaVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(hexaVTK->GetCellType(), hexaVTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 10) // 10-node tetrahedron
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ii=0; ii<10; ii++)
            tetra10VTK->GetPointIds()->SetId(ii, elems[ee]->nodeNums[ii]);

          uGridVTK->InsertNextCell(tetra10VTK->GetCellType(), tetra10VTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 18) // 18-node wedge
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ii=0; ii<18; ii++)
            wedge18VTK->GetPointIds()->SetId(wedge18nodemap[ii], elems[ee]->nodeNums[ii]);
          uGridVTK->InsertNextCell(wedge18VTK->GetCellType(), wedge18VTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 27) // 27-node hexahedron
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ii=0; ii<27; ii++)
            hexa27VTK->GetPointIds()->SetId(hexa27nodemap[ii], elems[ee]->nodeNums[ii]);

          uGridVTK->InsertNextCell(hexa27VTK->GetCellType(), hexa27VTK->GetPointIds());
        }
      }
      else
      {
        cout << "Wrong element type for plotGeom "  << endl;
        exit(1);
      }
    }

    uGridVTK->SetPoints(pointsVTK);

    // create a write object and write uGridVTK to it

    char fname[200];
    //sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-Geom.vtu");
    sprintf(fname,"%s%s%s%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(), "-Geom.vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    return;
}




void  GrowthFEM::postProcess()
{
    if(DEBUG) {cout << "     GrowthFEM::postProcess ... STARTED \n\n";}

    bool extrapolateFlag = false;
    int  vartype = 1;
    int  varindex = 4;
    int  index = 10;

    int  dd, ii, jj, kk, ll, nlocal, ind1, ind2, e, ee, count, gcount, ind;
    int  n1, n2, n3, n4, n5, n6;
    double vec[3], vec2[3], geom[3], xx, yy, zz, val, maxdisp = 0.0, maxstress = 0.0;
    //                        {0,1,2,3,4,5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    int  hexa27nodemap[27] =  {0,1,2,3,4,5,6,7,8,11,16, 9,17,10,18,19,12,15,13,14,24,22,20,21,23,25,26};
    //                        {0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15,16,17};
    int  wedge18nodemap[18] = {0,1,2,3,4,5,6,8,12,7,13,14, 9,11,10,15,17,16};

    vtkIdType pt[27];
    vector<int>  nodeNumsElem;

    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();
    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();

    vtkSmartPointer<vtkPolygon>              polygonVTK   =  vtkSmartPointer<vtkPolygon>::New();

    vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkWedge>                wedgeVTK     =  vtkSmartPointer<vtkWedge>::New();
    vtkSmartPointer<vtkHexahedron>           hexaVTK      =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK     =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK     =  vtkSmartPointer<vtkBiQuadraticQuad>::New();
    vtkSmartPointer<vtkQuadraticTetra>       tetra10VTK   =  vtkSmartPointer<vtkQuadraticTetra>::New();
    vtkSmartPointer<vtkBiQuadraticQuadraticWedge>  wedge18VTK = vtkSmartPointer<vtkBiQuadraticQuadraticWedge>::New();
    vtkSmartPointer<vtkTriQuadraticHexahedron>  hexa27VTK    =  vtkSmartPointer<vtkTriQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           dispVTK       = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           veloVTK       = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           SxxNodesVTK     =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           SyyNodesVTK     =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           SzzNodesVTK     =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           mpotVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTK   =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKmatltype  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKelemtype  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKBappl  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKBresi  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    vtkSmartPointer<vtkFloatArray>           vecVTK    =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           vecVTK2   =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK    =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK2   =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           strainVTK =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           stressVTK =  vtkSmartPointer<vtkFloatArray>::New();


    //cout << " jjjjjjjjjjjjjjjjjj " << endl;

    count = GeomData.NodePosOrig.size();

    vecVTK->SetName("disp");
    vecVTK2->SetName("velocity");
    //if(ndof > 1)
    //{
      vecVTK->SetNumberOfComponents(3);
      vecVTK->SetNumberOfTuples(count);
      vecVTK2->SetNumberOfComponents(3);
      vecVTK2->SetNumberOfTuples(count);
    //}

    //cout << " BBBBBBBBBBBB " << endl;

    int idd = SolnData.ElemProp[0]->id;


    if(ndim == 2)
    {
      vec[2] = 0.0;
      for(ii=0;ii<nNode;ii++)
      {
        if( midnodeData[ii][0] ) // if midnode
        {
          xx  = 0.5*GeomData.NodePosNew[ii][0];
          xx += 0.25*GeomData.NodePosNew[midnodeData[ii][1]][0];
          xx += 0.25*GeomData.NodePosNew[midnodeData[ii][2]][0];

          yy  = 0.5*GeomData.NodePosNew[ii][1];
          yy += 0.25*GeomData.NodePosNew[midnodeData[ii][1]][1];
          yy += 0.25*GeomData.NodePosNew[midnodeData[ii][2]][1];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          n1 = ii*ndof;
          n2 = midnodeData[ii][1]*ndof;
          n3 = midnodeData[ii][2]*ndof;

          vec[0]  = 0.50*SolnData.var1[n1];
          vec[0] += 0.25*SolnData.var1[n2];
          vec[0] += 0.25*SolnData.var1[n3];

          vec[1]  = 0.50*SolnData.var1[n1+1];
          vec[1] += 0.25*SolnData.var1[n2+1];
          vec[1] += 0.25*SolnData.var1[n3+1];

          vecVTK->InsertTuple(ii, vec);

          maxdisp = min(maxdisp, vec[1]);

          vec[0]  = 0.50*SolnData.var1Dot[n1];
          vec[0] += 0.25*SolnData.var1Dot[n2];
          vec[0] += 0.25*SolnData.var1Dot[n3];

          vec[1]  = 0.50*SolnData.var1Dot[n1+1];
          vec[1] += 0.25*SolnData.var1Dot[n2+1];
          vec[1] += 0.25*SolnData.var1Dot[n3+1];

          vecVTK2->InsertTuple(ii, vec);
        }
        else
        {
          xx = GeomData.NodePosNew[ii][0];
          yy = GeomData.NodePosNew[ii][1];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          kk = ii*ndof;

          vec[0] = SolnData.var1[kk];
          vec[1] = SolnData.var1[kk+1];

          vecVTK->InsertTuple(ii, vec);

          vec[0] = SolnData.var1Dot[kk];
          vec[1] = SolnData.var1Dot[kk+1];

          vecVTK2->InsertTuple(ii, vec);
        }
      }

      // center nodes for 9-noded quadrilateral elements
      if( idd == 6003 )
      {
        for(ee=0;ee<nElem;ee++)
        {
          xx  = GeomData.NodePosNew[elems[ee]->nodeNums[8]][0];
          yy  = GeomData.NodePosNew[elems[ee]->nodeNums[8]][1];

          pointsVTK->SetPoint(elems[ee]->nodeNums[8], xx, yy, 0.0);

          vec[0]  = SolnData.var1[ndof*elems[ee]->nodeNums[8]];
          vec[1]  = SolnData.var1[ndof*elems[ee]->nodeNums[8]+1];

          vecVTK->SetTuple(elems[ee]->nodeNums[8], vec);

          vec[0]  = SolnData.var1Dot[ndof*elems[ee]->nodeNums[8]];
          vec[1]  = SolnData.var1Dot[ndof*elems[ee]->nodeNums[8]+1];

          vecVTK2->SetTuple(elems[ee]->nodeNums[8], vec);
        }
      }

      if( elems[0]->nodeNums.size() == 2 ) // line
      {
        for(ee=0;ee<nElem;ee++)
        {
          lineVTK->GetPointIds()->SetId(0, elems[ee]->nodeNums[0]);
          lineVTK->GetPointIds()->SetId(1, elems[ee]->nodeNums[1]);

          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 3)
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
      }
      else if( elems[0]->nodeNums.size() == 6 ) // tria, 6 node
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<6;ll++)
            tria6VTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
      }
      else if( elems[0]->nodeNums.size() == 4 ) // quad, 4 node
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<4;ll++)
            quadVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
      }
      else if( elems[0]->nodeNums.size() == 9 ) // quad, 9 node
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<9;ll++)
            quad9VTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
      }
    }
    else //if(ndim == 3)
    {
      for(ii=0;ii<nNode;ii++)
      {
        if( midnodeData[ii][0] ) // if midnode
        {
          xx  = 0.50*GeomData.NodePosNew[ii][0];
          xx += 0.25*GeomData.NodePosNew[midnodeData[ii][1]][0];
          xx += 0.25*GeomData.NodePosNew[midnodeData[ii][2]][0];

          yy  = 0.50*GeomData.NodePosNew[ii][1];
          yy += 0.25*GeomData.NodePosNew[midnodeData[ii][1]][1];
          yy += 0.25*GeomData.NodePosNew[midnodeData[ii][2]][1];

          zz  = 0.50*GeomData.NodePosNew[ii][2];
          zz += 0.25*GeomData.NodePosNew[midnodeData[ii][1]][2];
          zz += 0.25*GeomData.NodePosNew[midnodeData[ii][2]][2];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          n1 = ii*ndof;
          n2 = midnodeData[ii][1]*ndof;
          n3 = midnodeData[ii][2]*ndof;

          vec[0]  = 0.50*SolnData.var1[n1];
          vec[0] += 0.25*SolnData.var1[n2];
          vec[0] += 0.25*SolnData.var1[n3];

          vec[1]  = 0.50*SolnData.var1[n1+1];
          vec[1] += 0.25*SolnData.var1[n2+1];
          vec[1] += 0.25*SolnData.var1[n3+1];

          vec[2]  = 0.50*SolnData.var1[n1+2];
          vec[2] += 0.25*SolnData.var1[n2+2];
          vec[2] += 0.25*SolnData.var1[n3+2];

          vecVTK->InsertTuple(ii, vec);

          vec[0]  = 0.50*SolnData.var1Dot[n1];
          vec[0] += 0.25*SolnData.var1Dot[n2];
          vec[0] += 0.25*SolnData.var1Dot[n3];

          vec[1]  = 0.50*SolnData.var1Dot[n1+1];
          vec[1] += 0.25*SolnData.var1Dot[n2+1];
          vec[1] += 0.25*SolnData.var1Dot[n3+1];

          vec[2]  = 0.50*SolnData.var1Dot[n1+2];
          vec[2] += 0.25*SolnData.var1Dot[n2+2];
          vec[2] += 0.25*SolnData.var1Dot[n3+2];

          vecVTK2->InsertTuple(ii, vec);
        }
        else
        {
          xx = GeomData.NodePosOrig[ii][0];
          yy = GeomData.NodePosOrig[ii][1];
          zz = GeomData.NodePosOrig[ii][2];

          //xx = GeomData.NodePosNew[ii][0];
          //yy = GeomData.NodePosNew[ii][1];
          //zz = GeomData.NodePosNew[ii][2];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          kk = ii*ndof;

          vec[0] = SolnData.var1[kk];
          vec[1] = SolnData.var1[kk+1];
          vec[2] = SolnData.var1[kk+2];

          vecVTK->InsertTuple(ii, vec);

          vec[0] = SolnData.var1Dot[kk];
          vec[1] = SolnData.var1Dot[kk+1];
          vec[2] = SolnData.var1Dot[kk+2];

          vecVTK2->InsertTuple(ii, vec);
        }
      }

        if( idd == 6053 ) // Wedge element
        {
          double  params[3][3] = {{0.5, 0.0, 0.0}, {0.0, 0.5, 0.0}, {0.5, 0.5, 0.0}};
          int  facenodes[6] = {15, 16, 17};

          for(ee=0;ee<nElem;ee++)
          {
            nodeNumsElem = elems[ee]->nodeNums;

            for(int ii=0; ii<3; ii++)
            {
              elems[ee]->computeGeomNew(params[ii], geom);

              pointsVTK->SetPoint(nodeNumsElem[facenodes[ii]], geom[0], geom[1], geom[2]);
              //cout << nodeNumsTemp[0] << '\t' << geom[0] << '\t' << geom[1] << '\t' << geom[2] << endl;
            }
          }
        }


        if( idd == 6054 ) // Hexa element
        {
          double  params[6][3] = {{0.0, 0.0, -1.0}, {0.0, -1.0, 0.0}, {-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
          int  facenodes[6] = {20, 21, 22, 23, 24, 25};

          for(ee=0;ee<nElem;ee++)
          {
            nodeNumsElem = elems[ee]->nodeNums;

            for(int ii=0; ii<6; ii++)
            {
              elems[ee]->computeGeomNew(params[ii], geom);

              pointsVTK->SetPoint(nodeNumsElem[facenodes[ii]], geom[0], geom[1], geom[2]);
              //cout << nodeNumsTemp[0] << '\t' << geom[0] << '\t' << geom[1] << '\t' << geom[2] << endl;
            }
          }
        }

      if(elems[0]->nodeNums.size() == 4) // tet
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<4;ll++)
            tetraVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 8) // hexa
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ll=0;ll<8;ll++)
            hexaVTK->GetPointIds()->SetId(ll, elems[ee]->nodeNums[ll]);

          uGridVTK->InsertNextCell(hexaVTK->GetCellType(), hexaVTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 10) // 10-node tetrahedron
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ii=0; ii<10; ii++)
            tetra10VTK->GetPointIds()->SetId(ii, elems[ee]->nodeNums[ii]);

          uGridVTK->InsertNextCell(tetra10VTK->GetCellType(), tetra10VTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 18) // 18-node wedge
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ii=0; ii<18; ii++)
            wedge18VTK->GetPointIds()->SetId(wedge18nodemap[ii], elems[ee]->nodeNums[ii]);

          uGridVTK->InsertNextCell(wedge18VTK->GetCellType(), wedge18VTK->GetPointIds());
        }
      }
      else if(elems[0]->nodeNums.size() == 27) // 27-node hexahedron
      {
        for(ee=0;ee<nElem;ee++)
        {
          for(ii=0; ii<27; ii++)
            hexa27VTK->GetPointIds()->SetId(hexa27nodemap[ii], elems[ee]->nodeNums[ii]);

          uGridVTK->InsertNextCell(hexa27VTK->GetCellType(), hexa27VTK->GetPointIds());
        }
      }
      else
      {
        cout << "Wrong element type for post-processing "  << endl;
        exit(1);
      }
    }
    
    cout << " iiiiiiiiiiiiii " << endl;

    cellDataVTK->SetName("pres");
    cellDataVTK->SetNumberOfTuples(nElem);
    cellDataVTKmatltype->SetName("matlType");
    cellDataVTKmatltype->SetNumberOfTuples(nElem);
    cellDataVTKelemtype->SetName("elemType");
    cellDataVTKelemtype->SetNumberOfTuples(nElem);

    for(ee=0;ee<nElem;ee++)
    {
      elems[ee]->elementContourplot(vartype, varindex, index);

      cellDataVTK->InsertTuple1(ee, elems[ee]->vals2project[0]);

      cellDataVTKmatltype->InsertTuple1(ee, elems[ee]->getMaterialTypeNumber() );

      cellDataVTKelemtype->InsertTuple1(ee, elems[ee]->getElementTypeNumber() );
    }

    if(intVarFlag)
    {
      scaVTK2->SetName("epstrn");
      scaVTK2->SetNumberOfTuples(count);
      projectStresses(1, 5, 6, index);
    }

    scaVTK->SetName("pres");
    scaVTK->SetNumberOfTuples(count);

    if(MIXED_ELEMENT)
    {
        if( (idd == 6003) || (idd == 6004) || (idd == 6053) || (idd == 6054) || (idd == 6055) )
        {
          for(ii=0;ii<nNode;ii++)
          {
            val = SolnData.var2[ii];

            if( midnodeData[ii][0] )
            {
              val  = 0.5*SolnData.var2[midnodeData[ii][1]];
              val += 0.5*SolnData.var2[midnodeData[ii][2]];
            }
            scaVTK->SetTuple1(ii, val);
          }
        }

        // face center node for the quadratic quad elements
        if( idd == 6003 )
        {
          for(ee=0;ee<nElem;ee++)
          {
            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[8], val);
          }
        }
        else if( idd == 6053 ) // Wedge element
        {
          for(ee=0;ee<nElem;ee++)
          {
            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[15], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[16], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[17], val);
          }
        }
        else if( idd == 6054 ) // Hexa element
        {
          for(ee=0;ee<nElem;ee++)
          {
            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[20], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[21], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[7]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[22], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[23], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[7]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[24], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[7]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[25], val);

            val  = 0.125*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[5]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[7]];

            scaVTK->SetTuple1(elems[ee]->nodeNums[26], val);
          }
        }
    }
    else
    {
        ///cout << vartype << '\t' << varindex << '\t' << index << endl;
        //projectStresses(extrapolateFlag, vartype, varindex, index);
        //cout << vartype << '\t' << varindex << '\t' << index << endl;
        //projectFromElemsToNodes(extrapolateFlag, vartype, varindex, index);
    }

    cout << " iiiiiiiiiiiiii " << endl;

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);
    uGridVTK->GetPointData()->AddArray(vecVTK2);
    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTKmatltype);
    uGridVTK->GetCellData()->AddArray(cellDataVTKelemtype);
    if(intVarFlag)
      uGridVTK->GetPointData()->AddArray(scaVTK2);

    // create a write object and write uGridVTK to it

    char fname[200];
    //sprintf(fname,"%s%s%06d%s",     files.Ofile.asCharArray(),"-",filecount, ".vtu");
    sprintf(fname,"%s%s%s%s%06d%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(),"-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    if(DEBUG) {cout << "     GrowthFEM::postProcess ... FINISHED \n\n";}

    return;
}



void  GrowthFEM::projectStrains(bool extrapolateFlag, int vartype, int varindex, int index)
{
    return;
}



void  GrowthFEM::projectStresses(bool extrapolateFlag, int vartype, int varindex, int index)
{
    int  ee, ii;

    //cout << " extrapolateFlag = " << extrapolateFlag << endl;

    // compute the stresses at element nodes
    for(ee=0; ee<nElem; ee++)
    {
      elems[ee]->projectToNodes(extrapolateFlag, vartype, varindex, index);
    }

    vector<int>  nodeNums;
    vector<double>  output(nNode, 0.0);
    double  volu, val;

    for(ee=0; ee<nElem; ee++)
    {
      nodeNums = elems[ee]->nodeNums;

      elems[ee]->computeVolume(false);

      volu = elems[ee]->getVolume();

      //cout << " volu= " << volu << endl;

      for(ii=0; ii<nodeNums.size(); ii++)
      {
        //output[nodeNums[ii]] += volu*elems[ee]->vals2project[ii];
        output[nodeNums[ii]] += elems[ee]->vals2project[ii];
      }
    }

    for(ee=0; ee<nNode; ee++)
    {
      //printVector(node_elem_conn[ee]);
      volu = 0.0;
      for(ii=0; ii<node_elem_conn[ee].size(); ii++)
      {
        volu += elems[node_elem_conn[ee][ii]]->getVolume();
      }

      //output[ee] /= volu;
      output[ee] /= node_elem_conn[ee].size();
    }

    // for Bernstein elements, the values for midnodes need to be interpolated
    // Only to be done when the values are extrapolated from Gauss points
    //if(extrapolateFlag)
    //{
      for(ii=0; ii<nNode; ii++)
      {
        if( midnodeData[ii][0] )
        {
          val  = 0.50*output[ii];
          val += 0.25*output[midnodeData[ii][1]];
          val += 0.25*output[midnodeData[ii][2]];

          scaVTK->SetTuple1(ii, val);
        }
        else
          scaVTK->SetTuple1(ii, output[ii]);
      }
    //}

    //cout << " extrapolateFlag = " << extrapolateFlag << endl;

    return;
}


void  GrowthFEM::projectInternalVariables(bool extrapolateFlag, int vartype, int varindex, int index)
{
    return;
}




void  GrowthFEM::projectFromElemsToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
{
    int  ee, ii;

    VectorXd  Flocal, solnTemp(nNode);

    if(rhsPostproc.rows() != nNode )
      rhsPostproc.resize(nNode);

    rhsPostproc.setZero();

    cout << " RhsToPostprocess " << endl;

    // compute the stresses at element nodes
    for(ee=0; ee<nElem; ee++)
    {
      elems[ee]->RhsToPostprocess(vartype, varindex, index, rhsPostproc);
    }

    cout << " RhsToPostprocess " << endl;

    SimplicialLDLT<SparseMatrix<double> > solver;

    //SuperLU<SparseMatrixXd> solver;

    //BiCGSTAB<SparseMatrixXd, IncompleteLUT<double> > solver;

    //solver.preconditioner().setDroptol(1.0e-3);
    //solver.preconditioner().setFillfactor(2);
    //solver.setMaxIterations(500);
    //solver.setTolerance(1.0e-8);

    solver.compute(spmtxPostproc);

    double rNorm = rhsPostproc.norm(), val;

    cout << " rNorm = " << rNorm << endl;

    if(rNorm > 1.0e-6)
    {
      solnTemp = solver.solve(rhsPostproc);
      //cout << solver.info() << '\t' << solver.error() << '\t' << solver.iterations() << endl;
    }
    else
      solnTemp.setZero();


    for(ii=0; ii<nNode; ii++)
    {
      if( midnodeData[ii][0] )
      {
        val  = 0.50*solnTemp[ii];
        val += 0.25*solnTemp[midnodeData[ii][1]];
        val += 0.25*solnTemp[midnodeData[ii][2]];

        scaVTK->SetTuple1(ii, val);
      }
      else
        scaVTK->SetTuple1(ii, solnTemp[ii]);
    }

    return;
}









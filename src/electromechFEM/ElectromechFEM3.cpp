
#include "ElectromechFEM.h"
#include "MyTime.h"
#include "Files.h"

//#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>


extern MyTime myTime;
extern Files files;



void ElectromechFEM::plotGeom(int a1, bool b1, int c1, bool d1, int* ffff)
{
    cout <<  " ElectromechFEM::plotGeom ... STARTED" <<  endl;
    //solverEigen->printMatrixPatternToFile();

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
    vtkSmartPointer<vtkTriQuadraticHexahedron>  hexa27VTK    =  vtkSmartPointer<vtkTriQuadraticHexahedron>::New();

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

      if( elems[0]->nodeNums.size() == 3 ) // tria
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

    sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-Geom.vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();

    cout <<  " ElectromechFEM::plotGeom ... STARTED" <<  endl;

    return;
}




void  ElectromechFEM::postProcess(int vartype, int varindex, int index, bool extrapolateFlag, double umin, double umax, int* resln)
{
    //elementDiffStiffTest(0.000001, 9, 6, 6, true);

    cout << " ElectromechFEM::postProcess " << endl;

    int  dd, ii, jj, kk, ll, nlocal, ind1, ind2, e, ee, count, gcount, ind;
    int  n1, n2, n3, n4, n5, n6;
    double vec[3], vec2[3], xx, yy, zz, val;
    //                        {0,1,2,3,4,5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    int  hexa27nodemap[27] =  {0,1,2,3,4,5,6,7,8,11,16, 9,17,10,18,19,12,15,13,14,24,22,20,21,23,25,26};
    //                        {0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15,16,17};
    int  wedge18nodemap[18] = {0,1,2,3,4,5,6,8,12,7,13,14, 9,11,10,15,17,16};

    vtkIdType pt[4];

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
    vtkSmartPointer<vtkTriQuadraticHexahedron>  hexa27VTK    =  vtkSmartPointer<vtkTriQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           dispVTK       = vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           veloVTK       = vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           epotVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           tempVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTK   =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTK2  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           strainVTK     =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           stressVTK     =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    int idd = SolnData.ElemProp[0]->id;

    dispVTK->SetNumberOfComponents(3);
    dispVTK->SetNumberOfTuples(nNode);
    veloVTK->SetNumberOfComponents(3);
    veloVTK->SetNumberOfTuples(nNode);

    epotVTK->SetNumberOfTuples(nNode);
    presVTK->SetNumberOfTuples(nNode);
    tempVTK->SetNumberOfTuples(nNode);

    //cellDataVTK->SetName("stress");
    //cellDataVTK->SetNumberOfTuples(nElem);
    dispVTK->SetName("disp");
    veloVTK->SetName("velocity");
    epotVTK->SetName("potential");
    presVTK->SetName("pressure");
    tempVTK->SetName("temperature");


    // prepare points and cells
    if(ndim == 2)
    {
      vec[2] = 0.0;
      for(ii=0;ii<nNode;ii++)
      {
        if( midnodeData[ii][0] ) // if midnode
        {
          /*
          xx  = 0.5*GeomData.NodePosOrig[ii][0];
          xx += 0.25*GeomData.NodePosOrig[midnodeData[ii][1]][0];
          xx += 0.25*GeomData.NodePosOrig[midnodeData[ii][2]][0];

          yy  = 0.5*GeomData.NodePosOrig[ii][1];
          yy += 0.25*GeomData.NodePosOrig[midnodeData[ii][1]][1];
          yy += 0.25*GeomData.NodePosOrig[midnodeData[ii][2]][1];
          */

          //
          xx  = 0.5*GeomData.NodePosNew[ii][0];
          xx += 0.25*GeomData.NodePosNew[midnodeData[ii][1]][0];
          xx += 0.25*GeomData.NodePosNew[midnodeData[ii][2]][0];

          yy  = 0.5*GeomData.NodePosNew[ii][1];
          yy += 0.25*GeomData.NodePosNew[midnodeData[ii][1]][1];
          yy += 0.25*GeomData.NodePosNew[midnodeData[ii][2]][1];
          //

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);


          n1 = ii*ndof;
          n2 = midnodeData[ii][1]*ndof;
          n3 = midnodeData[ii][2]*ndof;

          // displacement
          vec[0]  = 0.50*SolnData.var1[n1];
          vec[0] += 0.25*SolnData.var1[n2];
          vec[0] += 0.25*SolnData.var1[n3];

          vec[1]  = 0.50*SolnData.var1[n1+1];
          vec[1] += 0.25*SolnData.var1[n2+1];
          vec[1] += 0.25*SolnData.var1[n3+1];

          dispVTK->InsertTuple(ii, vec);

          // velocity
          vec[0]  = 0.50*SolnData.var1Dot[n1];
          vec[0] += 0.25*SolnData.var1Dot[n2];
          vec[0] += 0.25*SolnData.var1Dot[n3];

          vec[1]  = 0.50*SolnData.var1Dot[n1+1];
          vec[1] += 0.25*SolnData.var1Dot[n2+1];
          vec[1] += 0.25*SolnData.var1Dot[n3+1];

          veloVTK->InsertTuple(ii, vec);
        }
        else
        {
          //xx = GeomData.NodePosOrig[ii][0];
          //yy = GeomData.NodePosOrig[ii][1];

          xx = GeomData.NodePosNew[ii][0];
          yy = GeomData.NodePosNew[ii][1];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          kk = ii*ndof;

          // displacement
          vec[0] = SolnData.var1[kk];
          vec[1] = SolnData.var1[kk+1];

          dispVTK->InsertTuple(ii, vec);

          // velocity
          vec[0] = SolnData.var1Dot[kk];
          vec[1] = SolnData.var1Dot[kk+1];

          veloVTK->InsertTuple(ii, vec);
        }
      }

      if(elems[0]->nodeNums.size() == 3)
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
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
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

          // displacement
          vec[0]  = 0.50*SolnData.var1[n1];
          vec[0] += 0.25*SolnData.var1[n2];
          vec[0] += 0.25*SolnData.var1[n3];

          vec[1]  = 0.50*SolnData.var1[n1+1];
          vec[1] += 0.25*SolnData.var1[n2+1];
          vec[1] += 0.25*SolnData.var1[n3+1];

          vec[2]  = 0.50*SolnData.var1[n1+2];
          vec[2] += 0.25*SolnData.var1[n2+2];
          vec[2] += 0.25*SolnData.var1[n3+2];

          dispVTK->InsertTuple(ii, vec);

          // velocity
          vec[0]  = 0.50*SolnData.var1Dot[n1];
          vec[0] += 0.25*SolnData.var1Dot[n2];
          vec[0] += 0.25*SolnData.var1Dot[n3];

          vec[1]  = 0.50*SolnData.var1Dot[n1+1];
          vec[1] += 0.25*SolnData.var1Dot[n2+1];
          vec[1] += 0.25*SolnData.var1Dot[n3+1];

          vec[2]  = 0.50*SolnData.var1Dot[n1+2];
          vec[2] += 0.25*SolnData.var1Dot[n2+2];
          vec[2] += 0.25*SolnData.var1Dot[n3+2];

          veloVTK->InsertTuple(ii, vec);
        }
        else
        {
          //xx = GeomData.NodePosOrig[ii][0];
          //yy = GeomData.NodePosOrig[ii][1];
          //zz = GeomData.NodePosOrig[ii][2];

          xx = GeomData.NodePosNew[ii][0];
          yy = GeomData.NodePosNew[ii][1];
          zz = GeomData.NodePosNew[ii][2];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          kk = ii*ndof;

          // displacement
          vec[0] = SolnData.var1[kk];
          vec[1] = SolnData.var1[kk+1];
          vec[2] = SolnData.var1[kk+2];

          dispVTK->InsertTuple(ii, vec);

          // velocity
          vec[0] = SolnData.var1Dot[kk];
          vec[1] = SolnData.var1Dot[kk+1];
          vec[2] = SolnData.var1Dot[kk+2];

          veloVTK->InsertTuple(ii, vec);
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

    // pressure
    ////////////////////////////

    if(presDegree == 1)
    {
        for(ii=0;ii<nNode;ii++)
        {
          val = SolnData.var2[ii];

          if( midnodeData[ii][0] )
          {
            val  = 0.5*SolnData.var2[midnodeData[ii][1]];
            val += 0.5*SolnData.var2[midnodeData[ii][2]];
          }
          presVTK->SetTuple1(ii, val);
        }

        // face center node for the quadratic quad elements
        if( idd == 1012 )
        {
          for(ee=0;ee<nElem;ee++)
          {
            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];

            presVTK->SetTuple1(elems[ee]->nodeNums[8], val);
          }
        }
        else if( idd == 1052 ) // Hexa element
        {
          for(ee=0;ee<nElem;ee++)
          {
            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];

            presVTK->SetTuple1(elems[ee]->nodeNums[20], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[21], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[7]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];

            presVTK->SetTuple1(elems[ee]->nodeNums[22], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[23], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[7]];

            presVTK->SetTuple1(elems[ee]->nodeNums[24], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[7]];

            presVTK->SetTuple1(elems[ee]->nodeNums[25], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[10]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[12]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[14]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[15]];

            presVTK->SetTuple1(elems[ee]->nodeNums[26], val);
          }
        }
    }
    else
    {
      //cout << vartype << '\t' << varindex << '\t' << index << endl;
      projectStresses(extrapolateFlag, vartype, varindex, index);
      //projectStresses(extrapolateFlag, 4, 10, index);
    }
    //projectStresses(extrapolateFlag, 4, 20, index);


    cellDataVTK->SetName("pressure");
    cellDataVTK->SetNumberOfTuples(nElem);
    for(ee=0;ee<nElem;ee++)
    {
      elems[ee]->elementContourplot(vartype, varindex, index);

      cellDataVTK->InsertTuple1(ee, elems[ee]->vals2project[0]);
    }

    // electric potential
    ////////////////////////////

    if(epotDOF !=  0)
    {
      if(epotDegree == dispDegree)
      {
        for(ii=0; ii<nNode; ii++)
        {
          val = SolnData.var3[ii];

          if( midnodeData[ii][0] )
          {
            val  = 0.50*SolnData.var3[ii];
            val += 0.25*SolnData.var3[midnodeData[ii][1]];
            val += 0.25*SolnData.var3[midnodeData[ii][2]];
          }

          epotVTK->SetTuple1(ii, val);
        }
      }
      else if(epotDegree == (dispDegree-1))
      {
        for(ii=0; ii<nNode; ii++)
        {
          val = SolnData.var3[ii];

          //  for the mid nodes the value is computed as the average of the end nodes
          //  this is only for post-processing purpose using quadratic elements
          if( midnodeData[ii][0] )
          {
            val  = 0.5*SolnData.var3[midnodeData[ii][1]];
            val += 0.5*SolnData.var3[midnodeData[ii][2]];
          }

          epotVTK->SetTuple1(ii, val);
        }
      }
      else 
      {
        cerr <<  " 'dispDegree' and 'epotDegree' not compatible " << endl;
      }
    }


    // temperature
    ////////////////////////////

    if(tempDOF !=  0)
    {
      if(tempDegree == dispDegree)
      {
        for(ii=0; ii<nNode; ii++)
        {
          val = SolnData.var4[ii];

          if( midnodeData[ii][0] )
          {
            val  = 0.50*SolnData.var4[ii];
            val += 0.25*SolnData.var4[midnodeData[ii][1]];
            val += 0.25*SolnData.var4[midnodeData[ii][2]];
          }

          tempVTK->SetTuple1(ii, val);
        }
      }
      else if(tempDegree == (dispDegree-1))
      {
        for(ii=0; ii<nNode; ii++)
        {
          val = SolnData.var4[ii];

          //  for the mid nodes the value is computed as the average of the end nodes
          //  this is only for post-processing purpose using quadratic elements
          if( midnodeData[ii][0] )
          {
            val  = 0.5*SolnData.var4[midnodeData[ii][1]];
            val += 0.5*SolnData.var4[midnodeData[ii][2]];
          }

          tempVTK->SetTuple1(ii, val);
        }
      }
      else 
      {
        cerr <<  " 'dispDegree' and 'tempDegree' not compatible " << endl;
      }
    }


    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetVectors(dispVTK);
    //uGridVTK->GetPointData()->AddArray(veloVTK);
    uGridVTK->GetPointData()->SetScalars(presVTK);

    if(epotDOF !=  0) uGridVTK->GetPointData()->AddArray(epotVTK);
    if(tempDOF !=  0) uGridVTK->GetPointData()->AddArray(tempVTK);

    uGridVTK->GetCellData()->SetScalars(cellDataVTK);


    // create a write object and write uGridVTK to it
    char fname[200];
    sprintf(fname,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);

#if VTK_MAJOR_VERSION == 5
    writerUGridVTK->SetInput(uGridVTK);
#else
    writerUGridVTK->SetInputData(uGridVTK);
#endif

    writerUGridVTK->Write();

    return;
}


void  ElectromechFEM::projectFromElemsToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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

        //scaVTK->SetTuple1(ii, val);
      }
      //else
        //scaVTK->SetTuple1(ii, solnTemp[ii]);
    }

    return;
}



void  ElectromechFEM::projectStrains(bool extrapolateFlag, int vartype, int varindex, int index)
{
    return;
}



void  ElectromechFEM::projectStresses(bool extrapolateFlag, int vartype, int varindex, int index)
{
    int  ee, ii;

    //cout << " extrapolateFlag = " << extrapolateFlag << endl;

    // compute the stresses at element nodes
    for(ee=0; ee<nElem; ee++)
    {
      elems[ee]->projectToNodes(extrapolateFlag, vartype, varindex, index);
    }

    // add element nodal contributions to the global list of nodes
    //if( stressVTK->GetNumberOfTuples() != nNode );
    //{
      //stressVTK->SetNumberOfComponents(9);
      //stressVTK->SetNumberOfTuples(nNode);
    //}

    vector<int>  nodeNums;
    vector<double>  output(nNode, 0.0);
    double  volu, val;

    for(ee=0; ee<nElem; ee++)
    {
      nodeNums = elems[ee]->nodeNums;

      elems[ee]->computeVolume(true);

      volu = elems[ee]->getVolume();

      for(ii=0; ii<nodeNums.size(); ii++)
      {
        output[nodeNums[ii]] += volu*elems[ee]->vals2project[ii];
      }
    }

    for(ee=0; ee<nNode; ee++)
    {
      volu = 0.0;
      for(ii=0; ii<node_elem_conn[ee].size(); ii++)
      {
        volu += elems[node_elem_conn[ee][ii]]->getVolume();
      }

      output[ee] /= volu;
    }

    // 
    // for Bernstein elements, the values for midnodes need to be interpolated
    for(ii=0; ii<nNode; ii++)
    {
        if( midnodeData[ii][0] )
        {
          val  = 0.50*output[ii];
          val += 0.25*output[midnodeData[ii][1]];
          val += 0.25*output[midnodeData[ii][2]];

          //val  = 2.0*output[ii];
          //val -= 0.25*output[midnodeData[ii][1]];
          //val -= 0.25*output[midnodeData[ii][2]];

          presVTK->SetTuple1(ii, val);
        }
        else
          presVTK->SetTuple1(ii, output[ii]);
    }
    //

    return;
}


void  ElectromechFEM::projectInternalVariables(bool extrapolateFlag, int vartype, int varindex, int index)
{
    return;
}


// The arc length function solves the 2nd order equation w.r.t. lambIncrIter
// Returns two values as dlambda1 and dlambda2
inline  void   computeArclengthIncrement(double psi,  double dl, VectorXd& solnIncrStep, VectorXd& solnIncrIter1, VectorXd& solnIncrIter2, double lambIncrStep, VectorXd& forceVec, double* dlambda)

{
    double  psi2 = psi*psi;
    double  forceMagn = forceVec.dot(forceVec);

    //Calculate the coefficients of the polynomial

    double  a = solnIncrIter2.dot(solnIncrIter2) + psi2*forceMagn;

    VectorXd  vecTemp = solnIncrStep+solnIncrIter1;

    double  b = 2.0*vecTemp.dot(solnIncrIter2) + lambIncrStep*psi2*forceMagn;

    double  c = vecTemp.dot(vecTemp) + lambIncrStep*lambIncrStep*psi2*forceMagn - dl*dl;

    //print('%f %f %f' % (a, b, c))

    double  det = b*b-4.0*a*c;

    if(det>0.0)
    {
        dlambda[0] = (-b + sqrt(det))/2.0/a;
        dlambda[1] = (-b - sqrt(det))/2.0/a;
    }
    else
    {
        dlambda[0] = 0.0;
        dlambda[1] = 0.0;

        cout <<  "Possible issue in Arc Length equation " <<  endl;
    }

    return;
}


template <typename T> int signum(T val) {
    return (T(0) < val) - (val < T(0));
}


int  ElectromechFEM::solveWithArclength()
{
    double  tol = 1.0e-9;
    double  psi = 1.0;
    double  dl  = 0.05;

    VectorXd  forceVec(totalDOF), soln(totalDOF), solnIncrStep(totalDOF), solnIncrStep0(totalDOF);
    VectorXd  solnIncrIter(totalDOF), solnIncrIter1(totalDOF), solnIncrIter2(totalDOF);
    VectorXd  solnIncrTemp1(totalDOF), solnIncrTemp2(totalDOF);
    forceVec.setZero();
    soln = forceVec;

    double  dellam1, dellam2, dellam, lambIncrStep, lambIncrStep0, lambIncrIter, lamb=0.0, dlambda[2], det;
    double  aux1, aux2, aux3, aux4, dot1, dot2, fcheck, daomag;

    int  riks = 20000, maxiter = 10, iters=0;
    int  i, iloop;


for(i = 0; i<riks; i++)
    {
        if(lamb >= 1.0)
          break;


        // Increment starts; Set all variations=0


        solnIncrIter.setZero();
        solnIncrIter1.setZero();
        solnIncrIter2.setZero();
        solnIncrStep.setZero();

        dellam1 = 0.0;
        dellam2 = 0.0;
        dellam  = 0.0;
        lambIncrStep = 0.0;


        // calculate stiffness and residual

        // df,dfinv=dfcn((soln+solnIncrStep),th0,(lamb+lambIncrStep))

        calcStiffnessAndResidual();

        // compute solnIncrIter2

        //solnIncrIter2 = solverEigen->solve();              //(forceVec);


        computeArclengthIncrement(psi, dl, solnIncrStep, solnIncrIter1, solnIncrIter2, lambIncrStep, forceVec, dlambda);

        //det = np.linalg.det(df);

        if( signum(det) == signum(dlambda[0]) )
            lambIncrIter = dlambda[0];
        else
            lambIncrIter = dlambda[1];


        solnIncrIter = solnIncrIter1 + lambIncrIter*solnIncrIter2;

        solnIncrStep = solnIncrStep + solnIncrIter;
        lambIncrStep = lambIncrStep + lambIncrIter;

        //f = fcn((soln+solnIncrStep),th0,(lamb+lambIncrStep));

        //fcheck=sqrt(np.dot(f,f));

        if( fcheck<tol )
        {
            soln = soln + solnIncrStep;
            lamb = lamb + lambIncrStep;

            //write(str(soln[0])+' '+str(lamb)+"\n")

            solnIncrStep0 = solnIncrStep;
            lambIncrStep0 = lambIncrStep;
            iloop = 0;
        }
        else
        {
            iters = 0;

            while( fcheck>tol )
            {
                iters += 1;

                solnIncrStep = solnIncrStep;
                lambIncrStep = lambIncrStep;

                //f = fcn((soln+solnIncrStep),th0,(lamb+lambIncrStep));

                //df,dfinv = dfcn((soln+solnIncrStep),th0,(lamb+lambIncrStep));

                //solnIncrIter1 = -np.dot(dfinv,f);
                //solnIncrIter2 = np.dot(dfinv,forceVec);

                computeArclengthIncrement(psi, dl, solnIncrStep, solnIncrIter1, solnIncrIter2, lambIncrStep, forceVec, dlambda);

                solnIncrTemp1 = solnIncrIter1 + dlambda[0]*solnIncrIter2;
                solnIncrTemp2 = solnIncrIter1 + dlambda[1]*solnIncrIter2;

                //det=np.linalg.det(df);

                daomag = solnIncrStep0.dot(solnIncrStep0);
                
                if( daomag==0 )
                {
                    if( signum(lambIncrStep+dlambda[0]) == signum(det) )
                    {
                        solnIncrIter = solnIncrTemp1;
                        lambIncrIter = dlambda[0];
                    }
                    else
                    {
                        solnIncrIter = solnIncrTemp2;
                        lambIncrIter = dlambda[1];
                    }
                }
                else
                {
                    aux1 = solnIncrStep0.dot(solnIncrStep+solnIncrTemp1);
                    aux2 = solnIncrStep0.dot(solnIncrStep+solnIncrTemp2);

                    aux3 = lambIncrStep*(lambIncrStep+dlambda[0])*forceVec.dot(forceVec);
                    aux4 = lambIncrStep*(lambIncrStep+dlambda[1])*forceVec.dot(forceVec);

                    dot1 = aux1 + psi*psi*aux3;
                    dot2 = aux2 + psi*psi*aux4;

                    if( dot1>dot2 )
                    {
                        solnIncrIter = solnIncrTemp1;
                        lambIncrIter = dlambda[0];
                    }
                    else
                    {
                        solnIncrIter = solnIncrTemp2;
                        lambIncrIter = dlambda[1];
                    }

                    if( dlambda[0] == dlambda[1] )
                    {
                        solnIncrIter = solnIncrTemp1;
                        lambIncrIter = dlambda[0];
                    }

                    solnIncrStep = solnIncrStep+solnIncrIter;
                    lambIncrStep = lambIncrStep+lambIncrIter;

                    //f = fcn((soln+solnIncrStep),th0,(lamb+lambIncrStep));

                    //fcheck=np.linalg.norm(f)


                    if( iters>maxiter )
                    {
                        iters=maxiter+1;
                        break;
                    }
                }

                if( iters>maxiter )
                {
                    cout <<  "Convergence cannot achieved within " << maxiter << " iterations" << endl;
                    cout <<  "Program stops" <<  endl;
                    exit(-2);
                }
                else
                {
                    soln += solnIncrStep;
                    solnIncrStep0 = solnIncrStep;

                    lamb += lambIncrStep;
                    lambIncrStep0 = lambIncrStep;

                    //iout.write(str(soln[0])+' '+str(lamb)+"\n")
                }
            }
            cout <<  "Riks increment " <<  i << "completed successfully with " << iters << endl;
        } // else of if( fcheck<tol )

    } // for(i = 0; i<riks; i++)

    cout << "The program completed successfully" << endl;

    return 0;
}





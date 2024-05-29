
#include "MagnetomechFEM.h"
#include "MyTime.h"
#include "Files.h"
#include "TimeFunction.h"


//#include <Eigen/SuperLUSupport>
#include <Eigen/SparseExtra>
#include <Eigen/IterativeSolvers>


extern MyTime myTime;
extern Files files;
extern List<TimeFunction> timeFunction;



void MagnetomechFEM::plotGeom()
{
    cout <<  " MagnetomechFEM::plotGeom ... STARTED" <<  endl;
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

      cout << elems[0]->nodeNums.size() << endl;

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
        cout << "Wrong element type for MagnetomechFEM::plotGeom "  << endl;
        exit(1);
      }
    }

    uGridVTK->SetPoints(pointsVTK);


    cellDataVTKmatltype->SetName("matlType");
    cellDataVTKmatltype->SetNumberOfTuples(nElem);
    cellDataVTKelemtype->SetName("elemType");
    cellDataVTKelemtype->SetNumberOfTuples(nElem);

    for(ee=0;ee<nElem;ee++)
    {
        cellDataVTKmatltype->InsertTuple1(ee, elems[ee]->getMaterialTypeNumber() );
        cellDataVTKelemtype->InsertTuple1(ee, elems[ee]->getElementTypeNumber() );
    }

    uGridVTK->GetCellData()->AddArray(cellDataVTKmatltype);
    uGridVTK->GetCellData()->AddArray(cellDataVTKelemtype);

    if(!COUPLED_PROBLEM)
    {
      cellDataVTKBresi->SetName("Br");
      cellDataVTKBresi->SetNumberOfComponents(3);
      cellDataVTKBresi->SetNumberOfTuples(nElem);

      cellDataVTKBappl->SetName("Bappl");
      cellDataVTKBappl->SetNumberOfComponents(3);
      cellDataVTKBappl->SetNumberOfTuples(nElem);

      VectorXd  vecTemp(3);

      for(ee=0;ee<nElem;ee++)
      {
        vecTemp = elems[ee]->MatlData->getResiMagnfield();
        cellDataVTKBresi->InsertTuple(ee, &vecTemp(0) );

        vecTemp = elems[ee]->MatlData->getApplMagnfield();
        cellDataVTKBappl->InsertTuple(ee, &vecTemp(0) );
      }

      uGridVTK->GetCellData()->AddArray(cellDataVTKBresi);
      uGridVTK->GetCellData()->AddArray(cellDataVTKBappl);
    }

    // create a write object and write uGridVTK to it

    char fname[200];
    //sprintf(fname,"%s%s", files.Ofile.asCharArray(),"-Geom.vtu");
    sprintf(fname,"%s%s%s%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(), "-Geom.vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    cout <<  " MagnetomechFEM::plotGeom ... ENDED" <<  endl;

    return;
}




void  MagnetomechFEM::postProcess()
{
    if( (filecount % outputfrequency) != 0 )
      return;

    double  timerStart = MPI_Wtime();

    //elementDiffStiffTest(0.000001, 9, 6, 6, true);

    cout << " MagnetomechFEM::postProcess " << endl;

    bool extrapolateFlag = false;
    int  vartype = 4;
    int  varindex = 8;
    int  index = 1;


    int  dd, ii, jj, kk, ll, nlocal, ind1, ind2, e, ee, count, gcount, ind;
    int  n1, n2, n3, n4, n5, n6;
    double vec[3], vec2[3], xx, yy, zz, val, geom[3];
    //                        {0,1,2,3,4,5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    int  hexa27nodemap[27] =  {0,1,2,3,4,5,6,7,8,11,16, 9,17,10,18,19,12,15,13,14,24,22,20,21,23,25,26};
    //                        {0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15,16,17};
    int  wedge18nodemap[18] = {0,1,2,3,4,5,6,8,12,7,13,14, 9,11,10,15,17,16};

    vtkIdType pt[4];
    vector<int>  nodeNumsElem;

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

    int idd = SolnData.ElemProp[0]->id;
    cout << "idd = " << idd << endl;

    dispVTK->SetName("disp");
    dispVTK->SetNumberOfComponents(3);
    dispVTK->SetNumberOfTuples(nNode);

    veloVTK->SetName("velocity");
    veloVTK->SetNumberOfComponents(3);
    veloVTK->SetNumberOfTuples(nNode);

    mpotVTK->SetName("potential");
    mpotVTK->SetNumberOfTuples(nNode);

    //cellDataVTK->SetName("stress");
    //cellDataVTK->SetNumberOfTuples(nElem);

    presVTK->SetName("pressure");
    presVTK->SetNumberOfTuples(nNode);

    SxxNodesVTK->SetName("Sxx");
    SxxNodesVTK->SetNumberOfTuples(nNode);
    SyyNodesVTK->SetName("Syy");
    SyyNodesVTK->SetNumberOfTuples(nNode);
    SzzNodesVTK->SetName("Szz");
    SzzNodesVTK->SetNumberOfTuples(nNode);



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

        if( (idd == 2052) || (idd == 2055) )           // Wedge element
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
        if( (idd == 2053) || (idd == 2056) )           // Hexa element
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

    //cout << "pressure "  << endl;

    // pressure
    ////////////////////////////

    if( (presDegree == 1) && (varindex == 10) )
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
        if( (idd == 2002) || (idd == 2004) )
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
        else if( (idd == 2052) || (idd == 2055) )           // Wedge element
        {
          for(ee=0;ee<nElem;ee++)
          {
            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];

            presVTK->SetTuple1(elems[ee]->nodeNums[15], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[16], val);

            val  = 0.25*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.var2[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[17], val);
          }
        }
        else if( (idd == 2053) || (idd == 2056) )           // Hexa element
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

            val  = 0.125*SolnData.var2[elems[ee]->nodeNums[0]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[1]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[2]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[3]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[4]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[5]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[6]];
            val += 0.125*SolnData.var2[elems[ee]->nodeNums[7]];

            presVTK->SetTuple1(elems[ee]->nodeNums[26], val);
          }
        }
    }
    else
    {
      cout << vartype << '\t' << varindex << '\t' << index << endl;
      //projectStresses(extrapolateFlag, vartype, varindex, index);
      //projectStresses(extrapolateFlag, 4, 10, index);
    }
    //projectStresses(extrapolateFlag, 4, 20, index);


    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetVectors(dispVTK);
    uGridVTK->GetPointData()->AddArray(veloVTK);

    VectorXd  output1(nNode), output2(nNode), output3(nNode), output4(nNode);

    //projectStresses(extrapolateFlag, 4,  0, 1, output1);
    //projectStresses(extrapolateFlag, 4,  4, 1, output2);
    //projectStresses(extrapolateFlag, 4,  8, 1, output3);
    projectStresses(extrapolateFlag, 4, 10, 1, output4);

    for(ii=0; ii<nNode; ii++)
    {
        //SxxNodesVTK->SetTuple1(ii, output1[ii]);
        //SyyNodesVTK->SetTuple1(ii, output2[ii]);
        //SzzNodesVTK->SetTuple1(ii, output3[ii]);

        presVTK->SetTuple1(ii, output4[ii]);
    }

    //uGridVTK->GetPointData()->SetScalars(presVTK);
    //uGridVTK->GetPointData()->AddArray(SxxNodesVTK);
    //uGridVTK->GetPointData()->AddArray(SyyNodesVTK);
    uGridVTK->GetPointData()->AddArray(SzzNodesVTK);

    // magnetic potential
    ////////////////////////////

    if(mpotDOF !=  0)
    {
      if(mpotDegree == dispDegree)
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

          mpotVTK->SetTuple1(ii, val);
        }
      }
      else if(mpotDegree == -1)
      {
      }
      else 
      {
        cerr <<  " 'dispDegree' and 'mpotDegree' not compatible " << endl;
      }

      uGridVTK->GetPointData()->AddArray(mpotVTK);
    }



    /////////////////////////////
    // Cell data fields
    /////////////////////////////

    cellDataVTK->SetName("pressure");
    cellDataVTK->SetNumberOfTuples(nElem);

    cellDataVTKmatltype->SetName("matlType");
    cellDataVTKmatltype->SetNumberOfTuples(nElem);
    cellDataVTKelemtype->SetName("elemType");
    cellDataVTKelemtype->SetNumberOfTuples(nElem);

    for(ee=0;ee<nElem;ee++)
    {
        elems[ee]->elementContourplot(4, 10, 1);

        cellDataVTK->InsertTuple1(ee, elems[ee]->vals2project[0]);

        cellDataVTKmatltype->InsertTuple1(ee, elems[ee]->getMaterialTypeNumber() );
        cellDataVTKelemtype->InsertTuple1(ee, elems[ee]->getElementTypeNumber() );
    }

    uGridVTK->GetCellData()->SetScalars(cellDataVTK);
    uGridVTK->GetCellData()->AddArray(cellDataVTKmatltype);
    //uGridVTK->GetCellData()->AddArray(cellDataVTKelemtype);


    if(!COUPLED_PROBLEM)
    {
      cellDataVTKBresi->SetName("Br");
      cellDataVTKBresi->SetNumberOfComponents(3);
      cellDataVTKBresi->SetNumberOfTuples(nElem);

      cellDataVTKBappl->SetName("Bappl");
      cellDataVTKBappl->SetNumberOfComponents(3);
      cellDataVTKBappl->SetNumberOfTuples(nElem);

      VectorXd  vecTemp;

      double  fact = timeFunction[0].prop;

      for(ee=0;ee<nElem;ee++)
      {
        vecTemp = elems[ee]->MatlData->getResiMagnfield();
        cellDataVTKBresi->InsertTuple(ee, &vecTemp(0) );

        vecTemp = fact*elems[ee]->MatlData->getApplMagnfield();
        cellDataVTKBappl->InsertTuple(ee, &vecTemp(0) );
      }

      uGridVTK->GetCellData()->AddArray(cellDataVTKBresi);
      uGridVTK->GetCellData()->AddArray(cellDataVTKBappl);
    }



    // create a write object and write uGridVTK to it
    char fname[200];
    //sprintf(fname,"%s%s%06d%s", files.Ofile.asCharArray(),"-",filecount, ".vtu");
    sprintf(fname,"%s%s%s%s%06d%s", files.projDir.asCharArray(), "/", files.Ofile.asCharArray(),"-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    double  timerEnd = MPI_Wtime();
    printf("\n\n Elapsed time = %f seconds \n\n", timerEnd - timerStart );


    return;
}


void  MagnetomechFEM::projectStrains(bool extrapolateFlag, int vartype, int varindex, int index, VectorXd&  output)
{
    return;
}



void  MagnetomechFEM::projectStresses(bool extrapolateFlag, int vartype, int varindex, int index, VectorXd&  output)
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
    double  volu, val;

    if(output.rows() != nNode)  output.resize(nNode);
    output.setZero();

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

          output[ii] = val;
          //presVTK->SetTuple1(ii, val);
        }
        //else
          //presVTK->SetTuple1(ii, output[ii]);
    }
    //

    return;
}


void  MagnetomechFEM::projectInternalVariables(bool extrapolateFlag, int vartype, int varindex, int index, VectorXd&  output)
{
    return;
}





void  MagnetomechFEM::projectFromElemsToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
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





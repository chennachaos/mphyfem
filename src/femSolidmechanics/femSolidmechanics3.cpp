
#include "femSolidmechanics.h"
#include "MyTime.h"
#include "ElementBase.h"
#include "TimeFunction.h"


extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime                myTime;
extern bool debug;



void femSolidmechanics::plotGeom()
{
    if(debug) PetscPrintf(PETSC_COMM_WORLD, "\n     femSolidmechanics::plotGeom()  .... STARTED ...\n");
    PetscPrintf(PETSC_COMM_WORLD, "\n     femSolidmechanics::plotGeom()  .... STARTED ...\n");

    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();
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

    vtkSmartPointer<vtkLagrangeTriangle>    triaHigherVTK     =  vtkSmartPointer<vtkLagrangeTriangle>::New();
    vtkSmartPointer<vtkLagrangeTetra>    tetraHigherVTK     =  vtkSmartPointer<vtkLagrangeTetra>::New();

    vtkSmartPointer<vtkFloatArray>           dispVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           veloVTK      =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK2      =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           cellDataVTK  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTK2 =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           strainVTK    =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           stressVTK    =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           procIdVTKcell =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           timeloadstamp  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    int  dd, ii, jj, kk, ll, nlocal, index, ind1, ind2, e, ee, count, gcount, ind, npElem;
    double xx, yy, zz;

    //                       {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int  tria10nodemap[10] = {0, 2, 1, 7, 8, 6, 5, 4, 3, 9};

    //                        {0,1,2,3,4,5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    int  hexa27nodemap[27] =  {0,1,2,3,4,5,6,7,8,11,16, 9,17,10,18,19,12,15,13,14,24,22,20,21,23,25,26};

    //                        {0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15,16,17};
    int  wedge18nodemap[18] = {0,1,2,3,4,5,6,8,12,7,13,14, 9,11,10,15,17,16};

    //                        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int  tetra10nodemap[10] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 8};

    //                        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19};
    int  tetra20nodemap[20] = {0, 2, 3, 1, 9, 8, 14, 15,11,10, 5, 4,12,13,6,7,18,19,16,17};

    vtkIdType pt[50];

    vector<int> nodeNums;


    if(ndim == 2)
    {
      for(ii=0;ii<nNode_global;ii++)
      {
        /*
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
        */
          xx = GeomData.NodePosOrig[ii][0];
          yy = GeomData.NodePosOrig[ii][1];
        //}

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);
      }

      for(ee=0;ee<nElem_global;ee++)
      {
        nodeNums = elems[ee]->nodeNums;
        //printVector(nodeNums);
        npElem = nodeNums.size();

        if( npElem == 2 ) // 2-noded line
        {
          lineVTK->GetPointIds()->SetId(0, nodeNums[0]);
          lineVTK->GetPointIds()->SetId(1, nodeNums[1]);

          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
        }
        else if( npElem == 3 )  // 3-noded tria
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
        else if( npElem == 6 ) // 6-noded tria
        {
          for(ll=0;ll<6;ll++)
            tria6VTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
        else if( npElem == 4 ) // quad, 4 node
        {
          for(ll=0;ll<4;ll++)
            quadVTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
        else if( npElem == 9 ) // quad, 9 node
        {
          for(ll=0;ll<9;ll++)
            quad9VTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
        else if( npElem == 10 ) // tria, 10 node
        {
          triaHigherVTK->GetPointIds()->SetNumberOfIds(10);
          triaHigherVTK->GetPoints()->SetNumberOfPoints(10);
          triaHigherVTK->Initialize();

          for(ll=0;ll<10;ll++)
            triaHigherVTK->GetPointIds()->SetId(tria10nodemap[ll], nodeNums[ll]);

          uGridVTK->InsertNextCell(triaHigherVTK->GetCellType(), triaHigherVTK->GetPointIds());
        }
      }
    }
    else
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
        /*
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
        */
          xx = GeomData.NodePosOrig[ii][0];
          yy = GeomData.NodePosOrig[ii][1];
          zz = GeomData.NodePosOrig[ii][2];
        //}

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);
      }

      for(ee=0; ee<nElem_global; ee++)
      {
        nodeNums = elems[ee]->nodeNums;
        printVector(nodeNums);
        npElem = nodeNums.size();

        if(npElem == 4) // tet, 4-noded
        {
          for(ll=0;ll<4;ll++)
            tetraVTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
        else if(npElem == 8) // hexa, 8-noded
        {
          for(ll=0;ll<8;ll++)
            hexaVTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(hexaVTK->GetCellType(), hexaVTK->GetPointIds());
        }
        else if(npElem == 10) // 10-node tetrahedron
        {
          for(ii=0; ii<10; ii++)
            tetra10VTK->GetPointIds()->SetId(ii, nodeNums[ii]);

          uGridVTK->InsertNextCell(tetra10VTK->GetCellType(), tetra10VTK->GetPointIds());
        }
        else if(npElem == 18) // 18-node wedge
        {
          for(ii=0; ii<18; ii++)
            wedge18VTK->GetPointIds()->SetId(wedge18nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(wedge18VTK->GetCellType(), wedge18VTK->GetPointIds());
        }
        else if(npElem == 20) // 20-node tetrahedron
        {
          tetraHigherVTK->GetPointIds()->SetNumberOfIds(20);
          tetraHigherVTK->GetPoints()->SetNumberOfPoints(20);
          tetraHigherVTK->Initialize();

          for(ii=0; ii<20; ii++)
            tetraHigherVTK->GetPointIds()->SetId(tetra20nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(tetraHigherVTK->GetCellType(), tetraHigherVTK->GetPointIds());
        }
        else if(npElem == 27) // 27-node hexahedron
        {
          for(ii=0; ii<27; ii++)
            hexa27VTK->GetPointIds()->SetId(hexa27nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(hexa27VTK->GetCellType(), hexa27VTK->GetPointIds());
        }
        else
        {
          cout << "Wrong 3D element type for plotGeom "  << endl;
          exit(1);
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);

    // create a write object and write uGridVTK to it

    char fname[200];
    sprintf(fname,"%s%s", inputfilename.c_str(),"-Geom.vtu");

    writerUGridVTK->SetFileName(fname);

    writerUGridVTK->SetInputData(uGridVTK);

    writerUGridVTK->Write();

    if(debug) PetscPrintf(PETSC_COMM_WORLD, "\n     femSolidmechanics::plotGeom()  .... ENDED ...\n");

    return;
}




void  femSolidmechanics::postProcess()
{
    if(debug) cout << " femSolidmechanics::postProcess ... STARTED " << endl;

    int vartype=0, varindex=0, index=0;
    bool extrapolateFlag=false;

    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();
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

    vtkSmartPointer<vtkFloatArray>           dispVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           veloVTK      =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           presVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK2      =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           matlIdVTKcell  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           elemIdVTKcell =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           strainVTK    =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           stressVTK    =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           procIdVTKcell =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkLagrangeTriangle>    triaHigherVTK     =  vtkSmartPointer<vtkLagrangeTriangle>::New();
    vtkSmartPointer<vtkLagrangeTetra>    tetraHigherVTK     =  vtkSmartPointer<vtkLagrangeTetra>::New();

    vtkSmartPointer<vtkFloatArray>           timeloadstamp  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    int  dd, ii, jj, kk, ll, nlocal, ind1, ind2, e, ee, count, gcount, ind, npElem;
    int  n1, n2, n3, n4, n5, n6;
    double vec[3], vec2[3], xx, yy, zz, val, maxdisp = 0.0, maxstress = 0.0;

    //                       {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int  tria10nodemap[10] = {0, 2, 1, 7, 8, 6, 5, 4, 3, 9};

    //                        {0,1,2,3,4,5,6,7,8,9, 10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    int  hexa27nodemap[27] =  {0,1,2,3,4,5,6,7,8,11,16, 9,17,10,18,19,12,15,13,14,24,22,20,21,23,25,26};

    //                        {0,1,2,3,4,5,6,7, 8,9,10,11,12,13,14,15,16,17};
    int  wedge18nodemap[18] = {0,1,2,3,4,5,6,8,12,7,13,14, 9,11,10,15,17,16};

    //                        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
    int  tetra10nodemap[10] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 8};

    //                        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19};
    int  tetra20nodemap[20] = {0, 2, 3, 1, 9, 8, 14, 15,11,10, 5, 4,12,13,6,7,18,19,16,17};

    vtkIdType pt[50];
    vector<int> nodeNums;


    dispVTK->SetName("disp");
    veloVTK->SetName("velocity");
    dispVTK->SetNumberOfComponents(3);
    dispVTK->SetNumberOfTuples(nNode_global);
    veloVTK->SetNumberOfComponents(3);
    veloVTK->SetNumberOfTuples(nNode_global);

    matlIdVTKcell->SetName("MatlID");
    matlIdVTKcell->SetNumberOfTuples(nElem_global);
    elemIdVTKcell->SetName("ElemID");
    elemIdVTKcell->SetNumberOfTuples(nElem_global);

    //printVector(SolnData.disp);
    //printVector(SolnData.pres);

    timeloadstamp->SetName("TIME");
    timeloadstamp->SetNumberOfTuples(1);
    if(ARC_LENGTH)
      timeloadstamp->SetTuple1(0, loadFactor);
    else
      timeloadstamp->SetTuple1(0, myTime.cur);

    int  idd = ElementTypeDataList[0]->getElemTypeNameNum();

    if(ndim == 2)
    {
      vec[2] = 0.0;
      for(ii=0;ii<nNode_global;ii++)
      {
          xx = GeomData.NodePosNew[ii][0];
          yy = GeomData.NodePosNew[ii][1];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          kk = ii*ndof;

          vec[0] = SolnData.disp[kk];
          vec[1] = SolnData.disp[kk+1];

          dispVTK->InsertTuple(ii, vec);

          vec[0] = SolnData.velo[kk];
          vec[1] = SolnData.velo[kk+1];

          //vec[0] = SolnData.veloDotCur[kk];
          //vec[1] = SolnData.veloDotCur[kk+1];

          veloVTK->InsertTuple(ii, vec);
      }

      // center nodes for 9-noded quadrilateral elements
      if( (idd == 221) || (idd == 222) || (idd == 223) )
      {
        for(ee=0;ee<nElem_global;ee++)
        {
          //xx  = 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[0]][0];
          //xx += 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[1]][0];
          //xx += 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[2]][0];
          //xx += 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[3]][0];

          //yy  = 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[0]][1];
          //yy += 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[1]][1];
          //yy += 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[2]][1];
          //yy += 0.2500*GeomData.NodePosNew[elems[ee]->nodeNums[3]][1];

          xx  = GeomData.NodePosNew[elems[ee]->nodeNums[8]][0];
          yy  = GeomData.NodePosNew[elems[ee]->nodeNums[8]][1];

          pointsVTK->SetPoint(elems[ee]->nodeNums[8], xx, yy, 0.0);

          //vec[0]  = 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[0]];
          //vec[0] += 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[1]];
          //vec[0] += 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[2]];
          //vec[0] += 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[3]];

          //vec[1]  = 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[0]+1];
          //vec[1] += 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[1]+1];
          //vec[1] += 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[2]+1];
          //vec[1] += 0.2500*SolnData.disp[ndof*elems[ee]->nodeNums[3]+1];

          vec[0]  = SolnData.disp[ndof*elems[ee]->nodeNums[8]];
          vec[1]  = SolnData.disp[ndof*elems[ee]->nodeNums[8]+1];

          dispVTK->SetTuple(elems[ee]->nodeNums[8], vec);

          //vec[0]  = 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[0]];
          //vec[0] += 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[1]];
          //vec[0] += 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[2]];
          //vec[0] += 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[3]];

          //vec[1]  = 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[0]+1];
          //vec[1] += 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[1]+1];
          //vec[1] += 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[2]+1];
          //vec[1] += 0.2500*SolnData.velo[ndof*elems[ee]->nodeNums[3]+1];

          vec[0]  = SolnData.velo[ndof*elems[ee]->nodeNums[8]];
          vec[1]  = SolnData.velo[ndof*elems[ee]->nodeNums[8]+1];

          veloVTK->SetTuple(elems[ee]->nodeNums[8], vec);
        }
      }

      for(ee=0;ee<nElem_global;ee++)
      {
        matlIdVTKcell->SetTuple1(ee, elems[ee]->matType);
        elemIdVTKcell->SetTuple1(ee, elems[ee]->elmType);

        nodeNums = elems[ee]->nodeNums;
        //printVector(nodeNums);
        npElem = nodeNums.size();

        if( npElem == 2 ) // 2-noded line
        {
          lineVTK->GetPointIds()->SetId(0, nodeNums[0]);
          lineVTK->GetPointIds()->SetId(1, nodeNums[1]);

          uGridVTK->InsertNextCell(lineVTK->GetCellType(), lineVTK->GetPointIds());
        }
        else if( npElem == 3 )  // 3-noded tria
        {
          for(ll=0;ll<3;ll++)
            triaVTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
        else if( npElem == 6 ) // 6-noded tria
        {
          for(ll=0;ll<6;ll++)
            tria6VTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
        else if( npElem == 4 ) // quad, 4 node
        {
          for(ll=0;ll<4;ll++)
            quadVTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
        else if( npElem == 9 ) // quad, 9 node
        {
          for(ll=0;ll<9;ll++)
            quad9VTK->GetPointIds()->SetId(ll, nodeNums[ll]);

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
        else if( npElem == 10 ) // tria, 10 node
        {
          triaHigherVTK->GetPointIds()->SetNumberOfIds(10);
          triaHigherVTK->GetPoints()->SetNumberOfPoints(10);
          triaHigherVTK->Initialize();

          for(ll=0;ll<10;ll++)
            triaHigherVTK->GetPointIds()->SetId(tria10nodemap[ll], nodeNums[ll]);

          uGridVTK->InsertNextCell(triaHigherVTK->GetCellType(), triaHigherVTK->GetPointIds());
        }
      }

    }
    else //if(ndim == 3)
    {
      for(ii=0;ii<nNode_global;ii++)
      {
          xx = GeomData.NodePosNew[ii][0];
          yy = GeomData.NodePosNew[ii][1];
          zz = GeomData.NodePosNew[ii][2];

          pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

          vertexVTK->GetPointIds()->SetId(0, pt[0]);

          kk = ii*ndof;

          vec[0] = SolnData.disp[kk];
          vec[1] = SolnData.disp[kk+1];
          vec[2] = SolnData.disp[kk+2];

          dispVTK->InsertTuple(ii, vec);

          vec[0] = SolnData.velo[kk];
          vec[1] = SolnData.velo[kk+1];
          vec[2] = SolnData.velo[kk+2];

          veloVTK->InsertTuple(ii, vec);
      }


      for(ee=0; ee<nElem_global; ee++)
      {
        matlIdVTKcell->SetTuple1(ee, elems[ee]->matType);
        elemIdVTKcell->SetTuple1(ee, elems[ee]->elmType);

        nodeNums = elems[ee]->nodeNums;
        //printVector(nodeNums);
        npElem = nodeNums.size();

        if(npElem == 4) // tet, 4-noded
        {
          for(ii=0; ii<4; ii++)
            tetraVTK->GetPointIds()->SetId(ii, nodeNums[ii]);

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
        else if(npElem == 8) // hexa, 8-noded
        {
          for(ii=0; ii<8; ii++)
            hexaVTK->GetPointIds()->SetId(ii, nodeNums[ii]);

          uGridVTK->InsertNextCell(hexaVTK->GetCellType(), hexaVTK->GetPointIds());
        }
        else if(npElem == 10) // 10-node tetrahedron
        {
          for(ii=0; ii<10; ii++)
            tetra10VTK->GetPointIds()->SetId(tetra10nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(tetra10VTK->GetCellType(), tetra10VTK->GetPointIds());
        }
        else if(npElem == 11) // 11-node tetrahedron
        {
          for(ii=0; ii<10; ii++)
            tetra10VTK->GetPointIds()->SetId(tetra10nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(tetra10VTK->GetCellType(), tetra10VTK->GetPointIds());
        }
        else if(npElem == 18) // 18-node wedge
        {
          for(ii=0; ii<18; ii++)
            wedge18VTK->GetPointIds()->SetId(wedge18nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(wedge18VTK->GetCellType(), wedge18VTK->GetPointIds());
        }
        else if(npElem == 20) // 20-node tetrahedron
        {
          tetraHigherVTK->GetPointIds()->SetNumberOfIds(20);
          tetraHigherVTK->GetPoints()->SetNumberOfPoints(20);
          tetraHigherVTK->Initialize();

          for(ii=0; ii<20; ii++)
            tetraHigherVTK->GetPointIds()->SetId(tetra20nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(tetraHigherVTK->GetCellType(), tetraHigherVTK->GetPointIds());
        }
        else if(npElem == 27) // 27-node hexahedron
        {
          for(ii=0; ii<27; ii++)
            hexa27VTK->GetPointIds()->SetId(hexa27nodemap[ii], nodeNums[ii]);

          uGridVTK->InsertNextCell(hexa27VTK->GetCellType(), hexa27VTK->GetPointIds());
        }
        else
        {
          cout << "Wrong 3D element type for plotGeom "  << endl;
          exit(1);
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetFieldData()->AddArray(timeloadstamp);
    uGridVTK->GetPointData()->SetVectors(dispVTK);
    uGridVTK->GetPointData()->AddArray(veloVTK);
    uGridVTK->GetCellData()->AddArray(elemIdVTKcell);
    uGridVTK->GetCellData()->AddArray(matlIdVTKcell);

    //cellDataVTK->SetName("pres");
    //cellDataVTK->SetNumberOfTuples(nElem_global);
    /*
    for(ee=0;ee<nElem_global;ee++)
    {
      elems[ee]->elementContourplot(vartype, varindex, index);

      cellDataVTK->InsertTuple1(ee, elems[ee]->vals2project[0]);

      //maxstress = max(maxstress, elems[ee]->vals2project[0]);
    }
    */

    if(intVarFlag)
    {
      //scaVTK2->SetName("epstrn");
      //scaVTK2->SetNumberOfTuples(nNode_global);
      //projectStresses(1, 5, 6, index);

      //uGridVTK->GetPointData()->AddArray(scaVTK2);
    }

    presVTK->SetName("pressure");
    presVTK->SetNumberOfTuples(nNode_global);

    cout << vartype << '\t' << varindex << '\t' << index << endl;

    if( (idd == ELEM_SOLID_2D_TRIA7_MIXED1) ||
        (idd == ELEM_SOLID_2D_TRIA10_MIXED1) ||
        (idd == ELEM_SOLID_3D_TETRA11_MIXED1) ||
        (idd == ELEM_SOLID_3D_TETRA20_MIXED1) )
     {
      postProcessPressureForDiscontinuousElements();
    }
    else if(MIXED_ELEMENT)
    {
        if( (idd == 212) || (idd == 214) || (idd == 223) || (idd == 225) || (idd == 227) )
        {
          for(ii=0;ii<nNode_global;ii++)
          {
            val = SolnData.pres[ii];

            if( midnodeData[ii][0] )
            {
              val  = 0.5*SolnData.pres[midnodeData[ii][1]];
              val += 0.5*SolnData.pres[midnodeData[ii][2]];
            }
            presVTK->SetTuple1(ii, val);
          }
        }

        // face center node for the quadratic quad elements
        if( (idd == 223) )
        {
          for(ee=0;ee<nElem_global;ee++)
          {
            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[3]];

            presVTK->SetTuple1(elems[ee]->nodeNums[8], val);
          }
        }
        else if( (idd == 225) ) // Wedge element
        {
          for(ee=0;ee<nElem_global;ee++)
          {
            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[4]];

            presVTK->SetTuple1(elems[ee]->nodeNums[15], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[16], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[17], val);
          }
        }
        else if( (idd == 227) ) // Hexa element
        {
          for(ee=0;ee<nElem_global;ee++)
          {
            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[3]];

            presVTK->SetTuple1(elems[ee]->nodeNums[20], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[21], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[0]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[7]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[4]];

            presVTK->SetTuple1(elems[ee]->nodeNums[22], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[1]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[5]];

            presVTK->SetTuple1(elems[ee]->nodeNums[23], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[2]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[3]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[7]];

            presVTK->SetTuple1(elems[ee]->nodeNums[24], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[4]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[5]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[6]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[7]];

            presVTK->SetTuple1(elems[ee]->nodeNums[25], val);

            val  = 0.25*SolnData.pres[elems[ee]->nodeNums[10]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[12]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[14]];
            val += 0.25*SolnData.pres[elems[ee]->nodeNums[15]];

            presVTK->SetTuple1(elems[ee]->nodeNums[26], val);
          }
        }

        uGridVTK->GetPointData()->SetScalars(presVTK);
    }
    else
    {
        //projectStresses(extrapolateFlag, vartype, varindex, index);
        cout << vartype << '\t' << varindex << '\t' << index << endl;
        //projectFromElemsToNodes(extrapolateFlag, vartype, varindex, index);
    }




    char fname[200];
    sprintf(fname,"%s%s%06d%s", inputfilename.c_str(),"-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    if(debug) cout << " femSolidmechanics::postProcess ... ENDED " << endl;

    return;
}




void  femSolidmechanics::postProcessPressureForDiscontinuousElements()
{
    if(debug) cout << " femSolidmechanics::postProcessPressureForDiscontinuousElements ... STARTED " << endl;


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkLine>                 lineVTK      =  vtkSmartPointer<vtkLine>::New();
    vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkWedge>                wedgeVTK     =  vtkSmartPointer<vtkWedge>::New();
    vtkSmartPointer<vtkHexahedron>           hexaVTK      =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           presVTK      =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           matlIdVTKcell  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           elemIdVTKcell =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           timeloadstamp  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    int  dd, ii, jj, ll, nn, ind;
    double xx, yy, zz;

    vtkIdType pt[50];
    vector<int> nodeNums, nodeNumsPres;

    matlIdVTKcell->SetName("MatlID");
    matlIdVTKcell->SetNumberOfTuples(nElem_global);
    elemIdVTKcell->SetName("ElemID");
    elemIdVTKcell->SetNumberOfTuples(nElem_global);


    presVTK->SetName("pressure");
    presVTK->SetNumberOfTuples(nElem_global*elems[0]->nlbfP);

    timeloadstamp->SetName("TIME");
    timeloadstamp->SetNumberOfTuples(1);
    if(ARC_LENGTH)
      timeloadstamp->SetTuple1(0, loadFactor);
    else
      timeloadstamp->SetTuple1(0, myTime.cur);

  if(ndim == 2)
  {
    zz = 0.0;
    for(int ee=0; ee<nElem_global; ee++)
    {
      matlIdVTKcell->SetTuple1(ee, elems[ee]->matType);
      elemIdVTKcell->SetTuple1(ee, elems[ee]->elmType);

      nodeNums = elems[ee]->nodeNums;

      ind = ee*3;
      for(ll=0; ll<3; ll++)
      {
        nn = nodeNums[ll];

        xx = GeomData.NodePosNew[nn][0];
        yy = GeomData.NodePosNew[nn][1];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

        vertexVTK->GetPointIds()->SetId(0, pt[0]);

        triaVTK->GetPointIds()->SetId(ll, pt[0]);

        presVTK->SetTuple1(ind+ll, SolnData.pres[ind+ll]);
      }

      uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
    }
  }//if(ndim
  else
  {
    for(int ee=0; ee<nElem_global; ee++)
    {
      matlIdVTKcell->SetTuple1(ee, elems[ee]->matType);
      elemIdVTKcell->SetTuple1(ee, elems[ee]->elmType);

      nodeNums = elems[ee]->nodeNums;

      ind = ee*4;
      for(ll=0; ll<4; ll++)
      {
        nn = nodeNums[ll];

        xx = GeomData.NodePosNew[nn][0];
        yy = GeomData.NodePosNew[nn][1];
        zz = GeomData.NodePosNew[nn][2];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);

        vertexVTK->GetPointIds()->SetId(0, pt[0]);

        tetraVTK->GetPointIds()->SetId(ll, pt[0]);

        presVTK->SetTuple1(ind+ll, SolnData.pres[ind+ll]);
      }

      uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
    }

  }


    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetFieldData()->AddArray(timeloadstamp);
    uGridVTK->GetPointData()->SetScalars(presVTK);
    uGridVTK->GetCellData()->AddArray(elemIdVTKcell);
    uGridVTK->GetCellData()->AddArray(matlIdVTKcell);

    char fname[200];
    sprintf(fname,"%s%s%06d%s", inputfilename.c_str(),"-pressure-",filecount, ".vtu");

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    if(debug) cout << " femSolidmechanics::postProcessPressureForDiscontinuousElements ... ENDED " << endl;

  return;
}







void  femSolidmechanics::projectStrains(bool extrapolateFlag, int vartype, int varindex, int index)
{
    return;
}



void  femSolidmechanics::projectStresses(bool extrapolateFlag, int vartype, int varindex, int index)
{
    int  ee, ii;

    //cout << " extrapolateFlag = " << extrapolateFlag << endl;

    // compute the stresses at element nodes
    for(ee=0; ee<nElem_global; ee++)
    {
      elems[ee]->projectToNodes(extrapolateFlag, vartype, varindex, index);
    }

    vector<int>  nodeNums;
    vector<double>  output(nNode_global, 0.0);
    double  volu, val;

    for(ee=0; ee<nElem_global; ee++)
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

    for(ee=0; ee<nNode_global; ee++)
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
      for(ii=0; ii<nNode_global; ii++)
      {
        if( midnodeData[ii][0] )
        {
          val  = 0.50*output[ii];
          val += 0.25*output[midnodeData[ii][1]];
          val += 0.25*output[midnodeData[ii][2]];

          //scaVTK->SetTuple1(ii, val);
        }
        //else
          //scaVTK->SetTuple1(ii, output[ii]);
      }
    //}

    //cout << " extrapolateFlag = " << extrapolateFlag << endl;

    return;
}


void  femSolidmechanics::projectInternalVariables(bool extrapolateFlag, int vartype, int varindex, int index)
{
    return;
}




void  femSolidmechanics::projectFromElemsToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
{
    int  ee, ii;

    VectorXd  Flocal, solnTemp(nNode_global);

    cout << " RhsToPostprocess " << endl;
/*
    // compute the stresses at element nodes
    //for(ee=0; ee<nElem_global; ee++)
      //elems[ee]->RhsToPostprocess(vartype, varindex, index, rhsPostproc);

    cout << " RhsToPostprocess " << endl;

    for(ii=0; ii<nNode_global; ii++)
    {
      if( midnodeData[ii][0] )
      {
        val  = 0.50*solnTemp[ii];
        val += 0.25*solnTemp[midnodeData[ii][1]];
        val += 0.25*solnTemp[midnodeData[ii][2]];

        //scaVTK->SetTuple1(ii, val);
      }
      else
        //scaVTK->SetTuple1(ii, solnTemp[ii]);
    }
*/
    return;
}









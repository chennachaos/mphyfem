
#include "femMagnetomech.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include <boost/algorithm/string.hpp>

extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime myTime;
extern bool debug;



void femMagnetomech::plotGeom()
{
    cout <<  " femMagnetomech::plotGeom ... STARTED" <<  endl;
    //solverEigen->printMatrixPatternToFile();
/*
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
          xx = GeomData.NodePosOrig[ii][0];
          yy = GeomData.NodePosOrig[ii][1];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, 0.0);
      }

    }
    else
    {
      for(ii=0;ii<GeomData.NodePosOrig.size();ii++)
      {
          xx = GeomData.NodePosOrig[ii][0];
          yy = GeomData.NodePosOrig[ii][1];
          zz = GeomData.NodePosOrig[ii][2];

        pt[0] = pointsVTK->InsertNextPoint(xx, yy, zz);
      }

    }

    uGridVTK->SetPoints(pointsVTK);


    cellDataVTKmatltype->SetName("matlType");
    cellDataVTKmatltype->SetNumberOfTuples(nElem_global);
    cellDataVTKelemtype->SetName("elemType");
    cellDataVTKelemtype->SetNumberOfTuples(nElem_global);

    for(ee=0;ee<nElem_global;ee++)
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
      cellDataVTKBresi->SetNumberOfTuples(nElem_global);

      cellDataVTKBappl->SetName("Bappl");
      cellDataVTKBappl->SetNumberOfComponents(3);
      cellDataVTKBappl->SetNumberOfTuples(nElem_global);

      VectorXd  vecTemp(3);

      for(ee=0;ee<nElem_global;ee++)
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
*/
    cout <<  " femMagnetomech::plotGeom ... ENDED" <<  endl;

    return;
}




void  femMagnetomech::postProcess()
{
    if(debug) cout << " femSolidmechanics::postProcess ... STARTED " << endl;

    if( (filecount % outputfreq_vtk) != 0 )
    {
      filecount++;
      return;
    }

    bool extrapolateFlag = false;
    int  vartype = 4;
    int  varindex = 8;
    int  index = 1;


    int  dd, ii, jj, kk, ll, nlocal, ind1, ind2, e, ee, count, gcount, ind;
    int  n1, n2, n3, n4, n5, n6, npElem;
    double vec[3], vec2[3], xx, yy, zz, val, geom[3];

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


    vtkIdType pt[4];
    vector<int>  nodeNums;

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

    vtkSmartPointer<vtkFloatArray>           presVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           mpotVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTK   =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKmatltype  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKelemtype  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKBappl  =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           cellDataVTKBresi  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkFloatArray>           timeloadstamp  =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();

    dispVTK->SetName("disp");
    dispVTK->SetNumberOfComponents(3);
    dispVTK->SetNumberOfTuples(nNode_global);

    veloVTK->SetName("velocity");
    veloVTK->SetNumberOfComponents(3);
    veloVTK->SetNumberOfTuples(nNode_global);

    mpotVTK->SetName("potential");
    mpotVTK->SetNumberOfTuples(nNode_global);

    //cellDataVTK->SetName("stress");
    //cellDataVTK->SetNumberOfTuples(nElem_global);

    presVTK->SetName("pressure");
    presVTK->SetNumberOfTuples(nNode_global);

    SxxNodesVTK->SetName("Sxx");
    SxxNodesVTK->SetNumberOfTuples(nNode_global);
    SyyNodesVTK->SetName("Syy");
    SyyNodesVTK->SetNumberOfTuples(nNode_global);
    SzzNodesVTK->SetName("Szz");
    SzzNodesVTK->SetNumberOfTuples(nNode_global);

    timeloadstamp->SetName("TIME");
    timeloadstamp->SetNumberOfTuples(1);
    if(ARC_LENGTH)
      timeloadstamp->SetTuple1(0, loadFactor);
    else
      timeloadstamp->SetTuple1(0, myTime.cur);


    // prepare points and cells
    if(ndim == 2)
    {
      vec[2] = 0.0;

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

    //cout << "pressure "  << endl;

    // pressure
    ////////////////////////////

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetFieldData()->AddArray(timeloadstamp);
    uGridVTK->GetPointData()->SetVectors(dispVTK);
    uGridVTK->GetPointData()->AddArray(veloVTK);

    VectorXd  output1(nNode_global), output2(nNode_global), output3(nNode_global), output4(nNode_global);

    //projectStresses(extrapolateFlag, 4,  0, 1, output1);
    //projectStresses(extrapolateFlag, 4,  4, 1, output2);
    //projectStresses(extrapolateFlag, 4,  8, 1, output3);
    //projectStresses(extrapolateFlag, 4, 10, 1, output4);

    for(ii=0; ii<nNode_global; ii++)
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
        val = 0.0;
        for(ii=0; ii<nNode_global; ii++)
        {
          //val = SolnData.var3[ii];

          if( midnodeData[ii][0] )
          {
            //val  = 0.50*SolnData.mpot[ii];
            //val += 0.25*SolnData.var3[midnodeData[ii][1]];
            //val += 0.25*SolnData.var3[midnodeData[ii][2]];
          }

          mpotVTK->SetTuple1(ii, val);
        }
      }
      else if(mpotDegree == -1)
      {
        val = 0.0;
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
    cellDataVTK->SetNumberOfTuples(nElem_global);

    cellDataVTKmatltype->SetName("matlType");
    cellDataVTKmatltype->SetNumberOfTuples(nElem_global);
    cellDataVTKelemtype->SetName("elemType");
    cellDataVTKelemtype->SetNumberOfTuples(nElem_global);

    for(ee=0;ee<nElem_global;ee++)
    {
        //elems[ee]->elementContourplot(4, 10, 1);

        //cellDataVTK->InsertTuple1(ee, elems[ee]->vals2project[0]);
      cellDataVTK->InsertTuple1(ee, 0.0);

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
      cellDataVTKBresi->SetNumberOfTuples(nElem_global);

      cellDataVTKBappl->SetName("Bappl");
      cellDataVTKBappl->SetNumberOfComponents(3);
      cellDataVTKBappl->SetNumberOfTuples(nElem_global);

      VectorXd  vecTemp;

      double  fact = timeFunctions[0]->getValue();

      for(ee=0;ee<nElem_global;ee++)
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
    sprintf(fname,"%s%s%06d%s", inputfilename.c_str(),"-",filecount, ".vtu");
    filecount++;

    writerUGridVTK->SetFileName(fname);
    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    if(debug) cout << " femSolidmechanics::postProcess ... ENDED " << endl;

    return;
}


void  femMagnetomech::projectStrains(bool extrapolateFlag, int vartype, int varindex, int index, VectorXd&  output)
{
    return;
}



void  femMagnetomech::projectStresses(bool extrapolateFlag, int vartype, int varindex, int index, VectorXd&  output)
{
    int  ee, ii;

    //cout << " extrapolateFlag = " << extrapolateFlag << endl;

    // compute the stresses at element nodes
    for(ee=0; ee<nElem_global; ee++)
    {
      elems[ee]->projectToNodes(extrapolateFlag, vartype, varindex, index);
    }

    // add element nodal contributions to the global list of nodes
    //if( stressVTK->GetNumberOfTuples() != nNode_global );
    //{
      //stressVTK->SetNumberOfComponents(9);
      //stressVTK->SetNumberOfTuples(nNode_global);
    //}

    vector<int>  nodeNums;
    double  volu, val;

    if(output.rows() != nNode_global)  output.resize(nNode_global);
    output.setZero();

    for(ee=0; ee<nElem_global; ee++)
    {
      nodeNums = elems[ee]->nodeNums;

      elems[ee]->computeVolume(true);

      volu = elems[ee]->getVolume();

      for(ii=0; ii<nodeNums.size(); ii++)
      {
        output[nodeNums[ii]] += volu*elems[ee]->vals2project[ii];
      }
    }

    for(ee=0; ee<nNode_global; ee++)
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
    for(ii=0; ii<nNode_global; ii++)
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


void  femMagnetomech::projectInternalVariables(bool extrapolateFlag, int vartype, int varindex, int index, VectorXd&  output)
{
    return;
}





void  femMagnetomech::projectFromElemsToNodes(bool extrapolateFlag, int vartype, int varindex, int index)
{
    int  ee, ii;


    return;
}





void  femMagnetomech::writeResult()
{
    string  fname = inputfilename+".soln";

    ofstream fout(fname);

    if(fout.fail())
    {
       cout << " Could not open the Output file for writing the solution" << endl;
       exit(1);
    }

    int nn, ind;

    fout << "Dimension "  << ndim << endl;
    fout << "Time "       << myTime.cur << endl;
    fout << "Nodes "      << nNode_global << endl;
    fout << "Elments "    << nElem_global << endl;

    fout << "Displacement" << endl;
    if(ndim == 3)
    {
      for(nn=0; nn<nNode_global; nn++)
      {
        ind = nn*ndof;
        fout << nn+1 << setw(20) << SolnData.disp(ind) << setw(20) << SolnData.disp(ind+1) << setw(20) << SolnData.disp(ind+2) << endl;
      }
    }
    else
    {
      for(nn=0; nn<nNode_global; nn++)
      {
        ind = nn*ndof;
        fout << nn+1 << setw(20) << SolnData.disp(ind) << setw(20) << SolnData.disp(ind+1) << endl;
      }
    }

    fout << "Velocity" << endl;
    if(ndim == 3)
    {
      for(nn=0; nn<nNode_global; nn++)
      {
        ind = nn*ndof;
        fout << nn+1 << setw(20) << SolnData.velo(ind) << setw(20) << SolnData.velo(ind+1) << setw(20) << SolnData.velo(ind+2) << endl;
      }
    }
    else
    {
      for(nn=0; nn<nNode_global; nn++)
      {
        ind = nn*ndof;
        fout << nn+1 << setw(20) << SolnData.velo(ind) << setw(20) << SolnData.velo(ind+1) << endl;
      }
    }

    fout << "Acceleration" << endl;
    if(ndim == 3)
    {
      for(nn=0; nn<nNode_global; nn++)
      {
        ind = nn*ndof;
        fout << nn+1 << setw(20) << SolnData.acce(ind) << setw(20) << SolnData.acce(ind+1) << setw(20) << SolnData.acce(ind+2) << endl;
      }
    }
    else
    {
      for(nn=0; nn<nNode_global; nn++)
      {
        ind = nn*ndof;
        fout << nn+1 << setw(20) << SolnData.acce(ind) << setw(20) << SolnData.acce(ind+1) << endl;
      }
    }

    fout << "Pressure" << endl;
    for(nn=0; nn<nNode_global; nn++)
    {
        fout << nn+1 << setw(20) << SolnData.pres(nn) << endl;
    }

    fout.close();

    return;
}





void  femMagnetomech::readResult()
{
    if(restartfilename.size() == 0)
      return;

    cout << "restartfilename  = " << restartfilename << endl;

    string  fname  = "./inputs/"+restartfilename;


    ifstream infile(fname);

    if(infile.fail())
    {
       cout << " Could not open the input mesh file" << endl;
       exit(1);
    }

    string line;
    vector<string>  stringlist;
    int  nn, dd, ind;


    //read the dimension
    getline(infile,line);    boost::trim(line);
    //cout << line << endl;
    boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);
    assert(stoi(stringlist[1]) == ndim);

    //read time
    getline(infile,line);    boost::trim(line);
    //cout << line << endl;

    //read node count
    getline(infile,line);
    //cout << line << endl;
    boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);
    assert(stoi(stringlist[1]) == nNode_global);

    //read element count
    getline(infile,line);
    cout << line << endl;
    boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);
    assert(stoi(stringlist[1]) == nElem_global);

    //read Displacement tag
    getline(infile,line);

    //read Displacement values
    for(nn=0; nn<nNode_global; nn++)
    {
        getline(infile,line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        nn = stoi(stringlist[0])-1;
        ind = nn*ndof;

        for(dd=0; dd<ndim; dd++)
          SolnData.disp[ind+dd]   = stod(stringlist[dd+1]);
    }

    //read Velocity tag
    getline(infile,line);

    //read Displacement values
    for(nn=0; nn<nNode_global; nn++)
    {
        getline(infile,line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        nn = stoi(stringlist[0])-1;
        ind = nn*ndof;

        for(dd=0; dd<ndim; dd++)
          SolnData.velo[ind+dd]   = stod(stringlist[dd+1]);
    }

    //read Acceleration tag
    getline(infile,line);

    //read Displacement values
    for(nn=0; nn<nNode_global; nn++)
    {
        getline(infile,line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        nn = stoi(stringlist[0])-1;
        ind = nn*ndof;

        for(dd=0; dd<ndim; dd++)
          SolnData.acce[ind+dd]   = stod(stringlist[dd+1]);
    }

    //read Pressure tag
    getline(infile,line);

    //read Displacement values
    for(nn=0; nn<nNode_global; nn++)
    {
        getline(infile,line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        nn = stoi(stringlist[0])-1;

        SolnData.pres[nn] = stod(stringlist[1]);
    }

    infile.close();


    SolnData.saveSolution();

    SolnData.dispCur = SolnData.disp;

    for(nn=0; nn<nNode_global; nn++)
    {
      for(dd=0; dd<ndim; dd++)
      {
        // ndof is used here instead of ndim.
        // this is important for beam and shell elements
        ind = nn*ndof+dd;

        GeomData.NodePosNew[nn][dd] = GeomData.NodePosOrig[nn][dd] + SolnData.disp[ind];
        GeomData.NodePosCur[nn][dd] = GeomData.NodePosOrig[nn][dd] + SolnData.dispCur[ind];
      }
    }

    return;
}






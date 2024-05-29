
#include "femINSmixed.h"
#include "headersVTK.h"
#include "elementutilitiescfd.h"
#include "MyTime.h"
#include "TimeFunction.h"

extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime                 myTime;



int  femINSmixed::findElementNumber(myPoint& target_point)
{
    int elnum = -1;

    if( (npElemVelo == 6) || (npElemVelo == 7) ) // triangular elements
    {
      int nn;
      double  xNode[3], yNode[3];
      for(int ee=0; ee<nElem_global; ee++)
      {
        //printVector(elems[ee]->nodeNums);
        // get nodal coordinates
        nn = elems[ee]->nodeNums[0];
        xNode[0] = nodeCoords[nn][0];    yNode[0] = nodeCoords[nn][1];
        nn = elems[ee]->nodeNums[1];
        xNode[1] = nodeCoords[nn][0];    yNode[1] = nodeCoords[nn][1];
        nn = elems[ee]->nodeNums[2];
        xNode[2] = nodeCoords[nn][0];    yNode[2] = nodeCoords[nn][1];
        
        //cout << xNode[0] << '\t' << yNode[0] << endl;
        //cout << xNode[1] << '\t' << yNode[1] << endl;
        //cout << xNode[2] << '\t' << yNode[2] << endl;

        //cout << "pointInsideTria3Node    " << ee << endl;
        if( pointInsideTria3Node(xNode, yNode, &target_point[0]) )
        {
          //cout << "Element number = " <<  ee << endl;
          elnum = ee;

          break;
        }
      }
    }
    else //quadrilateral element
    {
      for(int ee=0; ee<nElem_global; ++ee)
      {
        if( elems[ee]->bbox.within(target_point) )
        {
          cout << "Element number = " <<  ee << endl;
          //ImmersedElements[ime]->elem = elems[ee];
          elnum = ee;

          break;
        }
      }
    }


    return elnum;
}





int  femINSmixed::postProcess()
{
    cout << " femINSmixed::postProcess " << endl;

    postProcessVelocity();
    cout << " postProcessVelocity " << endl;

    postProcessPressure();
    
    if(IB_MOVED)
    {
      for(int bb=0; bb<ImmersedBodyObjects.size(); bb++)
        ImmersedBodyObjects[bb]->postProcess(dirname, fileCount);
    }

    fileCount = fileCount+1;

    return 0;
}


int  femINSmixed::postProcessVelocity()
{
    cout << " femINSmixed::postProcessVelocity " << endl;

    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    //vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    //vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    //vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    //vtkSmartPointer<vtkHexahedron>           hexaVTK      =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK     =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK     =  vtkSmartPointer<vtkBiQuadraticQuad>::New();
    vtkSmartPointer<vtkQuadraticTetra>       tetra10VTK   =  vtkSmartPointer<vtkQuadraticTetra>::New();
    vtkSmartPointer<vtkQuadraticHexahedron>  hexa27VTK    =  vtkSmartPointer<vtkQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           vecVTK       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2;
    double  val, vec[3]={0.0,0.0,0.0};

    vecVTK->SetNumberOfTuples(nNode_Velo);
    vecVTK->SetNumberOfComponents(3);
    //scaVTK->SetNumberOfTuples(nNode);

    vecVTK->SetName("velocity");
    //scaVTK->SetName("pressure");

    vtkIdType pt[10];

    if(ndim == 2)
    {
      for(ii=0; ii<nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(nodeCoords[ii][0], nodeCoords[ii][1], 0.0);

        if( midNodeData[ii][0] )
        {
          jj = ii*ndim;
          n1 = midNodeData[ii][1]*ndim;
          n2 = midNodeData[ii][2]*ndim;

          vec[0] = 0.5*velo(jj)   + 0.25*velo(n1)   + 0.25*velo(n2);
          vec[1] = 0.5*velo(jj+1) + 0.25*velo(n1+1) + 0.25*velo(n2+1);

          //vec[0] = velo(jj)   ;
          //vec[1] = velo(jj+1) ;

          // pressure at the mid-nodes is computed as the average
          // only for the purpose of plotting
          //val = 0.5*pres(midNodeData[ii][1]) + 0.5*pres(midNodeData[ii][2]) ;

          vecVTK->InsertTuple(ii, vec);
          //scaVTK->SetTuple1(ii, val);
        }
        else
        {
          kk = ii*ndim;

          vec[0] = velo[kk];
          vec[1] = velo[kk+1];
          //val    = pres[ii];

          vecVTK->InsertTuple(ii, vec);
          //scaVTK->SetTuple1(ii, val);
        }
      }
      
      cout << " postProcess ... added points " << endl;

      for(ee=0; ee<nElem_global; ee++)
      {
        npElemVelo = elemConn[ee].size();

        if( (npElemVelo == 6) || (npElemVelo == 7) )
        {
          for(ii=0; ii<6; ii++)
            tria6VTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
        else if(npElemVelo == 9)
        {
          vec[0]  = 0.25*velo(elemConn[ee][8]*ndim);
          vec[0] += 0.0625*(velo(elemConn[ee][0]*ndim)+velo(elemConn[ee][1]*ndim)+velo(elemConn[ee][2]*ndim)+velo(elemConn[ee][3]*ndim));
          vec[0] += 0.1250*(velo(elemConn[ee][4]*ndim)+velo(elemConn[ee][5]*ndim)+velo(elemConn[ee][6]*ndim)+velo(elemConn[ee][7]*ndim));

          vec[1]  = 0.25*velo(elemConn[ee][8]*ndim+1);
          vec[1] += 0.0625*(velo(elemConn[ee][0]*ndim+1)+velo(elemConn[ee][1]*ndim+1)+velo(elemConn[ee][2]*ndim+1)+velo(elemConn[ee][3]*ndim+1));
          vec[1] += 0.1250*(velo(elemConn[ee][4]*ndim+1)+velo(elemConn[ee][5]*ndim+1)+velo(elemConn[ee][6]*ndim+1)+velo(elemConn[ee][7]*ndim+1));

          vecVTK->InsertTuple(elemConn[ee][8], vec);

          val = 0.25*(pres(elemConn[ee][0])+pres(elemConn[ee][1])+pres(elemConn[ee][2])+pres(elemConn[ee][3]));

          scaVTK->SetTuple1(elemConn[ee][8], val);

          for(ii=0; ii<npElemVelo; ii++)
            quad9VTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
      }
      cout << " postProcess ... added elements " << endl;
    }
    else // (ndim == 3)
    {
      for(ii=0; ii<nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(nodeCoords[ii][0], nodeCoords[ii][1], nodeCoords[ii][2]);

        if( midNodeData[ii][0] )
        {
          jj = ii*ndim;
          n1 = midNodeData[ii][1]*ndim;
          n2 = midNodeData[ii][2]*ndim;

          vec[0] = 0.5*velo(jj)   + 0.25*velo(n1)   + 0.25*velo(n2);
          vec[1] = 0.5*velo(jj+1) + 0.25*velo(n1+1) + 0.25*velo(n2+1);
          vec[2] = 0.5*velo(jj+2) + 0.25*velo(n1+2) + 0.25*velo(n2+2);

          // pressure at the mid-nodes is computed as the average
          // only for the purpose of plotting
          val   = 0.5*pres(midNodeData[ii][1]) + 0.5*pres(midNodeData[ii][1]) ;

          vecVTK->InsertTuple(ii, vec);
          scaVTK->SetTuple1(ii, val);
        }
        else
        {
          kk = ii*ndim;

          vec[0] = velo[kk];
          vec[1] = velo[kk+1];
          vec[2] = velo[kk+2];
          val    = pres[ii];

          vecVTK->InsertTuple(ii, vec);
          scaVTK->SetTuple1(ii, val);
        }
      }

      for(ee=0; ee<nElem_global; ee++)
      {
        npElemVelo = elemConn[ee].size();

        if(npElemVelo == 10)
        {
          for(ii=0; ii<npElemVelo; ii++)
            tetra10VTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

          uGridVTK->InsertNextCell(tetra10VTK->GetCellType(), tetra10VTK->GetPointIds());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    //uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->SetVectors(vecVTK);

    //Write the file.

    char VTKfilename[400];

    //sprintf(VTKfilename,"%s%s%06d%s", infilename.c_str(), "-",fileCount,".vtu");
    sprintf(VTKfilename,"%s%s%06d%s", dirname.c_str(), "-velocity-", fileCount,".vtu");

    writerUGridVTK->SetFileName(VTKfilename);

    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();
    cout << "iiiiiiiiiiiiiiiiiii " << endl;

    return 0;
}






int  femINSmixed::postProcessPressure()
{
    cout << " femINSmixed::postProcessPressure " << endl;

    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK     =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK    =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK    =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkTriangle>             triaVTK      =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK      =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK     =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkHexahedron>           hexaVTK      =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK     =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK     =  vtkSmartPointer<vtkBiQuadraticQuad>::New();
    vtkSmartPointer<vtkQuadraticTetra>       tetra10VTK   =  vtkSmartPointer<vtkQuadraticTetra>::New();
    vtkSmartPointer<vtkQuadraticHexahedron>  hexa27VTK    =  vtkSmartPointer<vtkQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           scaVTK       =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2;
    double  val, vec[3]={0.0,0.0,0.0};

    if(npElemVelo == 7)
      scaVTK->SetNumberOfTuples(nNode_Pres);
    else
      scaVTK->SetNumberOfTuples(nNode_Velo);

    scaVTK->SetName("pressure");

    vtkIdType pt[10];

    if(ndim == 2)
    {
      if(npElemVelo == 7)
      {
        for(ee=0; ee<nElem_global; ee++)
        {
          for(ii=0; ii<npElemPres; ii++)
          {
            nn = elems[ee]->nodeNums[ii];

            pt[0] = pointsVTK->InsertNextPoint(nodeCoords[nn][0], nodeCoords[nn][1], 0.0);

            triaVTK->GetPointIds()->SetId(ii, pt[0] );

            val    = pres[elems[ee]->nodeNumsPres[ii]];

            scaVTK->SetTuple1(pt[0], val);
          }

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
      }
      else
      {
        for(ii=0; ii<nNode_global; ii++)
        {
          pt[0] = pointsVTK->InsertNextPoint(nodeCoords[ii][0], nodeCoords[ii][1], 0.0);

          if( midNodeData[ii][0] )
          {
            // pressure at the mid-nodes is computed as the average
            // only for the purpose of plotting
            val = 0.5*pres(midNodeData[ii][1]) + 0.5*pres(midNodeData[ii][2]) ;

            scaVTK->SetTuple1(ii, val);
          }
          else
          {
            kk = ii*ndim;

            val    = pres[ii];

            scaVTK->SetTuple1(ii, val);
          }
        }

        for(ee=0; ee<nElem_global; ee++)
        {
          npElemVelo = elemConn[ee].size();

          if(npElemVelo == 6)
          {
            for(ii=0; ii<npElemVelo; ii++)
              tria6VTK->GetPointIds()->SetId(ii, elemConn[ee][ii] );

            uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
          }
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);

    //Write the file.

    char VTKfilename[200];

    sprintf(VTKfilename,"%s%s%06d%s", dirname.c_str(), "-pressure-",fileCount,".vtu");

    writerUGridVTK->SetFileName(VTKfilename);

    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    return 0;
}








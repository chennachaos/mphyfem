
#include "headersVTK.h"
#include "femNSCHstabilised.h"


int  femNSCHstabilised::postProcess()
{
    if(this_mpi_proc != 0)
      return 0;

    PetscPrintf(MPI_COMM_WORLD, " \n femNSCHstabilised::postProcess \n");

    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK      =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK     =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK     =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkTriangle>             triaVTK       =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK       =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK      =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkHexahedron>           hexaVTK       =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK      =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK      =  vtkSmartPointer<vtkBiQuadraticQuad>::New();
    vtkSmartPointer<vtkQuadraticTetra>       tetra10VTK    =  vtkSmartPointer<vtkQuadraticTetra>::New();
    vtkSmartPointer<vtkQuadraticHexahedron>  hexa27VTK     =  vtkSmartPointer<vtkQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           vecVTK        =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK        =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK2       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK3       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK4       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           procIdVTKnode =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           procIdVTKcell =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2, npElem;
    double  val, vec[3]={0.0,0.0,0.0}, *solnNS = &(NavierStokes.soln[0]), *solnCH = &(CahnHilliard.soln[0]);
    vector<int> nodeNums;

    vecVTK->SetName("velocity");
    vecVTK->SetNumberOfTuples(mesh->nNode_global);
    vecVTK->SetNumberOfComponents(3);

    scaVTK->SetName("pressure");
    scaVTK->SetNumberOfTuples(mesh->nNode_global);

    scaVTK2->SetName("phi");
    scaVTK2->SetNumberOfTuples(mesh->nNode_global);

    scaVTK3->SetName("eta");
    scaVTK3->SetNumberOfTuples(mesh->nNode_global);

    scaVTK4->SetName("rho");
    scaVTK4->SetNumberOfTuples(mesh->nNode_global);

    procIdVTKcell->SetName("procId");
    procIdVTKcell->SetNumberOfTuples(mesh->nElem_global);

    procIdVTKnode->SetName("procId");
    procIdVTKnode->SetNumberOfTuples(mesh->nNode_global);

    vtkIdType pt[10];

    if(mesh->ndim == 2)
    {
      for(ii=0; ii<mesh->nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(mesh->nodeCoordsOrig[ii][0], mesh->nodeCoordsOrig[ii][1], 0.0);

        procIdVTKnode->InsertTuple1(ii, mesh->node_proc_id[ii]);

        nn = mesh->node_map_get_new[ii];
        kk = nn*(ndim+1);

        // NavierStokes solution
        vec[0] = solnNS[kk];
        vec[1] = solnNS[kk+1];
        val    = solnNS[kk+2];

        vecVTK->InsertTuple(ii, vec);
        scaVTK->SetTuple1(ii, val);


        // NavierStokes solution
        scaVTK2->SetTuple1(ii, solnCH[nn*2]);
        scaVTK3->SetTuple1(ii, solnCH[nn*2+1]);

        scaVTK4->SetTuple1(ii, NavierStokes.rhoNodal[nn]);
      }

      for(ee=0; ee<mesh->nElem_global; ee++)
      {
        procIdVTKcell->SetTuple1(ee, mesh->elem_proc_id[ee]);

        nodeNums = mesh->elemConn[ee];
        npElem   = nodeNums.size();

        if(npElem == 3)
        {
          for(ii=0; ii<npElem; ii++)
            triaVTK->GetPointIds()->SetId(ii, mesh->node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
        else if(npElem == 6)
        {
          for(ii=0; ii<npElem; ii++)
            tria6VTK->GetPointIds()->SetId(ii, mesh->node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
        else if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            quadVTK->GetPointIds()->SetId(ii, mesh->node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
        else if(npElem == 9)
        {
          for(ii=0; ii<npElem; ii++)
            quad9VTK->GetPointIds()->SetId(ii, mesh->node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
      }
    }
    else // (ndim == 3)
    {
      for(ii=0; ii<mesh->nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(mesh->nodeCoordsOrig[ii][0], mesh->nodeCoordsOrig[ii][1], mesh->nodeCoordsOrig[ii][2]);

        procIdVTKnode->InsertTuple1(ii, mesh->node_proc_id[ii]);

        nn = mesh->node_map_get_new[ii];
        kk = nn*(ndim+1);

        vec[0] = solnNS[kk];
        vec[1] = solnNS[kk+1];
        vec[2] = solnNS[kk+2];
        val    = solnNS[kk+3];

        vecVTK->InsertTuple(ii, vec);
        scaVTK->SetTuple1(ii, val);
      }

      for(ee=0; ee<mesh->nElem_global; ee++)
      {
        procIdVTKcell->SetTuple1(ee, mesh->elem_proc_id[ee]);

        nodeNums = mesh->elemConn[ee];
        npElem   = nodeNums.size();

        if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            tetraVTK->GetPointIds()->SetId(ii, mesh->node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
        else if(npElem == 8)
        {
          for(ii=0; ii<npElem; ii++)
            hexaVTK->GetPointIds()->SetId(ii, mesh->node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(hexaVTK->GetCellType(), hexaVTK->GetPointIds());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->AddArray(procIdVTKnode);
    uGridVTK->GetPointData()->AddArray(scaVTK2);
    uGridVTK->GetPointData()->AddArray(scaVTK3);
    uGridVTK->GetPointData()->AddArray(scaVTK4);
    uGridVTK->GetPointData()->SetVectors(vecVTK);

    uGridVTK->GetCellData()->SetScalars(procIdVTKcell);

    //Write the file.

    char VTKfilename[200];

    cout << "infilename = " << infilename << endl;
    sprintf(VTKfilename,"%s%s%06d%s", infilename.c_str(), "-",fileCount,".vtu");

    writerUGridVTK->SetFileName(VTKfilename);

    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    fileCount++;

    PetscPrintf(MPI_COMM_WORLD, " \n femNSCHstabilised::postProcess \n");

    return 0;
}









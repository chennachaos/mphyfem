
#include "headersVTK.h"
#include "myMesh.h"
#include "util.h"
#include <chrono>
#include "metis.h"
#include "MyTime.h"
#include "TimeFunction.h"
//#include "metis.h"


extern  std::vector<unique_ptr<TimeFunction> > timeFunction;
extern  MyTime                 myTime;





int  myMesh::postProcess()
{
    if(this_mpi_proc != 0)
      return 0;

    PetscPrintf(MPI_COMM_WORLD, " \n myMesh::postProcess \n");
/*
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
    vtkSmartPointer<vtkFloatArray>           procIdVTKnode =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           procIdVTKcell =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2;
    double  val, vec[3]={0.0,0.0,0.0};
    vector<int> nodeNums;

    vecVTK->SetName("velocity");
    vecVTK->SetNumberOfTuples(nNode_global);
    vecVTK->SetNumberOfComponents(3);

    scaVTK->SetName("pressure");
    scaVTK->SetNumberOfTuples(nNode_global);

    procIdVTKcell->SetName("procId");
    procIdVTKcell->SetNumberOfTuples(nElem_global);

    procIdVTKnode->SetName("procId");
    procIdVTKnode->SetNumberOfTuples(nNode_global);

    vtkIdType pt[10];

    if(ndim == 2)
    {
      for(ii=0; ii<nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(nodeCoordsOrig[ii][0], nodeCoordsOrig[ii][1], 0.0);

        procIdVTKnode->InsertTuple1(ii, node_proc_id[ii]);

        nn = node_map_get_new[ii];
        kk = nn*ndof;

        vec[0] = soln[kk];
        vec[1] = soln[kk+1];
        val    = soln[kk+2];

        vecVTK->InsertTuple(ii, vec);
        scaVTK->SetTuple1(ii, val);
      }

      for(ee=0; ee<nElem_global; ee++)
      {
        procIdVTKcell->SetTuple1(ee, elem_proc_id[ee]);

        nodeNums = elemConn[ee];
        npElem   = nodeNums.size();

        if(npElem == 3)
        {
          for(ii=0; ii<npElem; ii++)
            triaVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
        else if(npElem == 6)
        {
          for(ii=0; ii<npElem; ii++)
            tria6VTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
        else if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            quadVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
        else if(npElem == 9)
        {
          for(ii=0; ii<npElem; ii++)
            quad9VTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
      }
    }
    else // (ndim == 3)
    {
      for(ii=0; ii<nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(nodeCoordsOrig[ii][0], nodeCoordsOrig[ii][1], nodeCoordsOrig[ii][2]);

        procIdVTKnode->InsertTuple1(ii, node_proc_id[ii]);

        nn = node_map_get_new[ii];
        kk = nn*ndof;

        vec[0] = soln[kk];
        vec[1] = soln[kk+1];
        vec[2] = soln[kk+2];
        val    = soln[kk+3];

        vecVTK->InsertTuple(ii, vec);
        scaVTK->SetTuple1(ii, val);
      }

      for(ee=0; ee<nElem_global; ee++)
      {
        procIdVTKcell->SetTuple1(ee, elem_proc_id[ee]);

        nodeNums = elemConn[ee];
        npElem   = nodeNums.size();

        if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            tetraVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
        else if(npElem == 8)
        {
          for(ii=0; ii<npElem; ii++)
            hexaVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(hexaVTK->GetCellType(), hexaVTK->GetPointIds());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->AddArray(procIdVTKnode);
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
*/
    PetscPrintf(MPI_COMM_WORLD, " \n myMesh::postProcess \n");

    return 0;
}









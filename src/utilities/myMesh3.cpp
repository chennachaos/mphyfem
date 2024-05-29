#include "headersVTK.h"
#include "myMesh.h"
#include "util.h"
#include <chrono>
#include "metis.h"




int myMesh::partitionMesh()
{
    PetscPrintf(MPI_COMM_WORLD, "\n     myMesh::partitionmyMesh()  .... STARTED ...\n");

    int  ee, ii, jj, kk, n1, n2;
    int  nparts = n_mpi_procs, subdomain=0;

    if(this_mpi_proc == 0)
    {

        /////////////////////////////////////////////////////////////////////////////
        //
        // Partition the mesh. Here, METIS is used.
        // 
        /////////////////////////////////////////////////////////////////////////////

        PetscInt  *eptr, *eind;

        errpetsc  = PetscMalloc1(nElem_global+1,  &eptr);CHKERRQ(errpetsc);
        errpetsc  = PetscMalloc1(nElem_global*npElem,  &eind);CHKERRQ(errpetsc);

        vector<int>  vecTemp2;

        eptr[0] = 0;
        kk = 0;
        for(ee=0; ee<nElem_global; ee++)
        {
            eptr[ee+1] = (ee+1)*npElem;

            vecTemp2 = elemConn[ee];
            //printVector(vecTemp2);

            for(ii=0; ii<npElem; ii++)
              eind[kk+ii] = vecTemp2[ii] ;

            kk += npElem;
        }

        int  ncommon_nodes;
        if(ndim == 2)
        {
          if( (npElem == 3) || (npElem == 4) )       // 3-noded tria or 4-noded quad elements
            ncommon_nodes = 2;
          else
            ncommon_nodes = 3;  // 6-noded tria or 9-noded quad elements
        }
        else
        {
          if(npElem == 4)          // 4-noded tetra element
            ncommon_nodes = 3;
          else if(npElem == 10)    // 10-noded tetra element
            ncommon_nodes = 6;
          else if(npElem == 6)     // 6-noded Wedge/Prism element
            ncommon_nodes = 4;
          else if(npElem == 8)     // 8-noded hexa element
            ncommon_nodes = 4;
          else if(npElem == 18)    // 18-noded Wedge/Prism element
            ncommon_nodes = 9;
          else
            ncommon_nodes = 9;     // 27-noded hexa element
        }

        idx_t objval;
        idx_t options[METIS_NOPTIONS];

        METIS_SetDefaultOptions(options);

        // Specifies the partitioning method.
        //options[METIS_OPTION_PTYPE] = METIS_PTYPE_RB;          // Multilevel recursive bisectioning.
        options[METIS_OPTION_PTYPE] = METIS_PTYPE_KWAY;        // Multilevel k-way partitioning.

        //options[METIS_OPTION_NSEPS] = 10;

        //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;     // Edge-cut minimization
        options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_VOL;     // Total communication volume minimization

        options[METIS_OPTION_NUMBERING] = 0;  // C-style numbering is assumed that starts from 0.

        cout << " Executing METIS subroutine " << endl;

        // METIS partition routine
        int ret = METIS_PartMeshNodal(&nElem_global, &nNode_global, eptr, eind, NULL, NULL, &nparts, NULL, options, &objval, &elem_proc_id[0], &node_proc_id[0]);
        //int ret = METIS_PartMeshDual(&nElem_global, &nNode_global, eptr, eind, NULL, NULL, &ncommon_nodes, &nparts, NULL, options, &objval, &elem_proc_id[0], &node_proc_id[0]);

        if(ret == METIS_OK)
          cout << " METIS partition routine successful "  << endl;
        else
          cout << " METIS partition routine FAILED "  << endl;

        errpetsc = PetscFree(eptr); CHKERRQ(errpetsc);
        errpetsc = PetscFree(eind); CHKERRQ(errpetsc);

        if( 1 < 0)
        {
          for(ee=0; ee<nNode_global; ee++)
            cout << ee << '\t' << node_proc_id[ee] << endl;
          cout << endl;  cout << endl;  cout << endl;

          for(ee=0; ee<nElem_global; ee++)
            cout << ee << '\t' << elem_proc_id[ee] << endl;
          cout << endl;  cout << endl;  cout << endl;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    errpetsc = MPI_Bcast(&elem_proc_id[0], nElem_global, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
    MPI_Barrier(MPI_COMM_WORLD);

    errpetsc = MPI_Bcast(&node_proc_id[0], nNode_global, MPI_INT, 0, MPI_COMM_WORLD); CHKERRQ(errpetsc);
    MPI_Barrier(MPI_COMM_WORLD);

    nElem_local = std::count(elem_proc_id.begin(), elem_proc_id.end(), this_mpi_proc);
    cout << " nElem_local =  " << nElem_local << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;


      nNode_owned = std::count(node_proc_id.begin(), node_proc_id.end(), this_mpi_proc);
      cout << " nNode_owned =  " << nNode_owned << '\t' << this_mpi_proc << '\t' << n_mpi_procs << endl;

      MPI_Barrier(MPI_COMM_WORLD);

      // find nodes local to each processor generate the list of locally owned nodes
      std::vector<int>  nodelist_owned(nNode_owned);

      kk=0;
      for(ii=0; ii<nNode_global; ii++)
      {
        if( node_proc_id[ii] == this_mpi_proc )
        {
          nodelist_owned[kk++] = ii;
        }
      }
      //cout << " Locally owned nodes " << '\t' << this_mpi_proc << endl;
      //printVector(nodelist_owned);

      MPI_Barrier(MPI_COMM_WORLD);

      // create the vector (of size n_mpi_procs)
      // consisting of nNode_owned from all the processors in the communication
      vector<int>  nNode_owned_vector(n_mpi_procs), nNode_owned_sum(n_mpi_procs);

      MPI_Allgather(&nNode_owned, 1, MPI_INT, &nNode_owned_vector[0], 1, MPI_INT, MPI_COMM_WORLD);
      //printVector(nNode_owned_vector);

      // compute the numbers of first and last nodes in the local processor
      nNode_owned_sum = nNode_owned_vector;
      for(ii=1; ii<n_mpi_procs; ii++)
      {
        nNode_owned_sum[ii] += nNode_owned_sum[ii-1];
      }
      //printVector(nNode_owned_sum);
      node_start = 0;
      if(this_mpi_proc > 0)
        node_start = nNode_owned_sum[this_mpi_proc-1];
      node_end   = nNode_owned_sum[this_mpi_proc]-1;

      cout << " node_start =  " << node_start << '\t' << node_end << '\t' << this_mpi_proc << endl;

      MPI_Barrier(MPI_COMM_WORLD);

      std::vector<int>  displs(n_mpi_procs);

      displs[0] = 0;
      for(ii=0; ii<n_mpi_procs-1; ii++)
        displs[ii+1] = displs[ii] + nNode_owned_vector[ii];

      // create a global list of nodelist_owned
      // which will serve as a mapping from NEW node numbers to OLD node numbers
      errpetsc = MPI_Allgatherv(&nodelist_owned[0], nNode_owned, MPI_INT, &node_map_get_old[0], &nNode_owned_vector[0], &displs[0], MPI_INT, MPI_COMM_WORLD);


      // create an array for mapping from OLD node numbers to NEW node numbers
      // Also, generate NodeTypeNew array for computing the local and global DOF size
      // as well as creating the element-wise array for element matrix/vector assembly
      for(ii=0; ii<nNode_global; ii++)
      {
        n1 = node_map_get_old[ii];
        node_map_get_new[n1] = ii;
      }

      // update elem<->node connectivity with new node numbers and delete arrays in other processors
      //vector<int>  vecTemp;
      for(ee=0; ee<nElem_global; ee++)
      {
          //if( (elem_proc_id[ee] == 0) || (elem_proc_id[ee] == this_mpi_proc) )
          //{
            for(ii=0; ii<npElem; ii++)
            {
              //vecTemp.push_back(elemConn[ee][ii]);
              elemConn[ee][ii] = node_map_get_new[elemConn[ee][ii]];
            }
          //}

          //if( !(elem_proc_id[ee] == 0) && !(elem_proc_id[ee] == this_mpi_proc) )
            //elemConn[ee].clear();
      }
      MPI_Barrier(MPI_COMM_WORLD);


      // update node numbers for boundary patches
      // update Dirichlet BC information with new node numbers
      //
      for(ii=0; ii<numBoundaryPatches; ++ii)
      {
        BoundaryPatches[ii]->updateNodeNumbers(node_map_get_new);
      }
      MPI_Barrier(MPI_COMM_WORLD);


    PetscPrintf(MPI_COMM_WORLD, "\n     myMesh::partitionmyMesh()  .... FINISHED ...\n");

    return 0;
}





int myMesh::prepareDataForParallel()
{

  return 0;
}






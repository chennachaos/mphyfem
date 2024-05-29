
#ifndef incl_myMesh_h
#define incl_myMesh_h

#include <string.h>
#include <vector>
#include <fstream>
#include <memory>
#include "BoundaryPatch.h"
#include "util.h"
#include "myMathFunction.h"
#include "petscksp.h"
#include "petscmat.h"


using namespace std;



class myMesh
{
    //private:

    public:

        PetscErrorCode  errpetsc;
        PetscInt  n_mpi_procs, this_mpi_proc;
        PetscInt  ndim, npElem, numBoundaryPatches;
        PetscInt  nElem_global, nNode_global, nElem_local, nNode_local, nNode_owned;
        PetscInt  node_start, node_end, elem_start, elem_end;
        
        bool  MESH_PREPARED;

        vector<int>              elem_proc_id, node_proc_id, node_map_get_old, node_map_get_new;
        vector<myPoint>          nodeCoordsOrig, nodeCoordsCur;               //!< coordinates of the nodes (or control points)
        vector<vector<int> >     elemConn;                  //!< element-node connectivity array

        vector<unique_ptr<BoundaryPatch> >    BoundaryPatches;

    public:

        myMesh();

        ~myMesh();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int getNDIM()
        {  return ndim;}

        int getnNodeLocal()
        {  return nNode_local;}

        int getnNodeGlobal()
        {  return nNode_global;}

        int getnElemLocal()
        {  return nElem_local;}

        int getnElemGlobal()
        {  return nElem_global;}

        int  readInputGMSH(string& fname);

        int  readBoundariesData(ifstream& infile, string& line);

        int  prepareInputData();

        int  readInputData(string& fname);

        int  partitionMesh();

        int  prepareDataForParallel();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  writeOutputData();

        int  postProcess();
};






#endif







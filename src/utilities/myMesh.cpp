
#include <algorithm>
#include <chrono>
#include "myMesh.h"
#include "util.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include "MyTime.h"
#include "TimeFunction.h"
#include "FunctionsProgram.h"
#include "mymapkeys.h"


extern  std::vector<unique_ptr<TimeFunction> > timeFunction;
extern  MyTime                 myTime;
extern  bool  debug;

using namespace std;


myMesh::myMesh()
{
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    nElem_global = 0; nNode_global = 0; npElem = 0;
    MESH_PREPARED = false;
}


myMesh::~myMesh()
{
    //cout << " myMesh::~myMesh() " << this_mpi_proc << endl;

    //deallocatePetscObjects();

    //cout << " myMesh::~myMesh() " << this_mpi_proc << endl;
}





int  myMesh::readInputGMSH(string& fname)
{
    PetscPrintf(MPI_COMM_WORLD, " myMesh::readInputGMSH \n");

    ifstream  infile(fname);

    if(infile.fail())
    {
       cout << " Could not open input file " << endl;
       exit(-1);
    }

    string line;
    vector<string>  stringlist;

    int  nPhysicalNames, ii, jj, j, ee, ind, index, count, tag;
    vector<string>        PhysicalNames;
    vector<int>           vecIntTemp, PhysicalNamesDims, PhysicalNamesTags, PhysicalNamesTagsDomains, PhysicalNamesTagsBoundaries;
    vector<vector<int> >  elemConnTemp;


    while(getline(infile,line))
    {
      boost::trim(line);

      // Physical names --- for BCs etc
      //
      if(line == "$PhysicalNames")
      {
        getline(infile,line);

        boost::trim(line);

        nPhysicalNames = stoi(line);

        PhysicalNames.resize(nPhysicalNames);
        PhysicalNamesDims.resize(nPhysicalNames);
        PhysicalNamesTags.resize(nPhysicalNames);

        for(ii=0; ii<nPhysicalNames; ++ii)
        {
          getline(infile,line);

          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

          boost::trim(stringlist[0]);
          PhysicalNamesDims[ii] = stoi(stringlist[0]);

          boost::trim(stringlist[1]);
          PhysicalNamesTags[ii] = stoi(stringlist[1]);

          boost::trim(stringlist[2]);
          PhysicalNames[ii] = stringlist[2].substr(1, stringlist[2].length()-2);
        }

        getline(infile,line);

        ndim = *max_element(PhysicalNamesDims.begin(), PhysicalNamesDims.end());

        for(ii=0; ii<nPhysicalNames; ++ii)
        {
            if(PhysicalNamesDims[ii] == ndim)
              PhysicalNamesTagsDomains.push_back(PhysicalNamesTags[ii]);
            else
              PhysicalNamesTagsBoundaries.push_back(PhysicalNamesTags[ii]);
        }

        PetscPrintf(MPI_COMM_WORLD, " PhysicalNamesTagsDomains ... \n");
        if(this_mpi_proc == 0) printVector(PhysicalNamesTagsDomains);
        PetscPrintf(MPI_COMM_WORLD, "\n\n\n");

        PetscPrintf(MPI_COMM_WORLD, " PhysicalNamesTagsBoundaries ... \n");
        if(this_mpi_proc == 0) printVector(PhysicalNamesTagsBoundaries);
        PetscPrintf(MPI_COMM_WORLD, "\n\n\n");

        // prepare boundary patches

        numBoundaryPatches = 0;
        for(ii=0; ii<nPhysicalNames; ++ii)
        {
            if(PhysicalNamesDims[ii] < ndim)
                numBoundaryPatches++;
        }

        for(ii=0; ii<nPhysicalNames; ++ii)
        {
            if(PhysicalNamesDims[ii] < ndim)
            {
                BoundaryPatches.push_back(make_unique<BoundaryPatch>(PhysicalNames[ii], PhysicalNamesTags[ii], ndim));
            }
        }
      }

      // nodes
      //
      if(line == "$Nodes")
      {
        getline(infile,line);        boost::trim(line);

        nNode_global = stoi(line);
        PetscPrintf(MPI_COMM_WORLD, " nNode_global = %d \n", nNode_global);

        nodeCoordsOrig.resize(nNode_global);

        for(ii=0; ii<nNode_global; ++ii)
        {
          getline(infile,line);
          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

          nodeCoordsOrig[ii][0] = stof(stringlist[1]);
          nodeCoordsOrig[ii][1] = stof(stringlist[2]);
          nodeCoordsOrig[ii][2] = stof(stringlist[3]);
        }

        getline(infile,line);
      }

      if(line == "$Elements")
      {
        getline(infile,line);
        boost::trim(line);
        boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

        count = stoi(stringlist[0]);

        elemConnTemp.resize(count);
        elemConn.resize(count);

        // elm-number elm-type number-of-tags < tag > … node-number-list
        //
        nElem_global = 0;
        for(ee=0; ee<count; ++ee)
        {
          getline(infile,line);
          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

          ind = stringlist.size();

          elemConnTemp[ee].resize(ind);

          for(j=0; j<ind; ++j)
            elemConnTemp[ee][j] = stoi(stringlist[j]);

          npElem = ind-5;
          vecIntTemp.resize(npElem);
          for(j=0; j<npElem; ++j)
          {
             vecIntTemp[j] = elemConnTemp[ee][j+5] - 1;
          }
          //printVector(vecIntTemp);

          tag = elemConnTemp[ee][3];
          if( std::count(PhysicalNamesTagsDomains.begin(), PhysicalNamesTagsDomains.end(), tag) )
          {
              elemConn[nElem_global++] = vecIntTemp;
           }
          else
          {
              for(j=0; j<numBoundaryPatches; ++j)
              {
                  if( tag == BoundaryPatches[j]->getTag() )
                  {
                      index = j;
                      break;
                  }
              }
              BoundaryPatches[index]->addElement(vecIntTemp);
          }
        }

        getline(infile,line);
      }
    }


    //
    // Process data
    //
    for(auto& bpt : BoundaryPatches)
    {
        bpt->processData();
        //bpt->printData();
    }

    //cout << " nElem_global = " << nElem_global << '\t' << elemConn.size() << endl;
    while(elemConn.size() > nElem_global)
        elemConn.pop_back();
    //cout << " nElem_global = " << nElem_global << '\t' << elemConn.size() << endl;

    //for(ii=0; ii<nElem_global; ++ii)
      //printVector(elemConn[ii]);
    //cout << endl;    cout << endl;

    PetscPrintf(MPI_COMM_WORLD, " myMesh::readInputGMSH \n");

    return 0;
}





int  myMesh::readBoundariesData(ifstream& infile, string& line)
{
    vector<string>  stringlist;
    string label;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            label = stringlist[0];
            
            int index=-1, ii;
            for(ii=0; ii<numBoundaryPatches; ++ii)
            {
                if( label == BoundaryPatches[ii]->getLabel() )
                {
                    index = ii;
                    break;
                }
            }

            if(index == -1)
            {
                throw runtime_error("Boundary label in myMesh::readBoundariesData");
            }

            BoundaryPatches[index]->readData(infile, line);
        }
    }

    return 0;
}








int myMesh::prepareInputData()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n  myMesh::prepareInputData()  .... STARTED ...\n\n");

    assert(ndim > 0 && ndim < 4);

    PetscPrintf(MPI_COMM_WORLD, "\n  numBoundaryPatches = %d \n", numBoundaryPatches);

    if(this_mpi_proc == 0)
    {
      cout << " myMesh statistics .....\n" << endl;
      cout << " nElem_global        = " << '\t' << nElem_global << endl;
      cout << " nNode_global        = " << '\t' << nNode_global  << endl;
      cout << " numBoundaryPatches  = " << '\t' << numBoundaryPatches  << endl;
    }


    node_map_get_old.resize(nNode_global, 0);
    node_map_get_new.resize(nNode_global, 0);

    elem_proc_id.resize(nElem_global, 0);
    node_proc_id.resize(nNode_global, 0);

    if(n_mpi_procs == 1)
    {
        elem_start  = 0;
        elem_end    = nElem_global-1;

        node_start  = 0;
        node_end    = nNode_global-1;

        nElem_local = nElem_global;
        nNode_local = nNode_global;
        nNode_owned = nNode_global;
        //ntotdofs_local  = ntotdofs_global;

        //row_start = 0;
        //row_end   = ntotdofs_global-1;

        for(int ii=0; ii<nNode_global; ii++)
        {
          node_map_get_old[ii] = ii;
          node_map_get_new[ii] = ii;
        }
    }
    else
    {
        PetscPrintf(MPI_COMM_WORLD, "\n\n Before partitionmyMesh ... \n\n");
        partitionMesh();
        PetscPrintf(MPI_COMM_WORLD, "\n\n After partitionmyMesh ... \n\n"); 
    }

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);

    MESH_PREPARED = true;

    PetscPrintf(MPI_COMM_WORLD, "\n\n  myMesh::prepareInputData()  .... FINISHED ...\n\n");

    return 0;
}













#include "femGrowth.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "NewElementMagnetoMech.h"
#include "NewMaterial.h"
#include <algorithm>
#include "utilitiesmaterial.h"
#include "QuadratureUtil.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>


extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime           myTime;
extern bool debug;

using namespace std;


int ElementBase::elemcount = 0;



femGrowth::femGrowth()
{

}


femGrowth::~femGrowth()
{

}



void femGrowth::printLogo()
{
  cout << "\n\n\n";
  cout << "   +----------------------------------------------------+  \n";
  cout << "   | ************************************************** |  \n";
  cout << "   | **                                              ** |  \n";
  cout << "   | **        MULTI-PHYSICS ANALYSIS PROGRAM        ** |  \n";
  cout << "   | **      ==================================      ** |  \n";
  cout << "   | **        Magneto Mechanics Problems            ** |  \n";
  cout << "   | **                                              ** |  \n";
  cout << "   | **                                              ** |  \n";
  cout << "   | **          Chennakesava Kadapa 2022            ** |  \n";
  cout << "   | **                                                 |  \n";
  cout << "   | ************************************************** |  \n";
  cout << "   +----------------------------------------------------|  \n";
  cout << "                                                           \n";

  return;
}








int femGrowth::prepareInputData()
{
    if(debug) PetscPrintf(PETSC_COMM_WORLD, "\n     femGrowth::prepareInputData()  .... STARTED ...\n");

    int  idd = ElementTypeDataList[0]->getElemTypeNameNum();
    cout << "idd = " << idd << endl;

    // ndof = ndim for continuum elements
    ndof = ndim;

    if( idd == ELEM_TRUSS_2D_NODES2 )
      ndof = 2;
    else if( idd == ELEM_TRUSS_3D_NODES2 )
      ndof = 3;
    if( (idd == ELEM_BEAM_2D_NODES2) || (idd == ELEM_BEAM_2D_NODES3) )
      ndof = 3;
    if( (idd == ELEM_BEAM_3D_NODES2) || (idd == ELEM_BEAM_3D_NODES3) )
      ndof = 6;
    else if( idd == ELEM_SHELL_FLAT_QUAD4 )
      ndof = 6;

    if(this_mpi_proc == 0)
    {
        cout << " ndim  = " << ndim         << endl;
        cout << " nElem = " << nElem_global << endl;
        cout << " nNode = " << nNode_global << endl;
        cout << " ndof  = " << ndof         << endl;
    }

    // ==================================================
    //
    // Partition the mesh if more than one processor
    //
    // ==================================================

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

        for(int ii=0; ii<nNode_global; ii++)
        {
          node_map_get_old[ii] = ii;
          node_map_get_new[ii] = ii;
        }
    }
    else
    {
        PetscPrintf(MPI_COMM_WORLD, "\n\n Before partitionMesh ... \n\n");
        PetscPrintf(MPI_COMM_WORLD, "\n\n n_mpi_procs %d ... \n\n", n_mpi_procs);

        //partitionMesh();

        PetscPrintf(MPI_COMM_WORLD, "\n\n After partitionMesh ... \n\n");
    }

    errpetsc = MPI_Barrier(MPI_COMM_WORLD);


    MIXED_ELEMENT      = false;
    MIXED_ELEMENT_P0   = false;
    MIXED_STAB_ELEMENT = false;

    if( (idd == ELEM_MAGNMECH_2D_HM_21)          ||
        (idd == ELEM_MAGNMECH_3D_HM_21)          )
    {
      MIXED_ELEMENT = true;

      for(auto& mat : MatlDataList)
        mat->MIXED_ELEMENT = true;
    }

    if( (idd == ELEM_MAGNMECH_2D_HM_10) || (idd == ELEM_MAGNMECH_3D_HM_10) )
    {
      MIXED_ELEMENT_P0 = true;

      for(auto& mat : MatlDataList)
        mat->MIXED_ELEMENT = true;
    }

    // set the flag for internal variables for elastoplastic materials
    //if( (elems[0]->nivGP) > 0)
      //intVarFlag = true;
    //else
      intVarFlag = false;


    cout << " nElem_global  = " << nElem_global  << endl;
    cout << " nNode_global  = " << nNode_global  << endl;
    cout << " ndim   = " << ndim << endl;
    cout << " ndof   = " << ndof << endl;
    cout << " intVarFlag   " << intVarFlag << endl;
    cout << " MIXED_ELEMENT      = " << MIXED_ELEMENT << endl;
    cout << " MIXED_ELEMENT_P0   = " << MIXED_ELEMENT_P0   << endl;
    cout << " MIXED_STAB_ELEMENT = " << MIXED_STAB_ELEMENT << endl;

    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    int ii, jj, kk, ee, nn, ind, n1, n2, dof, tag, index, matId=0, elemId=0, npElem;

    vector<int>  nodeNums;

    elems = new ElementBase* [nElem_global];

    for(ee=0;ee<nElem_global;ee++)
    {
      printVector(elemConn[ee]);

      npElem = elemConn[ee].size()-5;
      nodeNums.resize(npElem);
      for(ii=0; ii<npElem; ii++)
        nodeNums[ii] = elemConn[ee][5+ii]-1;

      tag = elemConn[ee][3];
      index = -1;
      for(ii=0; ii<Domains.size(); ii++)
      {
          if(tag == Domains[ii]->getTag())
          {
            index = ii;
            break;
          }
      }
      //cout << " index = " << index << endl;

      matId  = Domains[index]->matId;
      elemId = Domains[index]->elemId;

      //cout << " elemId = " << elemId << '\t' << ElementTypeDataList[elemId]->getElemTypeNameNum() << endl;

      elems[ee] = NewElementMagnetomechFEM(ElementTypeDataList[elemId]->getElemTypeNameNum());

      //cout << " matId  = " << matId << endl;

      elems[ee]->elenum       = ee;
      elems[ee]->nodeNums     = nodeNums;

      elems[ee]->elmType  =  elemId;
      elems[ee]->matType  =  matId;
      elems[ee]->secType  =  1;

      elems[ee]->SolnData     = &(SolnData);
      elems[ee]->GeomData     = &(GeomData);
      elems[ee]->MatlData     = MatlDataList[matId];
      elems[ee]->ElemTypeData = ElementTypeDataList[elemId];

      elems[ee]->prepareElemData();
    }

    PetscPrintf(MPI_COMM_WORLD, "\n elements are created and prepared \n");

    numBoundaryConditions = boundaryConitions.size();
    numTractionConditions = tractionConitions.size();

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    //SolnData.setTimeIncrementType(tis);
    //SolnData.setSpectralRadius(rhoInfty);
    SolnData.initialise(nNode_global*ndof, nNode_global*ndof, 0, 0);

    // node-to-element connectivity for post-processing

    vector<int>  vecTemp2;

    node_elem_conn.resize(nNode_global);

    for(ee=0;ee<nElem_global;ee++)
    {
      vecTemp2 = elems[ee]->nodeNums;
      for(ii=0; ii<vecTemp2.size(); ii++)
      {
        node_elem_conn[vecTemp2[ii]].push_back(ee);
      }
    }

    // find the unique values and sort the entries
    for(ee=0;ee<nNode_global;ee++)
    {
      findUnique(node_elem_conn[ee]);
    }

    // create the arrays for proce midnodes
    // index 0 --- 1 - midnode, 0 - not midnode
    // index 1 --- connecting node 1
    // index 2 --- connecting node 2
    midnodeData.resize(nNode_global);
    for(ii=0; ii<nNode_global; ii++)
    {
      midnodeData[ii].resize(3);
      // set the default value to 'not a mid-node'
      midnodeData[ii][0] = 0;
    }

    // adjust mid-nodes for quadratic Bezier triangle/tetrahedron element
    if( (idd==9) || (idd==20) || ((idd >= 201) && (idd <= 230)) )
    {
      //processForBernsteinElem_globalents();
    }

    cout << " MIXED_ELEMENT = " << MIXED_ELEMENT << endl;

    cout << " dispDegree  = " << dispDegree << endl;
    cout << " presDegree  = " << presDegree << endl;
    cout << " mpotDegree  = " << mpotDegree << endl;


    if(debug) PetscPrintf(PETSC_COMM_WORLD, "     femGrowth::prepareInputData()  .... FINISHED ...\n\n");

    return 0;
}





int femGrowth::processForBernsteinElements()
{
/*
    int  ii, jj, ee, n1, n2, nn, dof;
    double xx, yy, zz, disp[3];
    vector<int>  vecTemp;


    // create the arrays for process midnodes
    // index 0 --- 1 - midnode, 0 - not midnode
    // index 1 --- connecting node 1
    // index 2 --- connecting node 2
    midnodeData.resize(nNode);
    for(ii=0; ii<nNode; ii++)
    {
      midnodeData[ii].resize(3);
      // set the default value to 'not a mid-node'
      midnodeData[ii][0] = 0;
    }


    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;

    cout << " idd = " << idd << endl;

    if( (idd == 2010) || (idd == 2060) )
      return;


    // loop over all the elements to
    // a.) find the midnodes, and
    // b.) adjacent nodes to midnodes


    // 2D elements
    if(ndim == 2)
    {
      // triangular elements
      if( (idd == 2001) || (idd == 2003) )
      {
        int midnodemaptria6[6][2] = {{0,0}, {0,0}, {0,0},
                                     {0,1}, {1,2}, {2,0}};

        for(ee=0; ee<nElem; ee++)
        {
          vecTemp = elems[ee]->nodeNums;
          for(ii=3; ii<=5; ii++)
          {
            nn = vecTemp[ii];

            midnodeData[nn][0] = 1;
            midnodeData[nn][1] = vecTemp[midnodemaptria6[ii][0]];
            midnodeData[nn][2] = vecTemp[midnodemaptria6[ii][1]];
          }
        }
      }

      // quadrilateral elements
      if( (idd == 2002) || (idd == 2004) )
      {
        int midnodemapquad9[9][2] = {{0,0}, {0,0}, {0,0}, {0,0},
                                     {0,1}, {1,2}, {2,3}, {3,0}, {0,0}};

        for(ee=0; ee<nElem; ee++)
        {
          vecTemp = elems[ee]->nodeNums;
          for(ii=4; ii<=7; ii++)
          {
            nn = vecTemp[ii];

            midnodeData[nn][0] = 1;
            midnodeData[nn][1] = vecTemp[midnodemapquad9[ii][0]];
            midnodeData[nn][2] = vecTemp[midnodemapquad9[ii][1]];
          }
        }
      }

      // loop over the nodes and adjust nodal coordinates
      for(nn=0; nn<nNode; nn++)
      {
        if( midnodeData[nn][0] )
        {
          n1 = midnodeData[nn][1];
          n2 = midnodeData[nn][2];

          xx = 0.25*GeomData.NodePosOrig[n1][0] + 0.25*GeomData.NodePosOrig[n2][0];
          yy = 0.25*GeomData.NodePosOrig[n1][1] + 0.25*GeomData.NodePosOrig[n2][1];

          GeomData.NodePosOrig[nn][0] = 2.0*(GeomData.NodePosOrig[nn][0] - xx);
          GeomData.NodePosOrig[nn][1] = 2.0*(GeomData.NodePosOrig[nn][1] - yy);
        }
      }
    }
    else // 3D elements
    {
      // 10-noded tetrahedron element
      if( (idd == 2051) || (idd == 2054) )
      {
        int midnodemaptet10[10][2] = {{0,0}, {0,0}, {0,0}, {0,0},
                                      {0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3}};

        for(ee=0; ee<nElem; ee++)
        {
          vecTemp = elems[ee]->nodeNums;

          for(ii=4; ii<=9; ii++)
          {
            nn = vecTemp[ii];

            midnodeData[nn][0] = 1;
            midnodeData[nn][1] = vecTemp[midnodemaptet10[ii][0]];
            midnodeData[nn][2] = vecTemp[midnodemaptet10[ii][1]];
          }
        }
      }

      // 18-noded wedge element
      if( (idd == 2052) || (idd == 2055) )
      {
        int midnodemapwedge18[15][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0},
                                        {0,1}, {0,2}, {0,3}, {1,2}, {1,4}, {2,5}, {3,4}, {3,5}, {4,5}};

        // loop over all the elements to find the midnodes
        for(ee=0; ee<nElem; ee++)
        {
          vecTemp = elems[ee]->nodeNums;

          for(ii=6; ii<=14; ii++)
          {
            nn = vecTemp[ii];

            midnodeData[nn][0] = 1;
            midnodeData[nn][1] = vecTemp[midnodemapwedge18[ii][0]];
            midnodeData[nn][2] = vecTemp[midnodemapwedge18[ii][1]];
          }
        }
      }

      // 27-noded hexahedron element
      if( (idd == 2053) || (idd == 2056) )
      {
        int midnodemaphex27[20][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0},
                                      {0,1}, {0,3}, {0,4}, {1,2}, {1,5}, {2,3}, {2,6}, {3,7}, {4,5}, {4,7}, {5,6}, {6,7}};

        // loop over all the elements to find the midnodes
        for(ee=0; ee<nElem; ee++)
        {
          vecTemp = elems[ee]->nodeNums;

          for(ii=8; ii<=19; ii++)
          {
            nn = vecTemp[ii];

            midnodeData[nn][0] = 1;
            midnodeData[nn][1] = vecTemp[midnodemaphex27[ii][0]];
            midnodeData[nn][2] = vecTemp[midnodemaphex27[ii][1]];
          }
        }
      }

      // loop over the nodes and adjust nodal coordinates
      for(nn=0; nn<nNode; nn++)
      {
        //printVector(midnodeData[nn]);
        if( midnodeData[nn][0] )
        {
          n1 = midnodeData[nn][1];
          n2 = midnodeData[nn][2];

          xx = 0.25*GeomData.NodePosOrig[n1][0] + 0.25*GeomData.NodePosOrig[n2][0];
          yy = 0.25*GeomData.NodePosOrig[n1][1] + 0.25*GeomData.NodePosOrig[n2][1];
          zz = 0.25*GeomData.NodePosOrig[n1][2] + 0.25*GeomData.NodePosOrig[n2][2];

          GeomData.NodePosOrig[nn][0] = 2.0*(GeomData.NodePosOrig[nn][0] - xx);
          GeomData.NodePosOrig[nn][1] = 2.0*(GeomData.NodePosOrig[nn][1] - yy);
          GeomData.NodePosOrig[nn][2] = 2.0*(GeomData.NodePosOrig[nn][2] - zz);
        }
      }
    }
*/
    return 0;
}




int  femGrowth::prepareDataForMixedElements()
{
  int  idd = ElementTypeDataList[0]->getElemTypeNameNum();

  if( (idd == ELEM_SOLID_2D_TRIA7_MIXED1)   ||
      (idd == ELEM_SOLID_2D_TRIA10_MIXED1)  ||
      (idd == ELEM_SOLID_3D_TETRA11_MIXED1) ||
      (idd == ELEM_SOLID_3D_TETRA20_MIXED1) )
  {
      presDOF = nElem_global*elems[0]->nlbfP;

      assyForSolnPres.resize(presDOF);
      for(int ii=0; ii<presDOF; ii++)
      {
          assyForSolnPres[ii] = ii;
      }
      SolnData.pres.resize(presDOF);
  }
  else
  {
    // gather all the pressure nodes (dofs)
    //
    pressure_nodes.clear();

    int  ee, ii, jj, ind;
    vector<int>  nodeNumsPres, forAssyVecPres;

    for(ee=0; ee<nElem_global; ee++)
    {
        nodeNumsPres = elems[ee]->nodeNumsPres;

        for(ii=0; ii<nodeNumsPres.size(); ii++)
          pressure_nodes.push_back(nodeNumsPres[ii]);
    }

    findUnique(pressure_nodes);
    printVector(pressure_nodes);

    nNode_Pres = pressure_nodes.size();

    cout << " nNode_Pres = " << nNode_Pres << endl;

    vector<bool>  NodeType_pres(nNode_global,true); // all fixed
    vector<int>   ID_pres(nNode_global,-1);

    for(ii=0; ii<nNode_Pres; ii++)
    {
      NodeType_pres[pressure_nodes[ii]] = false;
    }

    //cout << "NodeType_pres" << endl;
    //for(ii=0;ii<nNode_global;ii++)
      //cout << ii << '\t' << NodeType_pres[ii] << endl;

    presDOF = 0;
    for(ii=0;ii<nNode_global;ii++)
    {
      if(!NodeType_pres[ii])
        ID_pres[ii] = presDOF++;
    }
    //cout << " presDOF = " << presDOF << endl;
    //cout << "ID_pres" << endl;
    //printVector(ID_pres);

    assyForSolnPres.resize(presDOF);
    int count = 0;
    for(ii=0;ii<nNode_global;ii++)
    {
      if( ID_pres[ii] != -1)
        assyForSolnPres[count++] = ii;
    }

    pressure_nodes_map_g2l.resize(nNode_global,-1);
    for(ii=0; ii<nNode_Pres; ii++)
    {
      pressure_nodes_map_g2l[pressure_nodes[ii]] = ii;
    }
    //cout << "pressure_nodes" << endl;
    //printVector(pressure_nodes);
    //cout << "pressure_nodes_map_g2l" << endl;
    //printVector(pressure_nodes_map_g2l);
    //cout << "assyForSolnPres" << endl;
    //printVector(assyForSolnPres);

    // specified pressure DOF
    //for(ii=0;ii<DirichletBCs_Pres.size();ii++)
    //{
      //SolnData.presApplied(pressure_nodes_map_g2l[DirichletBCs_Pres[ii][0]]) = DirichletBCs_Pres[ii][2];
    //}

    //printVector(SolnData.presApplied);

    // reassign the pressure nodes
    for(ee=0; ee<nElem_global; ee++)
    {
      //printVector(elems[ee]->nodeNumsPres);
      ind = elems[ee]->nodeNumsPres.size();

      for(ii=0; ii<ind; ii++)
      {
        elems[ee]->forAssyVecPres[ii] = ID_pres[elems[ee]->nodeNumsPres[ii]];
      }
      //printVector(elems[ee]->forAssyVecPres);
    }

    SolnData.pres.resize(max(nElem_global, nNode_global));
  }

  SolnData.pres.setZero();

  SolnData.presPrev = SolnData.pres;
  SolnData.presApplied = SolnData.pres;

  return 0;
}


int  femGrowth::prepareDataForPressure()
{
/*
    // this is the subroutine for the assumption that pressure is continuous
    // across material interfaces

    if( presDegree <= 0 )
    {
      presDOF = 0;

      return;
    }

    SolnData.var2.resize(max(nElem, nNode));
    SolnData.var2.setZero();

    SolnData.var2applied = SolnData.var2;
    SolnData.var2Prev    = SolnData.var2;
    SolnData.var2Prev2   = SolnData.var2;


    int  ee, ii, jj, ind;
    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;
    vector<bool>  NodeType_pres(nNode,true); // all fixed
    vector<int>   ID_pres(nNode,-1), vecTemp;

    // gather all the pressure nodes (dofs)
    //
    pressure_nodes.clear();

    if( (idd == 2001) || (idd == 2003) )                    // Tria6/3
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecPres.clear();

        for(ii=0; ii<3; ii++)
        {
            jj = elems[ee]->nodeNums[ii];

            pressure_nodes.push_back(jj);

            elems[ee]->forAssyVecPres.push_back(jj);
        }
      }
    }
    else if(
         (idd == 2002) || (idd == 2004)                     // Quad9/4
      || (idd == 2051) || (idd == 2054)                     // Tet10/4
      )
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecPres.clear();

        for(ii=0; ii<4; ii++)
        {
            jj = elems[ee]->nodeNums[ii];

            pressure_nodes.push_back(jj);

            elems[ee]->forAssyVecPres.push_back(jj);
        }
      }
    }
    else if( (idd == 2052) || (idd == 2055) )               // Wedge18/6
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecPres.clear();

        for(ii=0; ii<6; ii++)
        {
            jj = elems[ee]->nodeNums[ii];

            pressure_nodes.push_back(jj);

            elems[ee]->forAssyVecPres.push_back(jj);
        }
      }
    }
    else if( (idd == 2053) || (idd == 2056) )               // Hexa27/8
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecPres.clear();

        for(ii=0; ii<8; ii++)
        {
            jj = elems[ee]->nodeNums[ii];

            pressure_nodes.push_back(jj);

            elems[ee]->forAssyVecPres.push_back(jj);
        }
      }
    }

    //printVector(pressure_nodes);

    findUnique(pressure_nodes);
    nNode_Pres = pressure_nodes.size();
    //printVector(pressure_nodes);

    cout << " nNode_Pres = " << nNode_Pres << endl;

    // fix mid nodes
    for(ii=0;ii<nNode_Pres;ii++)
    {
      NodeType_pres[pressure_nodes[ii]] = false;
    }

    // specified pressure DOF
    for(ii=0;ii<DirichletBCs_Pres.size();ii++)
    {
      NodeType_pres[DirichletBCs_Pres[ii][0]] = true;
    }

    //cout << "NodeType_pres" << endl;
    //for(ii=0;ii<nNode;ii++)
      //cout << ii << '\t' << midnodeData[ii][0] << '\t' <<  NodeType_pres[ii] << endl;

    presDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if(!NodeType_pres[ii])
        ID_pres[ii] = presDOF++;
    }
    cout << " presDOF = " << presDOF << endl;
    //cout << "ID_pres" << endl;
    //printVector(ID_pres);

    assyForSolnPres.resize(presDOF);
    int count = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if( ID_pres[ii] != -1)
        assyForSolnPres[count++] = ii;
    }

    //cout << "pressure_nodes" << endl;
    //printVector(pressure_nodes);
    //cout << "assyForSolnPres" << endl;
    //printVector(assyForSolnPres);

    pressure_nodes_map_g2l.resize(nNode,-1);
    for(ii=0; ii<nNode_Pres; ii++)
    {
      pressure_nodes_map_g2l[pressure_nodes[ii]] = ii;
    }
    //cout << "pressure_nodes_map_g2l" << endl;
    //printVector(pressure_nodes_map_g2l);

    // specified pressure DOF
    for(ii=0;ii<DirichletBCs_Pres.size();ii++)
    {
      SolnData.var2applied(pressure_nodes_map_g2l[DirichletBCs_Pres[ii][0]]) = DirichletBCs_Pres[ii][2];
    }

    // reassign the pressure nodes
    for(ee=0; ee<nElem; ee++)
    {
      //printVector(elems[ee]->forAssyVecPres);
      jj = elems[ee]->forAssyVecPres.size();

      for(ii=0; ii<jj; ii++)
      {
        elems[ee]->forAssyVecPres[ii] = ID_pres[elems[ee]->forAssyVecPres[ii]];
      }
      //printVector(elems[ee]->forAssyVecPres);
    }
*/

    return 0;
}


/*
void  femGrowth::prepareDataForPressure()
{
    // this is the subroutine for the assumption that pressure is discontinuous
    // across material interfaces

    if( presDegree <= 0 )
    {
      presDOF = 0;

      return;
    }


    SolnData.var2.resize(max(nElem, nNode));
    SolnData.var2.setZero();

    SolnData.var2applied = SolnData.var2;
    SolnData.var2Prev    = SolnData.var2;
    SolnData.var2Prev2   = SolnData.var2;

    // gather all the pressure nodes (dofs)
    //

    vector<vector<int> >  pressure_nodes(numMaterials);
    //pressure_nodes.resize(numMaterials);

    int  ee, mm, ii, jj, matNum;
    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;
    vector<int>  vecIntTemp, nNode_Pres_Vec(numMaterials, 0);

    cout << " prepareDataForPressure ... idd = " << idd << endl;

    if( (idd == 2001) || (idd == 2003) )                    // Bezier Tria6/3
    {
      for(ee=0; ee<nElem; ++ee)
      {
        matNum = elems[ee]->matType;
        vecIntTemp = elems[ee]->nodeNums;

        pressure_nodes[matNum].push_back(vecIntTemp[0]);
        pressure_nodes[matNum].push_back(vecIntTemp[1]);
        pressure_nodes[matNum].push_back(vecIntTemp[2]);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[0]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[1]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[2]);
      }
    }
    else if( (idd == 2002) || (idd == 2004) ||              // Bezier Quad9/4
             (idd == 2051) || (idd == 2054)                 // Bezier Tet10/4
           )
    {
      for(ee=0; ee<nElem; ++ee)
      {
        matNum = elems[ee]->matType;
        vecIntTemp = elems[ee]->nodeNums;

        pressure_nodes[matNum].push_back(vecIntTemp[0]);
        pressure_nodes[matNum].push_back(vecIntTemp[1]);
        pressure_nodes[matNum].push_back(vecIntTemp[2]);
        pressure_nodes[matNum].push_back(vecIntTemp[3]);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[0]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[1]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[2]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[3]);
      }
    }
    else if( (idd == 2052) || (idd == 2055) )               // Wedge18/6
    {
      for(ee=0; ee<nElem; ++ee)
      {
        matNum = elems[ee]->matType;
        vecIntTemp = elems[ee]->nodeNums;

        pressure_nodes[matNum].push_back(vecIntTemp[0]);
        pressure_nodes[matNum].push_back(vecIntTemp[1]);
        pressure_nodes[matNum].push_back(vecIntTemp[2]);
        pressure_nodes[matNum].push_back(vecIntTemp[3]);
        pressure_nodes[matNum].push_back(vecIntTemp[4]);
        pressure_nodes[matNum].push_back(vecIntTemp[5]);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[0]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[1]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[2]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[3]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[4]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[5]);
      }
    }
    else if( (idd == 2053) || (idd == 2056) )               // Hexa27/8
    {
      for(ee=0; ee<nElem; ++ee)
      {
        matNum = elems[ee]->matType;
        vecIntTemp = elems[ee]->nodeNums;

        pressure_nodes[matNum].push_back(vecIntTemp[0]);
        pressure_nodes[matNum].push_back(vecIntTemp[1]);
        pressure_nodes[matNum].push_back(vecIntTemp[2]);
        pressure_nodes[matNum].push_back(vecIntTemp[3]);
        pressure_nodes[matNum].push_back(vecIntTemp[4]);
        pressure_nodes[matNum].push_back(vecIntTemp[5]);
        pressure_nodes[matNum].push_back(vecIntTemp[6]);
        pressure_nodes[matNum].push_back(vecIntTemp[7]);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[0]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[1]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[2]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[3]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[4]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[5]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[6]);
        elems[ee]->forAssyVecPres.push_back(vecIntTemp[7]);
      }
    }
    else
    {
      cerr << " Wrong element type in 'femGrowth::prepareDataForPressure()' " << endl;
      //exit(999);
    }

    nNode_Pres = 0;
    for(mm=0; mm<numMaterials; mm++)
    {
      //printVector(pressure_nodes);
      findUnique(pressure_nodes[mm]);

      nNode_Pres_Vec[mm] = pressure_nodes[mm].size();

      nNode_Pres += nNode_Pres_Vec[mm];
      //printVector(pressure_nodes);
    }

    cout << " nNode_Pres = " << nNode_Pres << endl;

    // all fixed
    vector<vector<bool> >  NodeType_pres(numMaterials,  vector<bool>(nNode,true));
    vector<vector<int> >   ID_pres(numMaterials, vector<int>(nNode,-1));

    // fix mid nodes

    for(mm=0; mm<numMaterials; mm++)
    {
      for(ii=0;ii<nNode_Pres_Vec[mm];ii++)
      {
        NodeType_pres[mm][pressure_nodes[mm][ii]] = false;
      }
    }

    // specified pressure DOF
    for(mm=0; mm<numMaterials; mm++)
    {
      for(ii=0;ii<DirichletBCs_Pres.size();ii++)
      {
        NodeType_pres[mm][DirichletBCs_Pres[ii][0]] = true;
      }
    }

    //cout << "NodeType_pres" << endl;
    //for(ii=0;ii<nNode;ii++)
      //cout << ii << '\t' << NodeType_pres[ii] << endl;

    presDOF = 0;
    for(mm=0; mm<numMaterials; mm++)
    {
      for(ii=0;ii<nNode;ii++)
      {
        if(!NodeType_pres[mm][ii])
          ID_pres[mm][ii] = presDOF++;
      }
    }
    cout << " presDOF = " << presDOF << endl;
    //cout << "ID_pres" << endl;
    //printVector(ID_pres);

    assyForSolnPres.resize(presDOF);
    int count = 0;
    for(mm=0; mm<numMaterials; mm++)
    {
      for(ii=0;ii<nNode;ii++)
      {
        if( ID_pres[mm][ii] != -1)
          assyForSolnPres[count++] = ii;
      }
    }

    cout << "pressure_nodes" << endl;

    pressure_nodes_map_g2l.resize(nNode,-1);
    for(mm=0; mm<numMaterials; mm++)
    {
      for(ii=0; ii<nNode_Pres_Vec[mm]; ii++)
      {
        pressure_nodes_map_g2l[pressure_nodes[mm][ii]] = ii;
      }
    }
    cout << "pressure_nodes" << endl;
    //printVector(pressure_nodes);
    //cout << "pressure_nodes_map_g2l" << endl;
    //printVector(pressure_nodes_map_g2l);
    //cout << "assyForSolnPres" << endl;
    //printVector(assyForSolnPres);

    // specified pressure DOF
    for(ii=0;ii<DirichletBCs_Pres.size();ii++)
    {
      SolnData.var2applied(pressure_nodes_map_g2l[DirichletBCs_Pres[ii][0]]) = DirichletBCs_Pres[ii][2];
    }

    //printVector(SolnData.var2applied);

    // reassign the pressure nodes
    for(ee=0; ee<nElem; ee++)
    {
      matNum = elems[ee]->matType;

      //printVector(elems[ee]->forAssyVecPres);

      jj = elems[ee]->forAssyVecPres.size();

      for(ii=0; ii<jj; ii++)
      {
        elems[ee]->forAssyVecPres[ii] = ID_pres[matNum][elems[ee]->forAssyVecPres[ii]];
      }
      //printVector(elems[ee]->forAssyVecPres);
    }

    return;
}
*/




int  femGrowth::prepareDataForMagneticPotential()
{
/*
    if( mpotDegree == -1)
        return;

    SolnData.var3.resize(nNode);
    SolnData.var3.setZero();
    SolnData.var3Cur     = SolnData.var3;
    SolnData.var3Prev    = SolnData.var3;
    SolnData.var3Prev2   = SolnData.var3;
    SolnData.var3applied = SolnData.var3;


    int  ee, ii, jj;
    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;

    vector<bool>  NodeType_mpot(nNode,false); // all free
    vector<int>   ID_mpot(nNode,-1);

    // gather all the magnetic potential field nodes
    //
    magnpote_nodes.clear();

    if( mpotDegree == dispDegree )
    {
      nNode_Mpot = nNode;

      magnpote_nodes.resize(nNode);
      for(ii=0; ii<nNode; ii++)
        magnpote_nodes[ii] = ii;

      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecMpot = elems[ee]->nodeNums;
      }

      // specified/fixed electric potential DOF
      for(ii=0;ii<DirichletBCs_Mpot.size();ii++)
      {
        NodeType_mpot[DirichletBCs_Mpot[ii][0]] = true;
      }
    }
    else
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecMpot.clear();

        for(ii=0; ii<=ndim; ii++)
        {
          jj = elems[ee]->nodeNums[ii];

          magnpote_nodes.push_back(jj);

          elems[ee]->forAssyVecMpot.push_back(jj);
        }
      }

      //printVector(magnpote_nodes);

      findUnique(magnpote_nodes);
      nNode_Mpot = magnpote_nodes.size();
      //printVector(magnpote_nodes);

      // fix mid nodes
      for(ii=0;ii<nNode;ii++)
      {
        if(midnodeData[ii][0])
        {
          NodeType_mpot[ii] = true;
        }
      }

      // specified/fixed electric potential DOF
      for(ii=0;ii<DirichletBCs_Mpot.size();ii++)
      {
        if(!midnodeData[DirichletBCs_Mpot[ii][0]][0])
        {
          NodeType_mpot[DirichletBCs_Mpot[ii][0]] = true;
        }
      }
    }

    cout << " nNode_Mpot = " << nNode_Mpot << endl;
    magnpote_nodes_map_g2l.resize(nNode,-1);
    for(ii=0; ii<nNode_Mpot; ii++)
    {
      magnpote_nodes_map_g2l[magnpote_nodes[ii]] = ii;
    }
    //cout << "magnpote_nodes" << endl;
    //printVector(magnpote_nodes);
    //cout << "magnpote_nodes_map_g2l" << endl;
    //printVector(magnpote_nodes_map_g2l);

    //cout << "NodeType_mpot" << endl;
    //for(ii=0;ii<nNode;ii++)
      //cout << ii << '\t' << NodeType_mpot[ii] << endl;

    mpotDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if(!NodeType_mpot[ii])
        ID_mpot[ii] = mpotDOF++;
    }
    cout << " mpotDOF = " << mpotDOF << endl;
    //cout << "ID_mpot" << endl;
    //printVector(ID_mpot);

    assyForSolnMpot.resize(mpotDOF);
    int count = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if( ID_mpot[ii] != -1)
        assyForSolnMpot[count++] = ii;
    }

    //cout << "assyForSolnMpot" << endl;
    //printVector(assyForSolnMpot);

    // reassign the pressure nodes
    for(ee=0; ee<nElem; ee++)
    {
      //printVector(elems[ee]->forAssyVecMpot);
      jj = elems[ee]->forAssyVecMpot.size();

      for(ii=0; ii<jj; ii++)
      {
        elems[ee]->forAssyVecMpot[ii] = ID_mpot[elems[ee]->forAssyVecMpot[ii]];
      }
      //printVector(elems[ee]->forAssyVecMpot);
    }
*/
    return 0;
}






int  femGrowth::prepareDataForTemperature()
{
/*
    if(tempDegree ==  -1)
    {
        tempDOF = 0;
        return;
    }

    SolnData.var4.resize(nNode);
    SolnData.var4.setZero();
    SolnData.var4Cur      = SolnData.var4;
    SolnData.var4Dot      = SolnData.var4;
    SolnData.var4DotPrev  = SolnData.var4;
    SolnData.var4DotCur   = SolnData.var4;
    SolnData.var4applied  = SolnData.var4;


    int  ee, ii, jj;
    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;

    vector<bool>  NodeType_temp(nNode,false); // all free
    vector<int>   ID_temp(nNode,-1);

    // gather all the electric potential field nodes
    //
    temperature_nodes.clear();

    if( tempDegree == dispDegree )
    {
      nNode_Temp = nNode;

      temperature_nodes.resize(nNode);
      for(ii=0; ii<nNode; ii++)
        temperature_nodes[ii] = ii;

      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecTemp = elems[ee]->nodeNums;
      }

      // specified/fixed temperature DOF
      for(ii=0;ii<DirichletBCs_Temp.size();ii++)
      {
        NodeType_temp[DirichletBCs_Temp[ii][0]] = true;
      }
    }
    else
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecTemp.clear();

        for(ii=0; ii<=ndim; ii++)
        {
          jj = elems[ee]->nodeNums[ii];

          temperature_nodes.push_back(jj);

          elems[ee]->forAssyVecTemp.push_back(jj);
        }
      }

      //printVector(temperature_nodes);

      findUnique(temperature_nodes);
      nNode_Temp = temperature_nodes.size();
      //printVector(temperature_nodes);

      // fix mid nodes
      for(ii=0;ii<nNode;ii++)
      {
        if(midnodeData[ii][0])
        {
          NodeType_temp[ii] = true;
        }
      }

      // specified/fixed electric potential DOF
      for(ii=0;ii<DirichletBCs_Temp.size();ii++)
      {
        if(!midnodeData[DirichletBCs_Temp[ii][0]][0])
        {
          NodeType_temp[DirichletBCs_Temp[ii][0]] = true;
        }
      }
    }

    cout << " nNode_Temp = " << nNode_Temp << endl;
    temperature_nodes_map_g2l.resize(nNode,-1);
    for(ii=0; ii<nNode_Mpot; ii++)
    {
      temperature_nodes_map_g2l[magnpote_nodes[ii]] = ii;
    }
    //cout << "magnpote_nodes" << endl;
    //printVector(magnpote_nodes);
    //cout << "temperature_nodes_map_g2l" << endl;
    //printVector(temperature_nodes_map_g2l);

    //cout << "NodeType_temp" << endl;
    //for(ii=0;ii<nNode;ii++)
      //cout << ii << '\t' << NodeType_temp[ii] << endl;

    tempDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if(!NodeType_temp[ii])
        ID_temp[ii] = tempDOF++;
    }
    cout << " tempDOF = " << tempDOF << endl;
    //cout << "ID_temp" << endl;
    //printVector(ID_temp);

    assyForSolnTemp.resize(tempDOF);
    int count = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if( ID_temp[ii] != -1)
        assyForSolnTemp[count++] = ii;
    }

    //cout << "assyForSolnTemp" << endl;
    //printVector(assyForSolnTemp);

    // reassign the pressure nodes
    for(ee=0; ee<nElem; ee++)
    {
      //printVector(elems[ee]->forAssyVecTemp);
      jj = elems[ee]->forAssyVecTemp.size();

      for(ii=0; ii<jj; ii++)
      {
        elems[ee]->forAssyVecTemp[ii] = ID_temp[elems[ee]->forAssyVecTemp[ii]];
      }
      //printVector(elems[ee]->forAssyVecTemp);
    }
*/
    return 0;
}



int femGrowth::setInitialConditions()
{
    PetscPrintf(MPI_COMM_WORLD, "femGrowth::setInitialConditions() ... \n");
    if(debug) cout <<  " femGrowth::setInitialConditions() ... STARTED " <<  endl;

    /*
    double  x=0.0, y=0.0, z=0.0, eps = 0.02;
    string  expr;

    SolnData.dispPrev.setZero();

    for(int dof=0; dof<ndof; ++dof)
    {
      expr = initialConitions[dof];

      cout << expr << endl;

      if( !expr.empty() )
      {
        myMathFunction  mathfun;
        mathfun.initialise(initialConitions[dof]);

        exprtk::expression<double>   expression;

        exprtk::symbol_table<double> symbol_table;
        exprtk::parser<double>       parser;

        symbol_table.add_variable("x",x);
        symbol_table.add_variable("y",y);
        symbol_table.add_variable("z",z);
        symbol_table.add_variable("eps",eps);
        symbol_table.add_constants();
        expression.register_symbol_table(symbol_table);

        parser.compile(expr, expression);

        for(int ii=0; ii<nNode_global; ++ii)
        {
          x = nodeCoordsOrig[ii][0];
          y = nodeCoordsOrig[ii][1];
          z = nodeCoordsOrig[ii][2];

          SolnData.dispPrev(node_map_get_new[ii]*ndof + dof) = expression.value();
        }
      }
    }
    PetscPrintf(MPI_COMM_WORLD, "femGrowth::setInitialConditions() ... \n");

    // add boundary Conditions
    //setBoundaryConditions();
    //solnPrev += solnApplied;

    SolnData.disp = SolnData.dispPrev;
*/

    int jj, nn, dof, dd, ind;

    double  xx=0.0, yy=0.0, zz=0.0, fact;
    double  specVal, dt=myTime.dt;

    SolnData.velo.setZero();

    for(nn=0; nn<nNode_global; nn++)
    {
        xx = GeomData.NodePosOrig[nn][0];
        yy = GeomData.NodePosOrig[nn][1];
        zz = GeomData.NodePosOrig[nn][2];
        ind = nn*ndof;

      //SolnData.velo[ind]    =  5.0*yy/3.0;
      //SolnData.velo[ind+1]  =  0.0;
      //SolnData.velo[ind+2]  =  0.0;

      //SolnData.velo[ind]    =    0.0;
      //SolnData.velo[ind+1]  =    0.0;
      //SolnData.velo[ind+2]  = xx;
    }

    //printVector(SolnData.velo);

    /*
    // adjust the velocity values for the midnoes
    VectorXd  velTemp(SolnData.velo.rows());
    velTemp = SolnData.velo;

    for(nn=0; nn<nNode_global; nn++)
    {
      if( midnodeData[nn][0] ) // if midnode
      {
        for(dd=0; dd<ndof; dd++)
        {
          xx = 0.25*velTemp(midnodeData[nn][1]*ndof+dd) + 0.25*velTemp(midnodeData[nn][2]*ndof+dd);

          SolnData.velo[nn*ndof+dd] = 2.0*(velTemp(nn*ndof+dd) - xx);
        }
        //cout << velTemp(nn*ndof) << '\t' << SolnData.velo(nn*ndof) << endl;
      }
    }
    */

    //printVector(SolnData.var1);

    if(debug) cout <<  " femGrowth::setInitialConditions() ... ENDED " <<  endl;

    return 0;
}





int femGrowth::setBoundaryConditions()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femGrowth::setBoundaryConditions() ... STARTED \n");}

    int  ii, nn, dof, index=-1, timeFuncNum, ind;
    double xc, yc, zc, value, loadFactor;
    vector<int> nodeNums;

    SolnData.dispApplied.setZero();

    // set boundary conditions
    for(auto&  bc : boundaryConitions)
    {
        cout << " Patch                 = " << bc->label << endl;
        cout << " BCType                = " << bc->BCType << endl;
        cout << " dof_specified_string  = " << bc->dof_specified_string << endl;
        cout << " dof_specified_int     = " << bc->dof_specified_int << endl;
        //cout << " expression            = " << bc->expression << endl;
        //cout << " timeFunctionNum       = " << bc->timeFunctionNum << endl;
        //cout << endl;  cout << endl;

        // nodeNums contains new node numbers which are used for the solution
        // nodeCoordsOrig contains coordinates of nodes with numbers before domain decompositions
        PetscPrintf(MPI_COMM_WORLD, " timeFunctionNum = 1 ... load Factor = %12.6f \n", timeFunctions[0]->getValue());

        nodeNums = bc->bpt->nodeNums;

        if(bc->BCType == "specified")
        {
            timeFuncNum = bc->getTimeFunctionNumber();

            if(timeFuncNum == -1)
              loadFactor = 1.0;
            else
              loadFactor = timeFunctions[timeFuncNum]->getValue();

            PetscPrintf(MPI_COMM_WORLD, " bc->label = %10s ...  bc->timeFunctionNum = %d ... load Factor = %12.6f \n", bc->label.c_str(), bc->timeFunctionNum,  loadFactor);

            myMathFunction  mathfun;
            mathfun.initialise(bc->expression);

            dof = bc->dof_specified_int;

            for(ii=0; ii<nodeNums.size(); ++ii)
            {
                //nn = node_map_get_old[nodeNums[i]];
                nn = nodeNums[ii];

                xc = GeomData.NodePosOrig[nn][0];
                yc = GeomData.NodePosOrig[nn][1];
                zc = GeomData.NodePosOrig[nn][2];

                value = mathfun.getValue(xc, yc, zc) * loadFactor;

                //cout << xc << '\t' << yc << '\t' << zc << '\t' << loadFactor << '\t' << value << endl;

                SolnData.dispApplied[nn*ndof+dof] = value;
            }
        }
        else if(bc->BCType == "fixed")
        {
            for(ii=0; ii<nodeNums.size(); ++ii)
            {
              nn = nodeNums[ii];
              ind = nn*ndof;

              for(dof=0; dof<ndim; ++dof)
              {
                SolnData.dispApplied[ind+dof] = 0.0;
              }
            }
        }
        else
        {
            throw runtime_error("Patch type not available in femGrowth::setBoundaryConditions");
        }
    }

    SolnData.dispApplied -= SolnData.disp;
    //printVector(SolnData.dispApplied);

    /*
    for(ii=0; ii<DirichletBCs.size(); ii++)
    {
      nn   = (int) (DirichletBCs[ii][0]);
      dof  = (int) (DirichletBCs[ii][1]);

      ind = nn*ndof+dof;

      if( midnodeData[nn][0] )
      {
        xx = 0.25*solnTemp(midnodeData[nn][1]*ndof+dof) + 0.25*solnTemp(midnodeData[nn][2]*ndof+dof);

        SolnData.dispApplied[ind] = 2.0*(solnTemp(ind) - xx);
      }

      SolnData.dispApplied[ind] -= SolnData.disp[ind];
    }
    */

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femGrowth::setBoundaryConditions() ... ENDED \n");}

    return 0;
}






int femGrowth::addExternalForces()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femGrowth::addExternalForces() ... STARTED \n");}

    solverPetsc->Fext.setZero();

    if(this_mpi_proc > 0) return 0;


    // set boundary conditions
    if(numTractionConditions > 0)
    {
    for(auto&  bc : tractionConitions)
    {
    if(bc->BCType == "traction")
    {
        cout << " Patch                 = " << bc->label << endl;
        cout << " BCType                = " << bc->BCType << endl;
        cout << " timeFunctionNum       = " << bc->timeFunctionNum << endl;
        cout << endl;  cout << endl;

        int  ii, dim, nn, ee, npElem, nGP, gp, TI, TIp1, TIp2, ind1, ind2;
        double xc, yc, zc, value, loadFactor;
        vector<int> nodeNums, forAssyVec;
        vector<double> gpoints1, gpoints2, gweights;

        // nodeNums contains new node numbers which are used for the solution
        // nodeCoordsOrig contains coordinates of nodes with numbers before domain decompositions

        nodeNums = bc->bpt->nodeNums;

        if(bc->timeFunctionNum == -1)
          loadFactor = 1.0;
        else
          loadFactor = timeFunctions[bc->timeFunctionNum]->getValue();

        PetscPrintf(MPI_COMM_WORLD, " bc->label = %10s ...  bc->timeFunctionNum = %d ... load Factor = %12.6f \n", bc->label.c_str(), bc->timeFunctionNum,  loadFactor);


        double  xNode[27], yNode[27], zNode[27], param[3], N[27], normal[3], tangent1[3], tangent2[3], traction[3], Jac, dvol;
        vector<double>  specValues = bc->specValues;

        specValues[0] *= loadFactor;
        specValues[1] *= loadFactor;
        specValues[2] *= loadFactor;

        VectorXd  Flocal(100);

        if(ndim == 2)
        {
            npElem = bc->bpt->elemConn[0].size();
            nGP = npElem;
            getGaussPoints1D(nGP, gpoints1, gweights);

            for(ee=0; ee<bc->bpt->elemConn.size(); ++ee)
            {
                //cout << " ee= " << ee << endl;
                nodeNums = bc->bpt->elemConn[ee];
                //printVector(nodeNums);

                npElem = nodeNums.size();

                for(ii=0; ii<npElem; ii++)
                {
                    //nn = node_map_get_old[nodeNums[i]];
                    nn = nodeNums[ii];

                    xNode[ii] = GeomData.NodePosOrig[nn][0];
                    yNode[ii] = GeomData.NodePosOrig[nn][1];
                }

                Flocal.setZero();

                for(gp=0; gp<nGP; gp++)
                {
                    param[0] = gpoints1[gp];

                    LagrangeBasisFunsEdge2D(npElem, param, xNode, yNode, N, normal, Jac);
                    dvol = gweights[gp]*Jac;

                    tangent1[0] = -normal[1];
                    tangent1[1] =  normal[0];

                    traction[0] = ( specValues[0]*normal[0] + specValues[1]*tangent1[0] )*dvol;
                    traction[1] = ( specValues[0]*normal[1] + specValues[1]*tangent1[1] )*dvol;

                    for(ii=0; ii<npElem; ii++)
                    {
                        TI = ii*ndof;

                        Flocal(TI)   += N[ii]*traction[0];
                        Flocal(TI+1) += N[ii]*traction[1];
                    }
                }//for(gp=0; gp<nGP; gp++)

                //printVector(Flocal);

                for(ii=0; ii<npElem; ii++)
                {
                  ind1 = ii*ndof;
                  ind2 = nodeNums[ii]*ndof;

                  for(dim=0; dim<ndim; dim++)
                    solverPetsc->Fext[ind2+dim]   += Flocal[ind1+dim];
                }

            }//for(ee=0;
        } //if(ndim
        else
        {
            for(ee=0; ee<bc->bpt->elemConn.size(); ++ee)
            {
                nodeNums = bc->bpt->elemConn[ee];

                npElem = nodeNums.size();

                for(ii=0; ii<npElem; ii++)
                {
                    //nn = node_map_get_old[nodeNums[i]];
                    nn = nodeNums[ii];

                    xNode[ii] = GeomData.NodePosOrig[nn][0];
                    yNode[ii] = GeomData.NodePosOrig[nn][1];
                    zNode[ii] = GeomData.NodePosOrig[nn][2];
                }

                nGP = getGaussPoints2D(npElem, gpoints1, gpoints2, gweights);

                Flocal.setZero();
                for(gp=0; gp<nGP; gp++)
                {
                    param[0] = gpoints1[gp];
                    param[1] = gpoints2[gp];

                    LagrangeBasisFunsFace(npElem, param, xNode, yNode, zNode, N, normal, tangent1, tangent2, Jac);

                    dvol = gweights[gp]*Jac;

                    traction[0] = ( specValues[0]*normal[0] + specValues[1]*tangent1[0] )*dvol;
                    traction[1] = ( specValues[0]*normal[1] + specValues[1]*tangent1[1] )*dvol;
                    traction[2] = ( specValues[0]*normal[2] + specValues[1]*tangent1[2] )*dvol;

                    //cout << "normals" << endl;
                    //cout <<   normal[0] << '\t' <<   normal[1] << '\t' <<   normal[2] << endl;
                    //cout << tangent1[0] << '\t' << tangent1[1] << '\t' << tangent1[2] << endl;
                    //cout << tangent2[0] << '\t' << tangent2[1] << '\t' << tangent2[2] << endl; cout << endl;

                    for(ii=0; ii<npElem; ii++)
                    {
                        TI = ii*ndim;

                        Flocal(TI)   += N[ii]*traction[0];
                        Flocal(TI+1) += N[ii]*traction[1];
                        Flocal(TI+2) += N[ii]*traction[2];
                    }
                }//for(gp=0; gp<nGP; gp++)


                for(ii=0; ii<npElem; ii++)
                {
                  ind1 = ii*ndim;
                  ind2 = nodeNums[ii]*ndim;

                  for(dim=0; dim<ndim; dim++)
                    solverPetsc->Fext[ind2+dim]   += Flocal[ind1+dim];
                }
            }//for(ee=0;
        }//else
    }
    }
    }

    if(nodalForces.count > 0)
    {
        int  nn, dof, ii;

        if(nodalForces.timeFunctionNum == -1)
          loadFactor = 1.0;
        else
          loadFactor = timeFunctions[nodalForces.timeFunctionNum]->getValue();

        PetscPrintf(MPI_COMM_WORLD, " Nodal Forces ...  nodalForces->timeFunctionNum = %d ... load Factor = %12.6f \n", nodalForces.timeFunctionNum,  loadFactor);

        // specified nodal forces
        for(ii=0;ii<nodalForces.count;++ii)
        {
              nn  = (int) (nodalForces.data[ii][0]);
              dof = (int) (nodalForces.data[ii][1]);

          solverPetsc->Fext[nn*ndof+dof] += loadFactor*nodalForces.data[ii][2];
        }
    }
    //printVector(solverPetsc->Fext);

    if(debug) {cout << endl; cout << endl; cout <<  " femGrowth::addExternalForces() ... ENDED " <<  endl;}

    return 0;
}






int femGrowth::timeUpdate()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, " femSolids::timeUpdate() ... STARTED \n");}

    myTime.update();

    // set parameters for the time integration scheme.
    // need to be done every time step to account for adaptive time stepping
    SolnData.setTimeParam();

    int  idd = 0;//SolnData.ElemProp[0]->id;
    // TODO: set MIXED_ELEMENT_P0
    if( MIXED_ELEMENT_P0 )
    {
        for(int ee=0; ee<nElem_global; ee++)
        {
            elems[ee]->presDOF     = elems[ee]->presDOFprev;
            elems[ee]->presDOFprev = elems[ee]->presDOFprev2;
        }
    }

    // copy internal variables intVar1 to intVar2
    //if(intVarFlag)
      //copyElemInternalVariables();

    // update time functions
    for(auto& tmf : timeFunctions)
      tmf->update();

    //SolnData.saveSolution();

    // set Dirichlet boundary conditions
    setBoundaryConditions();


    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n femSolids::timeUpdate() ... FINISHED \n\n"); }

    return 0;
}



int femGrowth::updateIterStep()
{
    if(debug)  cout << " femGrowth::updateIterStep() ... " << endl;

    SolnData.updateIterStep();

    cout << " aaaaaaaa " << endl;

    int  nn, dd, ind;
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

    if(debug)  cout << " femGrowth::updateIterStep() ... " << endl;

    return 0;
}


int femGrowth::reset()
{
  SolnData.reset();

  // copy internal variables intVar1 to intVar2
  if(intVarFlag)
    copyElemInternalVariables();

  return 0;
}



int femGrowth::saveSolution()
{
    SolnData.saveSolution();

    if( MIXED_ELEMENT_P0 )
    {
        for(int ee=0; ee<nElem_global; ee++)
        {
            elems[ee]->presDOFprev2 = elems[ee]->presDOFprev;
            elems[ee]->presDOFprev  = elems[ee]->presDOF;
        }
    }

    if(intVarFlag)
    {
        for(int ee=0; ee<nElem_global; ee++)
        {
          elems[ee]->ivar.saveSolution();
        }
    }

    return 0;
}



int femGrowth::copyElemInternalVariables()
{
  // copy internal variables intVar1 to intVar2

  for(int e=0; e<nElem_global; e++)
  {
    elems[e]->ivar.reset();
  }

  return 0;
}




bool femGrowth::converged()
{
  if (rhsNorm < conv_tol)
    return true;

  return false;
}




bool femGrowth::diverging(double factor)
{
  if (localStiffnessError != 0) return true;

  if (isnan(rhsNorm)) return true;

  if (rhsNormPrev > -0.1 && (rhsNorm / rhsNormPrev) > factor) return true;

  return false;
}




void femGrowth::writeNodalData()
{
    //cout << " femGrowth::writeNodalData() " << endl;
    int ii, jj, bb, ee, type, nn, n1, n2, n3, dof;
    double val;

    fout_nodaldata << myTime.cur;

    for(ee=0; ee<NodalDataOutput.size(); ee++)
    {
      type  = (int) (NodalDataOutput[ee][0]);
      nn    = (int) (NodalDataOutput[ee][1] - 1);
      dof   = (int) (NodalDataOutput[ee][2] - 1);

      switch(type)
      {
        case  1 : // total force on all the requested nodes

              val = 0.0;
              //for(ii=0; ii<(NodalDataOutput.size()-3); ii++)
                //val += SolnData.force[nn*ndof+dof];

        break;

        case  2 : // total reaction on all the requested nodes

              dof   = (int) (NodalDataOutput[ee][1] - 1);

              n1 = NodalDataOutput[ee].size()-2;
              val = 0.0;
              for(ii=0; ii<n1; ii++)
              {
                nn  =  (int) (NodalDataOutput[ee][2+ii]-1) ;
                val += SolnData.reac[nn*ndof+dof];
              }

        break;

        case  3 : // displacement

            /*
            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.disp[n1];
              val += 0.25*SolnData.disp[n2];
              val += 0.25*SolnData.disp[n3];
            }
            else
            {
            */
              val = SolnData.disp[nn*ndof+dof];
            //}

        break;

        case  4 : // velocity

            val = SolnData.velo[nn*ndof+dof];

        break;

        case  5 : // acceleration

            val = SolnData.acce[nn*ndof+dof];

        break;

        case  9: // loadfactor for the arc-length method

            val = loadFactor;

        break;

        default :

              cerr << "femSolidmechanics::writeNodalData() ... invalid value of 'type'!" << endl;
        break;
      }

      fout_nodaldata << setw(20) << val;
    }
    fout_nodaldata << endl;

    return;
}





int femGrowth::setSpecifiedDOFs_Displacement(vector<vector<bool> >&  NodeDofType)
{
  if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femGrowth::setSpecifiedDOFs_Displacement() ... STARTED \n");

  vector<int>  nodeNums;
  int ii, jj, nn, ind, dof;

  dofs_specified_disp.clear();


    // set specified DOFs
    for(auto&  bc : boundaryConitions)
    {
        cout << " Patch                 = " << bc->label << endl;
        cout << " BCType                = " << bc->BCType << endl;
        cout << " dof_specified_string  = " << bc->dof_specified_string << endl;
        cout << " dof_specified_int     = " << bc->dof_specified_int << endl;
        cout << " expression            = " << bc->expression << endl;
        cout << " timeFunctionNum       = " << bc->timeFunctionNum << endl;
        cout << endl;  cout << endl;
        //cout << " Patch                 = " << bc->bpt->getLabel() << endl;

        // nodeNums contains new node numbers which are used for the solution
        // nodeCoordsOrig contains coordinates of nodes with numbers before domain decompositions

        nodeNums = bc->bpt->nodeNums;

        if(bc->BCType == "specified")
        {
            dof = bc->dof_specified_int;

            for(ii=0; ii<nodeNums.size(); ++ii)
            {
                nn = nodeNums[ii];

                NodeDofType[nn][dof] = true;
                dofs_specified_disp.push_back(nn*ndof+dof);
            }
        }
        else if(bc->BCType == "fixed")
        {
            for(ii=0; ii<nodeNums.size(); ++ii)
            {
              nn = nodeNums[ii];
              ind = nn*ndof;
              for(dof=0; dof<ndof; ++dof)
              {
                NodeDofType[nn][dof] = true;
                dofs_specified_disp.push_back(ind+dof);
              }
            }
        }
        else
        {
            throw runtime_error("Patch type not available in femGrowth::setSpecifiedDOFs");
        }
    }

    findUnique(dofs_specified_disp);
    printVector(dofs_specified_disp);

    if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femGrowth::setSpecifiedDOFs_Displacement() ... ENDED \n");

    return 0;
}




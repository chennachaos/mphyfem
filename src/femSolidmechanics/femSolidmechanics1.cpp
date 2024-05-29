
#include "femSolidmechanics.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "NewLagrangeElement.h"
#include <algorithm>
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "QuadratureUtil.h"

extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime                myTime;
extern bool debug;

using namespace std;

int ElementBase::elemcount = 0;


femSolidmechanics::femSolidmechanics()
{
    // add new type
}


femSolidmechanics::~femSolidmechanics()
{


}





int femSolidmechanics::prepareInputData()
{
    if(debug) PetscPrintf(PETSC_COMM_WORLD, "\n     femSolidmechanics::prepareInputData()  .... STARTED ...\n");

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


    if( (idd == ELEM_SOLID_2D_TRIA7_MIXED1)    ||
        (idd == ELEM_SOLID_3D_TETRA11_MIXED1)  )
    {
       GeomData.addBubbleNodes(elemConn);
       nNode_global += nElem_global;
    }


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

    if( (idd == ELEM_SOLID_2D_MIXED1)          ||
        (idd == ELEM_SOLID_2D_TRIA7_MIXED1)    ||
        (idd == ELEM_SOLID_2D_TRIA10_MIXED1)   ||
        (idd == ELEM_SOLID_3D_MIXED1)          ||
        (idd == ELEM_SOLID_3D_TETRA11_MIXED1)  ||
        (idd == ELEM_SOLID_3D_TETRA20_MIXED1)  )
    {
      MIXED_ELEMENT = true;

      for(auto& mat : MatlDataList)
        mat->MIXED_ELEMENT = true;
    }

    if( (idd == ELEM_SOLID_2D_MIXED0) || (idd == ELEM_SOLID_3D_MIXED0) )
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

      elems[ee] = NewLagrangeElement(ElementTypeDataList[elemId]->getElemTypeNameNum());

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

    if(debug) PetscPrintf(PETSC_COMM_WORLD, "     femSolidmechanics::prepareInputData()  .... FINISHED ...\n\n");

    return 0;
}



int  femSolidmechanics::prepareDataForMixedElements()
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

  SolnData.presApplied = SolnData.pres;

  return 0;
}





int femSolidmechanics::processForBernsteinElements()
{
    int idd, ii, jj, ee, n1, n2, nn, dof;
    double xx, yy, zz, disp[3];
    vector<int>  vecTemp;

    VectorXd  solnTemp;

    idd = 0;// SolnData.ElemProp[elemConn[0][0]]->id;

    cout << " idd = " << idd << endl;

    // 2D elements
    if(ndim == 2)
    {
      // triangular elements
      if( (idd == 201) || (idd == 202) || (idd == 203) || (idd == 204) || (idd == 211) || (idd == 212))
      {
        // loop over all the elements to
        // a.) find the midnodes, and
        // b.) adjacent nodes to midnodes
        for(ee=0; ee<nElem_global; ee++)
        {
          vecTemp = elems[ee]->nodeNums;
          for(ii=0; ii<vecTemp.size(); ii++)
          {
            if( (ii >= 3) )
            {
              midnodeData[vecTemp[ii]][0] = 1;

              // edge 1-4-2
              if( ii == 3 )
              {
                midnodeData[vecTemp[ii]][1] = vecTemp[0];
                midnodeData[vecTemp[ii]][2] = vecTemp[1];
              }
              // edge 2-5-3
              if( ii == 4 )
              {
                midnodeData[vecTemp[ii]][1] = vecTemp[1];
                midnodeData[vecTemp[ii]][2] = vecTemp[2];
              }
              // edge 3-6-1
              if( ii == 5 )
              {
                midnodeData[vecTemp[ii]][1] = vecTemp[0];
                midnodeData[vecTemp[ii]][2] = vecTemp[2];
              }
            }
          }
        }
      }

      //
      // quadrilateral elements
      if( (idd == 222) || (idd == 223) )
      {
        /*
        // calculate the last node (center node) as the average of all the nodes
        for(ee=0; ee<nElem_global; ee++)
        {
          vecTemp = elems[ee]->nodeNums;

          xx  = (GeomData.NodePosOrig[vecTemp[0]][0] + GeomData.NodePosOrig[vecTemp[1]][0]);
          xx += (GeomData.NodePosOrig[vecTemp[2]][0] + GeomData.NodePosOrig[vecTemp[3]][0]);
          xx += (GeomData.NodePosOrig[vecTemp[4]][0] + GeomData.NodePosOrig[vecTemp[5]][0]);
          xx += (GeomData.NodePosOrig[vecTemp[6]][0] + GeomData.NodePosOrig[vecTemp[7]][0]);

          yy  = (GeomData.NodePosOrig[vecTemp[0]][1] + GeomData.NodePosOrig[vecTemp[1]][1]);
          yy += (GeomData.NodePosOrig[vecTemp[2]][1] + GeomData.NodePosOrig[vecTemp[3]][1]);
          yy += (GeomData.NodePosOrig[vecTemp[4]][1] + GeomData.NodePosOrig[vecTemp[5]][1]);
          yy += (GeomData.NodePosOrig[vecTemp[6]][1] + GeomData.NodePosOrig[vecTemp[7]][1]);

          GeomData.NodePosOrig[vecTemp[8]][0] = xx/8.0;
          GeomData.NodePosOrig[vecTemp[8]][1] = yy/8.0;
          //

          //
          // adjust the center node
          xx  = 0.0625*(GeomData.NodePosOrig[vecTemp[0]][0] + GeomData.NodePosOrig[vecTemp[1]][0] + GeomData.NodePosOrig[vecTemp[2]][0] + GeomData.NodePosOrig[vecTemp[3]][0]);
          xx += 0.125*(GeomData.NodePosOrig[vecTemp[4]][0] + GeomData.NodePosOrig[vecTemp[5]][0] + GeomData.NodePosOrig[vecTemp[6]][0] + GeomData.NodePosOrig[vecTemp[7]][0]);

          yy  = 0.0625*(GeomData.NodePosOrig[vecTemp[0]][1] + GeomData.NodePosOrig[vecTemp[1]][1] + GeomData.NodePosOrig[vecTemp[2]][1] + GeomData.NodePosOrig[vecTemp[3]][1]);
          yy += 0.125*(GeomData.NodePosOrig[vecTemp[4]][1] + GeomData.NodePosOrig[vecTemp[5]][1] + GeomData.NodePosOrig[vecTemp[6]][1] + GeomData.NodePosOrig[vecTemp[7]][1]);

          GeomData.NodePosOrig[vecTemp[8]][0] = 4.0*(GeomData.NodePosOrig[vecTemp[8]][0] - xx);
          GeomData.NodePosOrig[vecTemp[8]][1] = 4.0*(GeomData.NodePosOrig[vecTemp[8]][1] - yy);
          //
        }
        */

        int midnodemapquad9[9][2] = {{0,0}, {0,0}, {0,0}, {0,0},
                                     {0,1}, {1,2}, {2,3}, {3,0}, {0,0}};

        // loop over all the elements to
        // a.) find the midnodes, and
        // b.) adjacent nodes to midnodes
        for(ee=0; ee<nElem_global; ee++)
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
      //

      // loop over the nodes and adjust nodal coordinates
      for(nn=0; nn<nNode_global; nn++)
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

    // 3D elements
    if(ndim == 3)
    {
      // 10-noded tetrahedron element
      if( (idd == 205) || (idd == 206) || (idd == 207) || (idd == 208)  || (idd == 213) || (idd == 214) )
      {
        int midnodemaptet10[10][2] = {{0,0}, {0,0}, {0,0}, {0,0},
                                      {0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3}};

        // loop over all the elements to find the midnodes
        for(ee=0; ee<nElem_global; ee++)
        {
          vecTemp = elems[ee]->nodeNums;

          for(ii=4; ii <=9; ii++)
          {
            nn = vecTemp[ii];

            midnodeData[nn][0] = 1;
            midnodeData[nn][1] = vecTemp[midnodemaptet10[ii][0]];
            midnodeData[nn][2] = vecTemp[midnodemaptet10[ii][1]];
          }
        }
      }

      // 18-noded wedge element
      if( (idd == 224) || (idd == 225) )
      {
        int midnodemapwedge18[15][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0},
                                        {0,1}, {0,2}, {0,3}, {1,2}, {1,4}, {2,5}, {3,4}, {3,5}, {4,5}};

        // loop over all the elements to find the midnodes
        for(ee=0; ee<nElem_global; ee++)
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
      if( (idd == 226) || (idd == 227) )
      {
        int midnodemaphex27[20][2] = {{0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0}, {0,0},
                                      {0,1}, {0,3}, {0,4}, {1,2}, {1,5}, {2,3}, {2,6}, {3,7}, {4,5}, {4,7}, {5,6}, {6,7}};

        // loop over all the elements to find the midnodes
        for(ee=0; ee<nElem_global; ee++)
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

      //
      // loop over the nodes and adjust nodal coordinates
      for(nn=0; nn<nNode_global; nn++)
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
    return 0;
}



bool femSolidmechanics::converged()
{
  if (rhsNorm < conv_tol)
    return true;

  return false;
}




bool femSolidmechanics::diverging(double factor)
{
  if (localStiffnessError != 0) return true;

  if (isnan(rhsNorm)) return true;

  if (rhsNormPrev > -0.1 && (rhsNorm / rhsNormPrev) > factor) return true;

  return false;
}




int femSolidmechanics::printComputerTime(bool reset, int detailFlg)
{
  cout << "----------------------------------------------------" << endl;

  //printf("femSolidmechanics::calcStiffnessRes:  %7.3f sec ->%5.1f %\n", ctimCalcStiffRes, ctimCalcStiffRes*100.);
  //printf("femSolidmechanics::factoriseSolvUpdt: %7.3f sec ->%5.1f %\n", ctimFactSolvUpdt, ctimFactSolvUpdt*100.);

  if (reset)
  {
    ctimFactSolvUpdt = 0.;
    ctimCalcStiffRes = 0.;
  }

  return 0;
}



int femSolidmechanics::setInitialConditions()
{
    PetscPrintf(MPI_COMM_WORLD, "femSolidmechanics::setInitialConditions() ... \n");
    if(debug) cout <<  " femSolidmechanics::setInitialConditions() ... STARTED " <<  endl;

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
    PetscPrintf(MPI_COMM_WORLD, "femSolidmechanics::setInitialConditions() ... \n");

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
        if( midnodeData[nn][0] ) // if midnode
        {
            xx  = 0.50*GeomData.NodePosOrig[nn][0];
            xx += 0.25*GeomData.NodePosOrig[midnodeData[nn][1]][0];
            xx += 0.25*GeomData.NodePosOrig[midnodeData[nn][2]][0];

            yy  = 0.50*GeomData.NodePosOrig[nn][1];
            yy += 0.25*GeomData.NodePosOrig[midnodeData[nn][1]][1];
            yy += 0.25*GeomData.NodePosOrig[midnodeData[nn][2]][1];

            if(ndim == 3)
            {
              zz  = 0.50*GeomData.NodePosOrig[nn][2];
              zz += 0.25*GeomData.NodePosOrig[midnodeData[nn][1]][2];
              zz += 0.25*GeomData.NodePosOrig[midnodeData[nn][2]][2];
            }
        }
        else
        {
            xx = GeomData.NodePosOrig[nn][0];
            yy = GeomData.NodePosOrig[nn][1];

            if(ndim == 3)
              zz = GeomData.NodePosOrig[nn][2];
        }

        ind = nn*ndof;

      // square prism torsion only - Gil
      // (vx,vy,vz) = 1500 sin(pi*y/12) (z,0,-x) cm/s  (centimetre per second)
      fact = 1.5*sin(PI*yy/12.0); // cm/ms (centimetre per millisecond)
      //fact = 100.0*sin(PI*yy/12.0);
      //SolnData.velo[ind]    =  fact*zz;
      //SolnData.velo[ind+1]  =  0.0;
      //SolnData.velo[ind+2]  = -fact*xx;

      // triangular prism torsion only - Scovazzi
      //fact = 40.0*sin(PI*zz/12.0);
      //SolnData.velo[ind]    =  fact*yy;
      //SolnData.velo[ind+1]  = -fact*xx;
      //SolnData.velo[ind+2]  =  0.0;

      // triangular prism torsion and compression - Scovazzi
      //fact = zz;
      //SolnData.velo[ind]    =  fact*yy*0.5;
      //SolnData.velo[ind+1]  = -fact*xx*0.5;
      //SolnData.velo[ind+2]  = -fact*5.0;

      //SolnData.velo[ind]    =  5.0*yy/3.0;
      //SolnData.velo[ind+1]  =  0.0;
      //SolnData.velo[ind+2]  =  0.0;

      //SolnData.velo[ind]    =    0.0;
      //SolnData.velo[ind+1]  =    0.0;
      //SolnData.velo[ind+2]  = -227.0;
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

    if(debug) cout <<  " femSolidmechanics::setInitialConditions() ... ENDED " <<  endl;

    return 0;
}





int femSolidmechanics::setBoundaryConditions()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femSolidmechanics::setBoundaryConditions() ... STARTED \n");}

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
            throw runtime_error("Patch type not available in femSolidmechanics::setBoundaryConditions");
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

    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femSolidmechanics::setBoundaryConditions() ... ENDED \n");}

    return 0;
}






int femSolidmechanics::addExternalForces()
{
    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femSolidmechanics::addExternalForces() ... STARTED \n");}

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

                    cout << "normals" << endl;
                    cout <<   normal[0] << '\t' <<   normal[1] << '\t' <<   normal[2] << endl;
                    cout << tangent1[0] << '\t' << tangent1[1] << '\t' << tangent1[2] << endl;
                    cout << tangent2[0] << '\t' << tangent2[1] << '\t' << tangent2[2] << endl; cout << endl;

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

    if(debug) {cout << endl; cout << endl; cout <<  " femSolidmechanics::addExternalForces() ... ENDED " <<  endl;}

    return 0;
}






int femSolidmechanics::timeUpdate()
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
    if(intVarFlag)
      copyElemInternalVariables();

    // update time functions
    for(auto& tmf : timeFunctions)
      tmf->update();

    SolnData.saveSolution();

    // set Dirichlet boundary conditions
    setBoundaryConditions();


    if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n femSolids::timeUpdate() ... FINISHED \n\n"); }

    return 0;
}



int femSolidmechanics::updateIterStep()
{
    if(debug)  cout << " femSolidmechanics::updateIterStep() ... " << endl;

    SolnData.updateIterStep();

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

    if(debug)  cout << " femSolidmechanics::updateIterStep() ... " << endl;

    return 0;
}


int femSolidmechanics::reset()
{
  SolnData.reset();

  // copy internal variables intVar1 to intVar2
  //if(intVarFlag)
    //copyElemInternalVariables();

  return 0;
}



int femSolidmechanics::saveSolution()
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



int femSolidmechanics::copyElemInternalVariables()
{
  // copy internal variables intVar1 to intVar2

  for(int e=0; e<nElem_global; e++)
  {
    elems[e]->ivar.reset();
    //elems[e]->intVar = elems[e]->intVarPrev;
  }

  return 0;
}




int femSolidmechanics::writeNodalData()
{
    int ii, jj, bb, ee, type, nn, n1, n2, n3, dof;
    double val;

    for(ee=0; ee<OutputData.size(); ee++)
    {
      type  = (int) (OutputData[ee][0]);
      nn    = (int) (OutputData[ee][1] - 1);
      dof   = (int) (OutputData[ee][2] - 1);

      switch(type)
      {
        case  1 : // total force on all the requested nodes

              val = 0.0;
              for(ii=0; ii<(OutputData.size()-3); ii++)
                val += SolnData.force[OutputData[ee][3+ii]*ndof+dof];

        break;

        case  2 : // total reaction on all the requested nodes

              dof   = (int) (OutputData[ee][1] - 1);

              n1 = OutputData[ee].size()-2;
              val = 0.0;
              for(ii=0; ii<n1; ii++)
              {
                nn  =  (int) (OutputData[ee][2+ii]-1) ;
                val += SolnData.reac[nn*ndof+dof];
              }

        break;

        case  3 : // displacement

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
              val = SolnData.disp[nn*ndof+dof];
            }

        break;

        case  4 : // velocity

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.velo[n1];
              val += 0.25*SolnData.velo[n2];
              val += 0.25*SolnData.velo[n3];
            }
            else
            {
              val = SolnData.velo[nn*ndof+dof];
            }

        break;

        case  5 : // acceleration

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.acce[n1];
              val += 0.25*SolnData.acce[n2];
              val += 0.25*SolnData.acce[n3];
            }
            else
            {
              val = SolnData.acce[nn*ndof+dof];
            }

        break;

        case  9: // loadfactor for the arc-length method

            val = loadFactor;

        break;

        default : 

              cerr << "femSolidmechanics::writeNodalData() ... invalid value of 'type'!" << endl;
        break;
      }

      char      tmp[200];
      string    tmpStr;

      sprintf(tmp," \t %12.6E", val);

      tmpStr.append(tmp);

      //prgWriteToTFile(tmpStr);
    }

    //printData(0, 0);

    return 0;
}



int femSolidmechanics::writeNodalDataExplicitScheme()
{
    int ii, jj, bb, ee, type, nn, n1, n2, n3, dof;
    double val, fact=1.0;

    fout_explicit << myTime.cur << '\t' << timeFunctions[0]->getValue();

    for(ee=0; ee<OutputData.size(); ee++)
    {
      type  = (int) (OutputData[ee][0]);
      nn    = (int) (OutputData[ee][1] - 1);
      dof   = (int) (OutputData[ee][2] - 1);

      switch(type)
      {
        case  1 : // total force on all the requested nodes

              val = 0.0;
              for(ii=0; ii<(OutputData.size()-3); ii++)
                val += SolnData.force[OutputData[ee][3+ii]*ndof+dof];

        break;

        case  2 : // total reaction on all the requested nodes

              dof   = (int) (OutputData[ee][1] - 1);

              n1 = OutputData[ee].size()-2;
              val = 0.0;
              for(ii=0; ii<n1; ii++)
              {
                nn  =  (int) (OutputData[ee][2+ii]-1) ;
                val += SolnData.reac[nn*ndof+dof];
              }

        break;

        case  3 : // displacement

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
              val = SolnData.disp[nn*ndof+dof];
            }

        break;

        case  4 : // velocity

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.velo[n1];
              val += 0.25*SolnData.velo[n2];
              val += 0.25*SolnData.velo[n3];
            }
            else
            {
              val = SolnData.velo[nn*ndof+dof];
            }

        break;

        case  5 : // acceleration

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.acce[n1];
              val += 0.25*SolnData.acce[n2];
              val += 0.25*SolnData.acce[n3];
            }
            else
            {
              val = SolnData.acce[nn*ndof+dof];
            }

        break;

        default : 

               cerr << "femSolidmechanics::writeNodalData() ... invalid value of 'type'!" << endl;
        break;
      }

      fout_explicit << '\t' << val ;
    }
    fout_explicit << endl;

    return 0;
}



int femSolidmechanics::setSpecifiedDOFs_Displacement(vector<vector<bool> >&  NodeDofType)
{
  if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femSolidmechanics::setSpecifiedDOFs_Displacement() ... STARTED \n");

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
            throw runtime_error("Patch type not available in femSolidmechanics::setSpecifiedDOFs");
        }
    }

    findUnique(dofs_specified_disp);
    printVector(dofs_specified_disp);

    if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femSolidmechanics::setSpecifiedDOFs_Displacement() ... ENDED \n");

    return 0;
}











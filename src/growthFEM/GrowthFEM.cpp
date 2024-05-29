
#include "GrowthFEM.h"
#include "DataBlockTemplate.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "NewElementGrowthFEM.h"
#include "NewMaterial.h"
#include "List.h"
#include <algorithm>
#include "utilitiesmaterial.h"


extern List<TimeFunction> timeFunction;
extern MyTime           myTime;

using namespace std;

int ElementBase::elemcount = 0;


GrowthFEM::GrowthFEM()
{
}




GrowthFEM::~GrowthFEM()
{

}


void GrowthFEM::printLogo()
{
  cout << "\n\n\n";
  cout << "   +----------------------------------------------------+  \n";
  cout << "   | ************************************************** |  \n";
  cout << "   | **                                              ** |  \n";
  cout << "   | **        MULTI-PHYSICS ANALYSIS PROGRAM        ** |  \n";
  cout << "   | **      ==================================      ** |  \n";
  cout << "   | **        Growth Problems                       ** |  \n";
  cout << "   | **                                              ** |  \n";
  cout << "   | **                                              ** |  \n";
  cout << "   | **          Chennakesava Kadapa 2022            ** |  \n";
  cout << "   | **                                                 |  \n";
  cout << "   | ************************************************** |  \n";
  cout << "   +----------------------------------------------------|  \n";
  cout << "                                                           \n";

  return;
}







void  GrowthFEM::readInputData(std::ifstream &Ifile, MyString &line)
{
/*
  char fct[] = "GrowthFEM::readInputData";

  MyString tmpl, *word;

  char tmp[30];

  int nw, i, j, k, n, nn, ii, bb;
  
  double fact;

  MyStringList   sTmp;
  List<Vector<int> > lviTmp;
  List<Vector<double> > lvdTmp;
  vector<double>  dblTemp;

  DataBlockTemplate t1, t2;

  growthfem->key.addNew("element type", //0
                          "material type", //1
                          "nodes", //2
                          "elements", //3
                          "prescribed boundary conditions", //4
                          "nodal forces", //5
                          "element face loads", //6
                          "nodal data output", //7
                          "control parameters", //8
                          "initial conditions", //9
                          "face load elements"); //10

  switch (domain[GROWTHFEM].key.whichBegins(line))
  {
    case  0: //cout << "     GrowthFEM: reading 'element type' ...\n\n";

            //SolnData.ElemProp.add(new PropertyItem(ELEMENTTYPE));
            //SolnData.ElemProp[SolnData.ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            break;

    case  1: //cout << "     GrowthFEM: reading 'material type' ...\n\n";

            //SolnData.MatlProp.add(new PropertyItem(MATERIAL));
            //SolnData.MatlProp[SolnData.MatlProp.n-1].readInputData(Ifile,line,"input error in 'material type'!");

            break;

    case  2: //cout << "     GrowthFEM: reading 'nodes' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1, fct, "invalid specification of 'nodes' !");

            nodePosData.resize(lvdTmp.n);

            for(i=0; i<lvdTmp.n; i++)
            {
              for(j=0;j<(lvdTmp[i].n-1);j++)
                nodePosData[i][j] = lvdTmp[i][j+1];
            }

            break;

    case  3: //cout << "     GrowthFEM: reading 'elements' ...\n\n";

            if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
              prgError(1,fct,"invalid input in 'elements'!");

            elemConn.resize(lviTmp.n);

            for(i=0;i<lviTmp.n;i++)
            {
              //cout << lviTmp[i] << endl;
              elemConn[i].resize(lviTmp[i].n-1 );
              if(lviTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'elements' !");

              for(j=1;j<lviTmp[i].n;j++)
                elemConn[i][j-1] = lviTmp[i][j] - 1;
             }

            break;

    case  4: //cout << "     GrowthFEM: reading 'prescribed boundary conditions' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'prescribed boundary conditions'!");

            for(i=0;i<lvdTmp.n;i++)
            {
              if(lvdTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'prescribed boundary conditions' !");

              dblTemp.resize(lvdTmp[i].n);

              dblTemp[0] = lvdTmp[i][0]-1;
              dblTemp[1] = lvdTmp[i][1]-1;
              dblTemp[2] = lvdTmp[i][2];

              if(int(dblTemp[1]) < ndof)
              {
                DirichletBCs.push_back(dblTemp);
              }
              else
              {
                DirichletBCs_Pres.push_back(dblTemp);
              }
            }

            break;

    case  5: //cout << "     GrowthFEM: reading 'nodal forces' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'nodal forces'!");

             nodeForcesData.resize(lvdTmp.n);

             for(i=0;i<lvdTmp.n;i++)
             {
                nodeForcesData[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'nodal forces' !");

                for(j=0;j<lvdTmp[i].n;j++)
                  nodeForcesData[i][j] = lvdTmp[i][j];
             }

             break;

    case  6: //cout << "     GrowthFEM: reading 'element face loads' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'element face loads'!");

            NeumannBCs.resize(lvdTmp.n);

            for(i=0;i<lvdTmp.n;i++)
            {
              if(lvdTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'element face loads' !");

              NeumannBCs[i].resize(lvdTmp[i].n);

              NeumannBCs[i][0] = lvdTmp[i][0]-1; // element number
              NeumannBCs[i][1] = lvdTmp[i][1]-1; // edge/face number
              NeumannBCs[i][2] = lvdTmp[i][2]-1; // direction
              NeumannBCs[i][3] = lvdTmp[i][3];   // value
            }

            break;

    case  7: //cout << "     GrowthFEM: reading 'nodal data output' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'nodal data output'!");

             OutputData.resize(lvdTmp.n);

             for(i=0;i<lvdTmp.n;i++)
             {
                OutputData[i].resize(lvdTmp[i].n);
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'nodal data output' !");

                for(j=0;j<lvdTmp[i].n;j++)
                  OutputData[i][j] = lvdTmp[i][j];
             }

             //printVector(OutputData);

             break;

    case   8: //cout << "     GrowthFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( GrowthFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            rhoInfty = lvdTmp[0][2];

            if(tis < 100)
              IMPLICIT_SOLVER = true;
            else
              IMPLICIT_SOLVER = false;

            break;

    case  9: //cout << "     GrowthFEM: reading 'initial conditions' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'nodal data output'!");

             InitialConds.resize(lvdTmp.n);

             for(i=0;i<lvdTmp.n;i++)
             {
                if(lvdTmp[i].n < 1)
                   prgError(2, fct, "invalid number of 'nodal data output' !");

                InitialConds[i].resize(lvdTmp[i].n);
                for(j=0;j<lvdTmp[i].n;j++)
                  InitialConds[i][j] = lvdTmp[i][j];
             }

             //printVector(InitialConds);

             break;

    case  10: cout << "     GrowthFEM: reading 'face load elements' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'face load elements'!");

            ElemFaceLoadData.resize(lvdTmp.n);

            for(i=0; i<lvdTmp.n; i++)
            {
              if(lvdTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'face load elements' !");

              ElemFaceLoadData[i].resize(lvdTmp[i].n-1);

              for(j=0; j<lvdTmp[i].n-1; j++)
                ElemFaceLoadData[i][j] = lvdTmp[i][1+j];

              for(j=5; j<lvdTmp[i].n-1; j++)
                ElemFaceLoadData[i][j] = ElemFaceLoadData[i][j]-1;
            }

            cout << " aaaaaaaaaaaaa " << lvdTmp.n << endl;

            break;


    case -1: // go and inherit from DOMAIN

            this->Domain::readInputData(Ifile,line);

            break;
  }
*/
  return;
}





void GrowthFEM::prepareElemProp()
{
    if(SolnData.ElemProp.size() < 1)
      throw  runtime_error("\n\n Element property data missing! \n\n");


    for(int ii=0; ii<SolnData.ElemProp.size(); ii++)
    {
      string matkey(SolnData.ElemProp[ii]->name);

      SolnData.ElemProp[ii]->id = getElementID_Growth(matkey);

      cout << "SolnData.ElemProp[ii]->id = " << SolnData.ElemProp[ii]->id << endl;
    }

    return;
}







void GrowthFEM::prepareMatlProp()
{
    if(SolnData.MatlProp.size() < 1)
      throw  runtime_error("\n\n Material property data missing! \n\n");

    intVarFlag = false;

    for(int ii=0; ii<SolnData.MatlProp.size(); ii++)
    {
      // assign correct id of the material type (as stored in the database) based on the material name
      cout <<  SolnData.MatlProp[ii]->name <<  endl;

      std::string matkey(SolnData.MatlProp[ii]->name);

      int  idd = getMaterialID(matkey);

      cout << " idd = " << idd << endl;

      SolnData.MatlProp[ii]->id = idd;

      if( ( (idd >= 101) && (idd <= 106)) || (idd == 2006) || (idd == 2007) )
        intVarFlag = true;

      MatlDataList.push_back( NewMaterial(idd) );

      MatlDataList[ii]->matData.resize( SolnData.MatlProp[ii]->data.size() );

      MatlDataList[ii]->SolnData = &(SolnData);

      for(int jj = 0; jj<SolnData.MatlProp[ii]->data.size(); jj++)
        MatlDataList[ii]->matData[jj] = SolnData.MatlProp[ii]->data[jj];
    }

    numMaterials = MatlDataList.size();

    return;
}





void GrowthFEM::prepareInputData()
{
    //printf("\n     GrowthFEM::prepareInputData()  .... STARTED ...\n");

    int ii, jj, kk, ee, aa, bb, cc, nn, ind, count, gp, r;

    assert(ndim > 0 && ndim < 4);

    //cout << " ndim  = " << ndim << endl;
    //cout << " ndof = " << ndof << endl;

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();

    // Prepare patchElemProp data. Assign suitable id(in the database) based on the element name
    if(SolnData.ElemProp.size() > 0)
      prepareElemProp();

    // Prepare patchMatlProp data. Assign suitable id(in the database) based on the material name
    //if(SolnData.MatlProp.size() > 0)
      //prepareMatlProp();

    numMaterials = MatlDataList.size();

    cout << " numMaterials  = " << numMaterials  << endl;

    int idd = SolnData.ElemProp[elemConn[0][0]]->id;

    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    nNode = nodePosData.size();
    nElem = elemConn.size();
    npElem = elemConn[0].size()-3;

    elems = new ElementBase* [nElem];

    for(ee=0;ee<nElem;ee++)
    {
      elems[ee] = NewElementGrowthFEM(SolnData.ElemProp[elemConn[ee][0]]->id);

      elems[ee]->elenum   =  ee;

      elems[ee]->elmType  =  elemConn[ee][0];
      elems[ee]->matType  =  elemConn[ee][1];
      elems[ee]->secType  =  elemConn[ee][2];

      vector<int>  vecTemp;
      for(ii=0;ii<elemConn[ee].size()-3;ii++)
        vecTemp.push_back(elemConn[ee][3+ii]);

      elems[ee]->nodeNums = vecTemp;

      elems[ee]->SolnData = &(SolnData);
      elems[ee]->GeomData = &(GeomData);
      elems[ee]->MatlData = MatlDataList[elemConn[ee][1]];
    }
    cout << " AAAAAAAAAAAAAAAAA " << endl;

    // add Neumann/Traction BC values to the elements
    vector<double> vecTemp(3);
    for(ii=0; ii<NeumannBCs.size(); ii++)
    {
      ee  = (int) NeumannBCs[ii][0]; // element number

      vecTemp[0] = NeumannBCs[ii][1]; // edge/face number
      vecTemp[1] = NeumannBCs[ii][2]; // element number
      vecTemp[2] = NeumannBCs[ii][3]; // specified value

      elems[ee]->NeumannBCs.push_back(vecTemp);
    }

    // add face load elements
    nElemFaces = elemConnFaces.size();
    nElemFaces = ElemFaceLoadData.size();

    elemsFaces = new ElementBase* [nElemFaces];
    for(ee=0; ee<nElemFaces; ee++)
    {
      //elemsFaces[ee] = NewElement(elemConnFaces[ee][0]);
      cout << "idd = " << SolnData.ElemProp[elemConnFaces[ee][0]]->id <<  endl;

      elemsFaces[ee] = NewElementGrowthFEM(SolnData.ElemProp[elemConnFaces[ee][0]]->id);

      elemsFaces[ee]->elenum   =  ee;

      elemsFaces[ee]->elmType  =  elemConnFaces[ee][0];
      elemsFaces[ee]->matType  =  elemConnFaces[ee][1];
      elemsFaces[ee]->secType  =  elemConnFaces[ee][2];

      vector<int>  vecTemp;
      for(ii=0;ii<elemConnFaces[ee].size()-3;ii++)
        vecTemp.push_back(elemConnFaces[ee][3+ii]);

      elemsFaces[ee]->nodeNums = vecTemp;

      elemsFaces[ee]->SolnData = &(SolnData);
      elemsFaces[ee]->GeomData = &(GeomData);
      //elemsFaces[ee]->MatlData = MatlData;

      elemsFaces[ee]->MatlData = MatlDataList[elemConnFaces[ee][1]];

      elemsFaces[ee]->prepareElemData();
    }


    elemsFaces = new ElementBase* [nElemFaces];
    for(ee=0; ee<nElemFaces; ee++)
    {
      elemsFaces[ee] = NewElementGrowthFEM(ElemFaceLoadData[ee][0]);

      nn  = (int)  ElemFaceLoadData[ee][1];
      elemsFaces[ee]->tracForces[0] = ElemFaceLoadData[ee][2];
      elemsFaces[ee]->tracForces[1] = ElemFaceLoadData[ee][3];
      elemsFaces[ee]->tracForces[2] = ElemFaceLoadData[ee][4];

      vector<int>  vecTemp(nn);
      for(ii=0; ii<nn; ii++)
      {
        vecTemp[ii] = (int) ElemFaceLoadData[ee][5+ii];

      }
      elemsFaces[ee]->nodeNums = vecTemp;

      elemsFaces[ee]->SolnData = &(SolnData);
      elemsFaces[ee]->GeomData = &(GeomData);
    }

    ///////////////////////////////////////////////////////////////////
    //
    // set GeomData details
    //
    ///////////////////////////////////////////////////////////////////

    GeomData.setDimension(ndim);
    GeomData.setNdof(ndof);

    GeomData.build();
    GeomData.setNodalPositions(nodePosData);

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    SolnData.initialise(nNode*ndof, nNode, 0, 0);

    // node-to-element connectivity for post-processing

    vector<int>  vecTemp2;

    node_elem_conn.resize(nNode);

    for(ee=0;ee<nElem;ee++)
    {
      vecTemp2 = elems[ee]->nodeNums;
      for(ii=0; ii<vecTemp2.size(); ii++)
      {
        node_elem_conn[vecTemp2[ii]].push_back(ee);
      }
    }

    // find the unique values and sort the entries
    for(ee=0;ee<nNode;ee++)
    {
      findUnique(node_elem_conn[ee]);
    }

    // create the arrays for proce midnodes
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

    cout << " idd = " << idd << endl;

    // adjust mid-nodes for quadratic Bezier triangle/tetrahedron element
    if( (idd==6003) || (idd==6004) || (idd==6005) || (idd==6053) || (idd==6054) )
    {
      processForBernsteinElements();
    }

    MIXED_ELEMENT = false;
    for(ii=0; ii<MatlDataList.size(); ii++)
      MatlDataList[ii]->MIXED_ELEMENT = false;

    if( (idd == 6003) || (idd==6004) || (idd == 6053) || (idd == 6054) || (idd == 6055) )
    {
      MIXED_ELEMENT = true;
    }

    if( (idd == 6002) || (idd == 6003) || (idd==6004) || (idd==6005) || (idd == 6052) || (idd == 6053) || (idd == 6054)  || (idd == 6055))
    {
      for(ii=0; ii<MatlDataList.size(); ii++)
        MatlDataList[ii]->MIXED_ELEMENT = true;
    }

    cout << " idd  = " << idd << endl;
    cout << " ndim  = " << ndim << endl;
    cout << " ndof = " << ndof << endl;
    cout << " MIXED_ELEMENT = " << MIXED_ELEMENT << endl;


    SolnData.TRULY_INCOMPRESSIBLE = false;
    for(ii=0; ii<MatlDataList.size(); ii++)
    {
      if(CompareDoubles(MatlDataList[ii]->getKinv(),0.0))
        SolnData.TRULY_INCOMPRESSIBLE = true;
    }


    if(nodePosDataTarget.size() > 0)
    {
        SolnData.dispTarget.resize(nNode*3);
        int r1, r2, r3;
        for(int ii=0; ii<nNode; ii++)
        {
            r1 = 3*ii;
            r2 = r1+1;
            r3 = r2+1;

            SolnData.dispTarget[r1] = nodePosDataTarget[ii][0] - nodePosData[ii][0];
            SolnData.dispTarget[r2] = nodePosDataTarget[ii][1] - nodePosData[ii][1];
            SolnData.dispTarget[r3] = nodePosDataTarget[ii][2] - nodePosData[ii][2];
        }
    }
    
    
    if(DEBUG) {cout <<  " GrowthFEM::prepareInputData() ... FINISHED \n\n" <<  endl;}

    return;
}



//
void  GrowthFEM::prepareDataForPressure()
{
    // this is the subroutine for the assumption that pressure is dis-continuous
    // across material interfaces

    SolnData.var2.resize(max(nElem, nNode));
    SolnData.var2.setZero();

    SolnData.var2applied = SolnData.var2;

    // gather all the pressure nodes (dofs)
    //

    vector<vector<int> >  pressure_nodes(numMaterials);
    //pressure_nodes.resize(numMaterials);

    int  ee, mm, ii, jj, matNum;
    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;
    vector<int>  vecIntTemp, nNode_Pres_Vec(numMaterials, 0);

    cout << " prepareDataForPressure ... idd = " << idd << endl;

    if( idd == 6003 ) // Bezier Quad9/4
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
    if( idd==6004 ) // Bezier Tria6/3
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
    else if(idd == 6053)   // Bezier Wedge18/6
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
    else if(idd == 6054)   // Bezier Hex27/8
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
    else if(idd == 6055)   // Bezier Tet10/4
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
    else
    {
      cerr << " Wrong element type in 'GrowthFEM::prepareDataForPressure()' " << endl;
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
//




/*
void  GrowthFEM::prepareDataForPressure()
{
    // this is the subroutine for the assumption that pressure is continuous
    // across material interfaces

    SolnData.var2.resize(max(nElem, nNode));
    SolnData.var2.setZero();

    SolnData.var2applied = SolnData.var2;

    // gather all the pressure nodes (dofs)
    //
    vector<int>  pressure_nodes;

    int  ee, ii, jj;
    int  idd = SolnData.ElemProp[elemConn[0][0]].id;

    cout << " prepareDataForPressure ... idd = " << idd << endl;

    if( (idd == 6003) || (idd==6004) ) // Bezier Quad9/4
    {
      for(ee=0; ee<nElem; ++ee)
      {
        pressure_nodes.push_back(elemConn[ee][3]);
        pressure_nodes.push_back(elemConn[ee][4]);
        pressure_nodes.push_back(elemConn[ee][5]);
        pressure_nodes.push_back(elemConn[ee][6]);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][3]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][4]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][5]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][6]);
      }
    }
    else if(idd == 6053)   // Bezier Wedge18/6
    {
      for(ee=0; ee<nElem; ++ee)
      {
        pressure_nodes.push_back(elemConn[ee][3]);
        pressure_nodes.push_back(elemConn[ee][4]);
        pressure_nodes.push_back(elemConn[ee][5]);
        pressure_nodes.push_back(elemConn[ee][6]);
        pressure_nodes.push_back(elemConn[ee][7]);
        pressure_nodes.push_back(elemConn[ee][8]);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][3]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][4]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][5]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][6]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][7]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][8]);
      }
    }
    else if(idd == 6054)   // Bezier Hex27/8
    {
      for(ee=0; ee<nElem; ++ee)
      {
        pressure_nodes.push_back(elemConn[ee][3]);
        pressure_nodes.push_back(elemConn[ee][4]);
        pressure_nodes.push_back(elemConn[ee][5]);
        pressure_nodes.push_back(elemConn[ee][6]);
        pressure_nodes.push_back(elemConn[ee][7]);
        pressure_nodes.push_back(elemConn[ee][8]);
        pressure_nodes.push_back(elemConn[ee][9]);
        pressure_nodes.push_back(elemConn[ee][10]);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][3]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][4]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][5]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][6]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][7]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][8]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][9]);
        elems[ee]->forAssyVecPres.push_back(elemConn[ee][10]);
      }
    }
    else
    {
      cerr << " Wrong element type in 'GrowthFEM::prepareDataForPressure()' " << endl;
      //exit(999);
    }

    //printVector(pressure_nodes);

    findUnique(pressure_nodes);
    nNode_Pres = pressure_nodes.size();
    //printVector(pressure_nodes);

    cout << " nNode_Pres = " << nNode_Pres << endl;

    vector<bool>  NodeType_pres(nNode,true); // all fixed
    vector<int>   ID_pres(nNode,-1);

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
      //cout << ii << '\t' << NodeType_pres[ii] << endl;

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

    pressure_nodes_map_g2l.resize(nNode,-1);
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
    for(ii=0;ii<DirichletBCs_Pres.size();ii++)
    {
      SolnData.var2applied(pressure_nodes_map_g2l[DirichletBCs_Pres[ii][0]]) = DirichletBCs_Pres[ii][2];
    }

    //printVector(SolnData.var2applied);

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

    return;
}
*/




void GrowthFEM::processForBernsteinElements()
{
    int idd, ii, jj, ee, n1, n2, nn, dof;
    double xx, yy, zz, disp[3];
    vector<int>  vecTemp;

    VectorXd  solnTemp;

    idd = SolnData.ElemProp[elemConn[0][0]]->id;

    cout << " processForBernsteinElements ... idd = " << idd << '\t' << ndim << endl;

    // 2D elements
    if(ndim == 2)
    {
      // triangular element
      if( idd == 6004 )
      {
        int midnodemapquad9[6][2] = {{0,0}, {0,0}, {0,0},
                                     {0,1}, {1,2}, {2,0}};

        // loop over all the elements to
        // a.) find the midnodes, and
        // b.) adjacent nodes to midnodes
        for(ee=0; ee<nElem; ee++)
        {
          vecTemp = elems[ee]->nodeNums;
          for(ii=3; ii<=5; ii++)
          {
            nn = vecTemp[ii];

            midnodeData[nn][0] = 1;
            midnodeData[nn][1] = vecTemp[midnodemapquad9[ii][0]];
            midnodeData[nn][2] = vecTemp[midnodemapquad9[ii][1]];
          }
        }
      }

      // quadrilateral element
      if( idd == 6003 )
      {
        int midnodemapquad9[9][2] = {{0,0}, {0,0}, {0,0}, {0,0},
                                     {0,1}, {1,2}, {2,3}, {3,0}, {0,0}};

        // loop over all the elements to
        // a.) find the midnodes, and
        // b.) adjacent nodes to midnodes
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

    // 3D elements
    if(ndim == 3)
    {
      // 10-noded tetrahedron element
      if( (idd == 205) || (idd == 206) || (idd == 207) || (idd == 208)  || (idd == 213) || (idd == 214) )
      {
        int midnodemaptet10[10][2] = {{0,0}, {0,0}, {0,0}, {0,0},
                                      {0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3}};

        // loop over all the elements to find the midnodes
        for(ee=0; ee<nElem; ee++)
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
      if( idd == 6053 )
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
      if( idd == 6054 )
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

      //
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
    return;
}





bool GrowthFEM::converged()
{
  char fct[] = "GrowthFEM::converged";


  if (rNorm < tol && localStiffnessError == 0)
  {
    return true;
  }

  return false;
}




bool GrowthFEM::diverging(double factor)
{
  if (rNormPrev > -0.1 && (rNorm / rNormPrev) > factor) return true;

  if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;

  return false;
}



void GrowthFEM::setInitialConditions()
{
    cout <<  " GrowthFEM::setInitialConditions() ... STARTED" <<  endl;

    SolnData.var1Dot.setZero();

    SolnData.var1.setZero();
    double  ii, jj, ind;
    double  xx,  yy, rad;
    for(ii=0; ii<nNode; ii++)
    {
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
        xx = GeomData.NodePosOrig[ii][0];
        yy = GeomData.NodePosOrig[ii][1];
      }

      rad = sqrt(xx*xx+yy*yy);

      //SolnData.var1[ii*ndof+1] = -xx*0.05;
      //SolnData.var1[ii*ndof+2] = rad*rad*0.01;
      SolnData.var1[ii*ndof+1] = -0.02*sin(PI*xx);
    }
    updateIterStep();

    cout <<  " GrowthFEM::setInitialConditions() ... ENDED" <<  endl;

    return;
}


void GrowthFEM::assignBoundaryConditions()
{
    double timeFact=1.0, specVal;
    timeFact = timeFunction[0].prop;

    int idd, aa, ee, ii, jj, n1, n2, nn, dof, ind;

    double  xx=0.0, yy=0.0, zz=0.0, theta, disp[3];

    VectorXd  solnTemp(nNode*ndof);

    double  tCur = myTime.cur;

    //cout << "SolnData.var1applied.rows() =" << SolnData.var1applied.rows() << endl;
    SolnData.var1applied.setZero();
    solnTemp.setZero();

    for(ii=0; ii<DirichletBCs.size(); ii++)
    {
      nn       =  (int) (DirichletBCs[ii][0]);
      dof      =  (int) (DirichletBCs[ii][1]);

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

      ind = nn*ndof+dof;

      specVal  =  DirichletBCs[ii][2]*timeFact;

      solnTemp[ind] =  specVal;

      SolnData.var1applied[ind] = specVal;

      //cout << nn << '\t' << xx << '\t' << yy << '\t' << zz << '\t' <<  specVal << endl;
    }

    //printVector(SolnData.var1);

    for(ii=0; ii<DirichletBCs.size(); ii++)
    {
      nn   = (int) (DirichletBCs[ii][0]);
      dof  = (int) (DirichletBCs[ii][1]);

      ind = nn*ndof+dof;

      if( midnodeData[nn][0] )
      {
        xx = 0.25*solnTemp(midnodeData[nn][1]*ndof+dof) + 0.25*solnTemp(midnodeData[nn][2]*ndof+dof);

        SolnData.var1applied[ind] = 2.0*(solnTemp(ind) - xx);
      }

      SolnData.var1applied[ind] -= SolnData.var1[ind];
    }

    //printVector(SolnData.var1applied);

    return;
}



void GrowthFEM::setTimeParam()
{
  SolnData.setTimeParam();
  // solve for the initial profile

  return;
}



void GrowthFEM::timeUpdate()
{
    if(DEBUG) {cout << endl; cout << endl; cout <<  " GrowthFEM::timeUpdate() ... STARTED " <<  endl;}

    //update time parameters for the solution algorithm
    SolnData.timeUpdate();

    int  idd = SolnData.ElemProp[0]->id;
    if( (idd == 2010) || (idd == 2060) )
    {
        for(int ee=0; ee<nElem; ee++)
        {
            elems[ee]->presDOF     = elems[ee]->presDOFprev;
            elems[ee]->presDOFprev = elems[ee]->presDOFprev2;
        }
    }

    // copy internal variables intVar1 to intVar2
    if(intVarFlag)
      copyElemInternalVariables();

    // update time functions
    timeFunction[0].update();

    localStiffnessError = 0;
    iterCount = 1;
    filecount++;

    // set Dirichlet boundary conditions
    assignBoundaryConditions();

    if(DEBUG) {cout <<  " GrowthFEM::timeUpdate() ... FINISHED " <<  endl; cout << endl; cout << endl;}

    return;
}



void GrowthFEM::updateIterStep()
{
    SolnData.updateIterStep();

    int  ii, ind, bb;

    for(bb=0;bb<nNode;bb++)
    {
      for(ii=0;ii<ndim;ii++)
      {
        ind = bb*ndof+ii;

        GeomData.NodePosNew[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1[ind];
        GeomData.NodePosCur[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1Cur[ind];
      }
    }

    return;
}


void GrowthFEM::updateGeometry()
{
    //update variables at (n+af)
    SolnData.var1Cur        =  SolnData.var1;
    SolnData.var1DotCur     =  SolnData.var1Dot;
    SolnData.var1DotDotCur  =  SolnData.var1DotDot;

    int  bb, ii, ind;
    //update the coordinates
    for(bb=0;bb<nNode;bb++)
    {
      for(ii=0;ii<ndim;ii++)
      {
        ind = bb*ndof+ii;

        GeomData.NodePosNew[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1[ind];
        GeomData.NodePosCur[bb][ii] = GeomData.NodePosOrig[bb][ii] + SolnData.var1Cur[ind];
      }
    }

    return;
}




void GrowthFEM::saveSolution()
{
    SolnData.saveSolution();

    int  idd = SolnData.ElemProp[0]->id;
    if( (idd == 2010) || (idd == 2060) )
    {
        for(int ee=0; ee<nElem; ee++)
        {
            elems[ee]->presDOFprev2 = elems[ee]->presDOFprev;
            elems[ee]->presDOFprev  = elems[ee]->presDOF;
        }
    }

    if(intVarFlag)
    {
        for(int ee=0; ee<nElem; ee++)
        {
          elems[ee]->ivar.saveSolution();
        }
    }

    return;
}



void GrowthFEM::reset()
{
  SolnData.reset();

  // copy internal variables intVar1 to intVar2
  if(intVarFlag)
    copyElemInternalVariables();

  return;
}


void GrowthFEM::copyElemInternalVariables()
{
  // copy internal variables intVar1 to intVar2

  for(int e=0; e<nElem; e++)
  {
    elems[e]->ivar.reset();
    //elems[e]->intVar = elems[e]->intVarPrev;
  }

  return;
}




void GrowthFEM::writeNodalData()
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

              val  = 0.50*SolnData.var1[n1];
              val += 0.25*SolnData.var1[n2];
              val += 0.25*SolnData.var1[n3];
            }
            else
            {
              val = SolnData.var1[nn*ndof+dof];
            }

        break;

        case  4 : // velocity

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.var1Dot[n1];
              val += 0.25*SolnData.var1Dot[n2];
              val += 0.25*SolnData.var1Dot[n3];
            }
            else
            {
              val = SolnData.var1Dot[nn*ndof+dof];
            }

        break;

        case  5 : // acceleration

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.var1DotDot[n1];
              val += 0.25*SolnData.var1DotDot[n2];
              val += 0.25*SolnData.var1DotDot[n3];
            }
            else
            {
              val = SolnData.var1DotDot[nn*ndof+dof];
            }

        break;

        case  9: // loadfactor for the arc-length method

            val = loadfactor;

        break;

        default : 

              prgError(1,"GrowthFEM::writeNodalData()","invalid value of 'type'!");
        break;
      }

      char        tmp[200];
      MyString    tmpStr;

      sprintf(tmp," \t %12.6E", val);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
    }

    //printData(0, 0);

    return;
}



void GrowthFEM::writeNodalDataExplicitScheme()
{
    int ii, jj, bb, ee, type, nn, n1, n2, n3, dof;
    double val, fact=1.0;

    fout_explicit << myTime.cur << '\t' << timeFunction[0].prop ;

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

              val  = 0.50*SolnData.var1[n1];
              val += 0.25*SolnData.var1[n2];
              val += 0.25*SolnData.var1[n3];
            }
            else
            {
              val = SolnData.var1[nn*ndof+dof];
            }

        break;

        case  4 : // velocity

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.var1Dot[n1];
              val += 0.25*SolnData.var1Dot[n2];
              val += 0.25*SolnData.var1Dot[n3];
            }
            else
            {
              val = SolnData.var1Dot[nn*ndof+dof];
            }

        break;

        case  5 : // acceleration

            if( midnodeData[nn][0] ) // if midnode
            {
              n1 = nn*ndof+dof;
              n2 = midnodeData[nn][1]*ndof+dof;
              n3 = midnodeData[nn][2]*ndof+dof;

              val  = 0.50*SolnData.var1DotDot[n1];
              val += 0.25*SolnData.var1DotDot[n2];
              val += 0.25*SolnData.var1DotDot[n3];
            }
            else
            {
              val = SolnData.var1DotDot[nn*ndof+dof];
            }

        break;

        default : 

               prgError(1,"GrowthFEM::writeNodalData()","invalid value of 'type'!");
        break;
      }

      fout_explicit << '\t' << val ;
    }
    fout_explicit << endl;

    return;
}




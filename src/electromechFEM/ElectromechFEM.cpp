
#include "ElectromechFEM.h"
#include "DataBlockTemplate.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "NewElement.h"
#include "NewMaterial.h"
#include "List.h"
#include <algorithm>


extern List<TimeFunction> timeFunction;
extern MyTime           myTime;

using namespace std;

int ElementBase::elemcount = 0;


ElectromechFEM::ElectromechFEM()
{
    ndof = nElem = 0;

    dispDOF = presDOF = epotDOF = tempDOF = totalDOF = 0;

    nNode = npElem = ndof = 0;

    firstIter = STAGGERED = true;
    IMPLICIT_SOLVER = true;
    MIXED_ELEMENT = false;
    ARCLENGTH = false;

    tol = -2.0;
    lamb = 0.0; lambPrev = 0.0; lambStep = 0.0; lambIter = 0.0;

    solverEigen  = NULL;
    //solverPetsc  = NULL;
    elems  = NULL;
    MatlData = NULL;

    localStiffnessError = 0;
    filecount = 0;
    IterNum = 1;

    scaVTK    =  vtkSmartPointer<vtkFloatArray>::New();
    presVTK   =  vtkSmartPointer<vtkFloatArray>::New();

}


ElectromechFEM::~ElectromechFEM()
{
  if(elems != NULL)
  {
    for(int ii=0;ii<nElem;ii++)
      delete elems[ii];

    delete [] elems;
    elems = NULL;
  }

  if(elemsFaces != NULL)
  {
    for(int ii=0;ii<nElemFaces;ii++)
      delete elemsFaces[ii];

    delete [] elemsFaces;
    elemsFaces = NULL;
  }

  if(solverEigen != NULL)
    delete solverEigen;
  solverEigen = NULL;

  if(MatlData != NULL)
    delete  MatlData;
  MatlData = NULL;
}




void  ElectromechFEM::readInputData(std::ifstream &Ifile, MyString &line)
{
/*
  char fct[] = "ElectromechFEM::readInputData";

  MyString tmpl, *word;

  char tmp[30];

  int nw, i, j, k, n, nn, ii, bb;

  double fact;

  MyStringList   sTmp;
  List<Vector<int> > lviTmp;
  List<Vector<double> > lvdTmp;
  vector<double>  dblTemp;

  DataBlockTemplate t1, t2;

    electromechfem->key.addNew("element type", //0
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

  switch (domain[ELECTROMECHFEM].key.whichBegins(line))
  {
    case  0: //cout << "     ElectromechFEM: reading 'element type' ...\n\n";

            //SolnData.ElemProp.add(new PropertyItem(ELEMENTTYPE));
            //SolnData.ElemProp[SolnData.ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            break;

    case  1: //cout << "     ElectromechFEM: reading 'material type' ...\n\n";

            //SolnData.MatlProp.add(new PropertyItem(MATERIAL));
            //SolnData.MatlProp[SolnData.MatlProp.n-1].readInputData(Ifile,line,"input error in 'material type'!");

            break;

    case  2: //cout << "     ElectromechFEM: reading 'nodes' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1, fct, "invalid specification of 'nodes' !");

            nodePosData.resize(lvdTmp.n);

            for(i=0; i<lvdTmp.n; i++)
            {
              for(j=0;j<(lvdTmp[i].n-1);j++)
                nodePosData[i][j] = lvdTmp[i][j+1];
            }

            break;

    case  3: //cout << "     ElectromechFEM: reading 'elements' ...\n\n";

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

    case  4: //cout << "     ElectromechFEM: reading 'prescribed boundary conditions' ...\n\n";

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

              if(dblTemp[1] < ndim)
              {
                DirichletBCs.push_back(dblTemp);
              }
              else if(dblTemp[1] == ndim)
              {
                DirichletBCs_Epot.push_back(dblTemp);
              }
              else if(dblTemp[1] == (ndim+2) )
              {
                DirichletBCs_Pres.push_back(dblTemp);
              }
              else
              {
                prgError(3, fct, "unknown DOF number in 'prescribed boundary conditions' !");
              }
            }

            break;

    case  5: //cout << "     ElectromechFEM: reading 'nodal forces' ...\n\n";

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

    case  6: //cout << "     ElectromechFEM: reading 'element face loads' ...\n\n";

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

    case  7: //cout << "     ElectromechFEM: reading 'nodal data output' ...\n\n";

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

    case   8: //cout << "     ElectromechFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( ElectromechFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            rhoInfty = lvdTmp[0][2];

            if(tis < 100)
              IMPLICIT_SOLVER = true;
            else
              IMPLICIT_SOLVER = false;

            break;

    case  9: //cout << "     ElectromechFEM: reading 'initial conditions' ...\n\n";

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

    case  10: cout << "     ElectromechFEM: reading 'face load elements' ...\n\n";

             if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
               prgError(1,fct,"invalid input in 'face load elements'!");

            ElemFaceLoadData.resize(lvdTmp.n);

            for(i=0; i<lvdTmp.n; i++)
            {
              if(lvdTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'face load elements' !");

              ElemFaceLoadData[i].resize(lvdTmp[i].n-1);

              for(j=0; j<lvdTmp[i].n-1; j++)
                ElemFaceLoadData[i][j] = lvdTmp[i][1+j]-1;
            }

            cout << " aaaaaaaaaaaaa " << lvdTmp.n << endl;

            break;


    case -1: // go and inherit from DOMAIN

            this->Domain::readInputData(Ifile,line);

            ndim  = ndim;
            ndof = ndf;

            break;
  }
*/
  return;
}




void ElectromechFEM::prepareElemProp()
{
/*
    char fct[] = "ElectromechFEM::prepareElemProp()";

    if(SolnData.ElemProp.n < 1)
      prgError(1,fct,"'element property data ' missing!");

    char *elmTypeNamesBFEM[] = ELECTROMECH_ELEMENT_TYPE_NAMES;

    int ii, idd;

    for(ii=0; ii<SolnData.ElemProp.n; ii++)
    {
      // assign correct id of the element type (as stored in the database) based on the element name
      idd = SolnData.ElemProp[ii].name.which(elmTypeNamesBFEM);

      if(idd < 0)
      {
        prgError(2,fct,"unknown element type name!");
      }

      if(ndim)
      idd += 1001;

      cout << "SolnData.ElemProp[ii].id = " << idd << endl;

      SolnData.ElemProp[ii].id = idd;
    }
*/
    return;
}





void ElectromechFEM::prepareMatlProp()
{
/*
    char fct[] = "ElectromechFEM::prepareMatlProp()";

    if(SolnData.MatlProp.n < 1)
      prgError(1,fct,"'patch material property data ' missing!");

    char *matlTypeNames[] = MATERIAL_TYPE_NAMES_ELECTROMECH;

    intVarFlag = false;

    int  ii, idd;
    for(int ii=0; ii<SolnData.MatlProp.n; ii++)
    {
      // assign correct id of the material type (as stored in the database) based on the material name
      cout <<  SolnData.MatlProp[ii].name <<  endl;

      idd = SolnData.MatlProp[ii].name.which(matlTypeNames);

      cout << SolnData.MatlProp[ii].name << endl;

      cout << " idd = " << idd << endl;

      if(idd < 0)
        prgError(2,fct,"unknown material type name!");

      idd += 1001;

      SolnData.MatlProp[ii].id = idd;

      //if( (idd == 3) || (idd == 4) || (idd == 5) )
        //intVarFlag = true;
    }

    MatlData = NewMaterial(SolnData.MatlProp[elemConn[0][1]].id);

    MatlData->matData.resize( SolnData.MatlProp[elemConn[0][1]].data.n );
    for(int ii = 0; ii<SolnData.MatlProp[elemConn[0][1]].data.n; ii++)
      MatlData->matData[ii] = SolnData.MatlProp[elemConn[0][1]].data[ii];
*/
    return;
}






void ElectromechFEM::prepareInputData()
{
    printf("\n     ElectromechFEM::prepareInputData()  .... STARTED ...\n");

    int ii, jj, kk, ee, aa, bb, cc, nn, ind, count, gp, r;

    assert(ndim > 0 && ndim < 4);

    cout << " ndim  = " << ndim << endl;
    cout << " ndim  = " << ndim << endl;
    cout << " ndof = " << ndof << endl;

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

    if(SolnData.MatlProp.size() > 0)
      prepareMatlProp();

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
      elems[ee] = NewElement(SolnData.ElemProp[elemConn[ee][0]]->id);

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
      elems[ee]->MatlData = MatlData;
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

    SolnData.setTimeIncrementType(tis);
    SolnData.setSpectralRadius(rhoInfty);
    SolnData.initialise(nNode*ndof, 0, 0, 0);

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


    // adjust boundary conditions for mid-nodes for quadratic triangle/tetrahedron element
    idd = SolnData.ElemProp[elemConn[0][0]]->id;
    if( (idd == 1002) || (idd == 1003) || (idd == 1004) ||
        (idd == 1005) || (idd == 1006) || (idd == 1007) ||
        (idd == 1012) || (idd == 1013) ||
        (idd == 1042) || (idd == 1043) || (idd == 1044) ||
        (idd == 1045) || (idd == 1046) || (idd == 1047) ||
        (idd == 1052) || (idd == 1053) )
    {
      processForBernsteinElements();
    }


    MIXED_ELEMENT = false;
    if( (idd == 1006) || (idd == 1007) || (idd == 1012) ||
        (idd == 1046) || (idd == 1047) || (idd == 1052) )
    {
      MIXED_ELEMENT = true;
    }

    MatlData->MIXED_ELEMENT = false;
    if( (idd == 1004) || (idd == 1005) || (idd == 1006) || (idd == 1007) || (idd == 1010) || (idd == 1012) || (idd == 1013) ||
        (idd == 1044) || (idd == 1045) || (idd == 1046) || (idd == 1047) || (idd == 1050) || (idd == 1052) || (idd == 1053) )
    {
      MatlData->MIXED_ELEMENT = true;
    }

    dispDegree = presDegree = epotDegree = tempDegree = -1;

    // set the polynomial degrees of different fields
    if( (idd == 1001) || (idd == 1008) || (idd == 1009) || (idd == 1041) || (idd == 1048) || (idd == 1049) )
    {
      dispDegree =  1;    epotDegree =  1;
    }
    else if( (idd == 1002) || (idd == 1003) || (idd == 1042)  )
    {
      dispDegree =  2;    epotDegree =  2;
    }
    else if( idd == 1004 )
    {
      dispDegree =  2;    epotDegree =  1;    presDegree =  0;
    }
    else if( (idd == 1005) || (idd == 1045) )
    {
      dispDegree =  2;    epotDegree =  2;    presDegree =  0;
    }
    else if( (idd == 1013) || (idd == 1053) )
    {
      dispDegree =  2;    epotDegree =  2;    presDegree =  0;
    }
    else if( idd == 1006 )
    {
      dispDegree =  2;    epotDegree =  1;    presDegree =  1;
    }
    else if( (idd == 1007)  || (idd == 1047) || (idd == 1012) || (idd == 1052) )
    {
      dispDegree =  2;    epotDegree =  2;    presDegree =  1;
    }
    else if( (idd == 1010) || (idd == 1011) || (idd == 1050) || (idd == 1051) )
    {
      dispDegree =  1;    epotDegree =  1;    presDegree =  0;
    }
    else
    {
        cerr <<  "Wrong element type to assess polynomial degrees " <<  endl;
        exit(-1);
    }


    cout << " idd  = " << idd << endl;
    cout << " MIXED_ELEMENT = " << MIXED_ELEMENT << endl;

    cout << " dispDegree  = " << dispDegree << endl;
    cout << " presDegree  = " << presDegree << endl;
    cout << " epotDegree  = " << epotDegree << endl;
    cout << " tempDegree  = " << tempDegree << endl;


    // add face load elements
    nElemFaces = ElemFaceLoadData.size();

    elemsFaces = new ElementBase* [nElemFaces];
    for(ee=0; ee<nElemFaces; ee++)
    {
      //elemsFaces[ee] = NewElement(ElemFaceLoadData[ee][0]);
      elemsFaces[ee] = NewElement(SolnData.ElemProp[ElemFaceLoadData[ee][0]]->id);

      elemsFaces[ee]->elenum   =  ee;

      elemsFaces[ee]->elmType  =  ElemFaceLoadData[ee][0];
      elemsFaces[ee]->matType  =  ElemFaceLoadData[ee][1];
      elemsFaces[ee]->secType  =  ElemFaceLoadData[ee][2];

      vector<int>  vecTemp;
      for(ii=0;ii<ElemFaceLoadData[ee].size()-3;ii++)
        vecTemp.push_back(ElemFaceLoadData[ee][3+ii]);

      elemsFaces[ee]->nodeNums = vecTemp;

      elemsFaces[ee]->SolnData = &(SolnData);
      elemsFaces[ee]->GeomData = &(GeomData);
      //elemsFaces[ee]->MatlData = MatlData;

      elemsFaces[ee]->prepareElemData();
    }

    if(CompareDoubles(SolnData.MatlProp[0]->data[0],0.0))
      SolnData.TRULY_INCOMPRESSIBLE = true;
    else
      SolnData.TRULY_INCOMPRESSIBLE = false;

    printf("     ElectromechFEM::prepareInputData()  .... FINISHED ...\n\n");

    return;
}





void ElectromechFEM::processForBernsteinElements()
{
    int idd, ii, jj, ee, n1, n2, nn, dof;
    double xx, yy, zz, disp[3];
    vector<int>  vecTemp;

    idd = SolnData.ElemProp[elemConn[0][0]]->id;

    cout << " idd = " << idd << endl;

    // 2D elements
    if(ndim == 2)
    {
      // triangular elements
      if( (idd == 1002) || (idd == 1003) || (idd == 1004) ||
          (idd == 1005) || (idd == 1006) || (idd == 1007) )
      {
        // loop over all the elements to
        // a.) find the midnodes, and
        // b.) adjacent nodes to midnodes

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
      if( (idd == 1012) || (idd == 1013) )
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
    else // 3D elements
    {
      // 10-noded tetrahedron element
      if( (idd == 1042) || (idd == 1043) || (idd == 1044) ||
          (idd == 1045) || (idd == 1046) || (idd == 1047) )
      {
        int midnodemaptet10[10][2] = {{0,0}, {0,0}, {0,0}, {0,0},
                                      {0,1}, {1,2}, {0,2}, {0,3}, {1,3}, {2,3}};

        // loop over all the elements to
        // a.) find the midnodes, and
        // b.) adjacent nodes to midnodes
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

      // 27-noded hexahedron element
      if( (idd == 1052) || (idd == 1053) )
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

    return;
}





void  ElectromechFEM::prepareDataForPressure()
{
    if( presDegree <= 0 )
    {
      presDOF = 0;

      return;
    }

    SolnData.var2.resize(max(nElem, nNode));
    SolnData.var2.setZero();

    SolnData.var2applied = SolnData.var2;

    int  ee, ii, jj, ind;
    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;
    vector<bool>  NodeType_pres(nNode,true); // all fixed
    vector<int>   ID_pres(nNode,-1), vecTemp;

    // gather all the pressure nodes (dofs)
    //
    pressure_nodes.clear();
/*
    // linear element for pressure
    if( presDegree == (dispDegree-1) )
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecPres.clear();

        for(ii=0; ii<ind; ii++)
        {
            jj = elems[ee]->nodeNums[ii];

            pressure_nodes.push_back(jj);

            elems[ee]->forAssyVecPres.push_back(jj);
        }
      }

      // free all the corner nodes
      for(ii=0;ii<nNode;ii++)
      {
        if(!midnodeData[ii][0])
        {
          NodeType_pres[ii] = false;
        }
      }
    }
    // element-wise constant for pressure
    else if( presDegree == (dispDegree-2) )
    {
      for(ee=0; ee<nElem; ++ee)
      {
        pressure_nodes.push_back(ee);

        elems[ee]->forAssyVecPres.clear();
        elems[ee]->forAssyVecPres.push_back(ee);
      }

      for(ii=0;ii<nElem;ii++)
      {
        NodeType_pres[ii] = false;
      }
    }
    else
    {
      cerr << " Error in 'ElectromechFEM::prepareCouplingDataDispPres()' " << endl;
      cerr << " presDegree and dispDegree combination not supported yet " << endl; 
      cerr << " dispDegree = " << dispDegree << endl;
      cerr << " presDegree = " << presDegree << endl;
    }
*/

    if(idd == 1007) // Tria6/3
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
         (idd == 1012) // Quad9/4
      || (idd == 1047) // Tet10/4
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
    else if(idd == 1052) // Hexa27/8
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

    return;
}






void  ElectromechFEM::prepareDataForElectricPotential()
{
    SolnData.var3.resize(nNode);
    SolnData.var3.setZero();
    SolnData.var3Cur = SolnData.var3;

    SolnData.var3applied = SolnData.var3;


    int  ee, ii, jj;
    int  idd = SolnData.ElemProp[elemConn[0][0]]->id;

    vector<bool>  NodeType_epot(nNode,false); // all free
    vector<int>   ID_epot(nNode,-1);

    // gather all the electric potential field nodes
    //
    elecpote_nodes.clear();

    if( epotDegree == dispDegree )
    {
      nNode_Epot = nNode;

      elecpote_nodes.resize(nNode);
      for(ii=0; ii<nNode; ii++)
        elecpote_nodes[ii] = ii;

      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecEpot = elems[ee]->nodeNums;
      }

      // specified/fixed electric potential DOF
      for(ii=0;ii<DirichletBCs_Epot.size();ii++)
      {
        NodeType_epot[DirichletBCs_Epot[ii][0]] = true;
      }
    }
    else
    {
      for(ee=0; ee<nElem; ++ee)
      {
        elems[ee]->forAssyVecEpot.clear();

        for(ii=0; ii<=ndim; ii++)
        {
          jj = elems[ee]->nodeNums[ii];

          elecpote_nodes.push_back(jj);

          elems[ee]->forAssyVecEpot.push_back(jj);
        }
      }

      //printVector(elecpote_nodes);

      findUnique(elecpote_nodes);
      nNode_Epot = elecpote_nodes.size();
      //printVector(elecpote_nodes);

      // fix mid nodes
      for(ii=0;ii<nNode;ii++)
      {
        if(midnodeData[ii][0])
        {
          NodeType_epot[ii] = true;
        }
      }

      // specified/fixed electric potential DOF
      for(ii=0;ii<DirichletBCs_Epot.size();ii++)
      {
        if(!midnodeData[DirichletBCs_Epot[ii][0]][0])
        {
          NodeType_epot[DirichletBCs_Epot[ii][0]] = true;
        }
      }
    }

    cout << " nNode_Epot = " << nNode_Epot << endl;
    elecpote_nodes_map_g2l.resize(nNode,-1);
    for(ii=0; ii<nNode_Epot; ii++)
    {
      elecpote_nodes_map_g2l[elecpote_nodes[ii]] = ii;
    }
    //cout << "elecpote_nodes" << endl;
    //printVector(elecpote_nodes);
    //cout << "elecpote_nodes_map_g2l" << endl;
    //printVector(elecpote_nodes_map_g2l);

    //cout << "NodeType_epot" << endl;
    //for(ii=0;ii<nNode;ii++)
      //cout << ii << '\t' << NodeType_epot[ii] << endl;

    epotDOF = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if(!NodeType_epot[ii])
        ID_epot[ii] = epotDOF++;
    }
    cout << " epotDOF = " << epotDOF << endl;
    //cout << "ID_epot" << endl;
    //printVector(ID_epot);

    assyForSolnEpot.resize(epotDOF);
    int count = 0;
    for(ii=0;ii<nNode;ii++)
    {
      if( ID_epot[ii] != -1)
        assyForSolnEpot[count++] = ii;
    }

    //cout << "assyForSolnEpot" << endl;
    //printVector(assyForSolnEpot);

    // reassign the pressure nodes
    for(ee=0; ee<nElem; ee++)
    {
      //printVector(elems[ee]->forAssyVecEpot);
      jj = elems[ee]->forAssyVecEpot.size();

      for(ii=0; ii<jj; ii++)
      {
        elems[ee]->forAssyVecEpot[ii] = ID_epot[elems[ee]->forAssyVecEpot[ii]];
      }
      //printVector(elems[ee]->forAssyVecEpot);
    }

    return;
}






void  ElectromechFEM::prepareDataForTemperature()
{
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
    for(ii=0; ii<nNode_Epot; ii++)
    {
      temperature_nodes_map_g2l[elecpote_nodes[ii]] = ii;
    }
    //cout << "elecpote_nodes" << endl;
    //printVector(elecpote_nodes);
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

    return;
}




bool ElectromechFEM::converged()
{
  char fct[] = "ElectromechFEM::converged";

  if (rNorm < tol && localStiffnessError == 0)
    return true;

  return false;
}




bool ElectromechFEM::diverging(double factor)
{
  if (rNormPrev > -0.1 && (rNorm / rNormPrev) > factor) return true;

  if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;

  return false;
}


void ElectromechFEM::printComputerTime(bool reset, int detailFlg)
{
  cout << "----------------------------------------------------" << endl;

  if (reset)
  {
    ctimFactSolvUpdt = 0.;
    ctimCalcStiffRes = 0.;
  }

  return;
}



void ElectromechFEM::setInitialConditions()
{
    int jj, nn, dof, dd, ind;

    double  xx=0.0, yy=0.0, zz=0.0, fact;
    double  specVal, dt=myTime.dt, mtwodt=-2.0*dt, mthreedt=-3.0*dt, omega=105.0;

    for(nn=0; nn<nNode; nn++)
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

      /*
      SolnData.var1[ind]    = analy.computeValue(0, xx, yy, 0.0);
      SolnData.var1[ind+1]  = analy.computeValue(1, xx, yy, 0.0);
      SolnData.var1[ind+2]  = analy.computeValue(2, xx, yy, 0.0);

      SolnData.var1Prev[ind]   = analy.computeValue(0, xx, yy, -dt);
      SolnData.var1Prev[ind+1] = analy.computeValue(1, xx, yy, -dt);
      SolnData.var1Prev[ind+2] = analy.computeValue(2, xx, yy, -dt);

      SolnData.var1Prev2[ind]   = analy.computeValue(0, xx, yy, mtwodt);
      SolnData.var1Prev2[ind+1] = analy.computeValue(1, xx, yy, mtwodt);
      SolnData.var1Prev2[ind+2] = analy.computeValue(2, xx, yy, mtwodt);

      SolnData.var1Prev3[ind]   = analy.computeValue(0, xx, yy, mthreedt);
      SolnData.var1Prev3[ind+1] = analy.computeValue(1, xx, yy, mthreedt);
      SolnData.var1Prev3[ind+2] = analy.computeValue(2, xx, yy, mthreedt);
      */

      //SolnData.var1Dot[ind+1]    = 500.0; // 2D tensile test

      //SolnData.var1Dot[ind]    = 10.0; // Finite strain beam by Gil
      //SolnData.var1Dot[ind+1]    =  5000.0; // 2D beam problem for the paper
      //SolnData.var1Dot[ind+1]    = 20000.0; // 3D connecting rod
      //SolnData.var1Dot[ind+1]    = 2.0; // 3D long beam
      //SolnData.var1Dot[ind+2]    = 1.0; // 3D long beam
      //SolnData.var1Dot[ind]    = -omega*yy;
      //SolnData.var1Dot[ind+1]  =  omega*xx;

      // square prism torsion only - Gil
      // (vx,vy,vz) = 1500 sin(pi*y/12) (z,0,-x) cm/s  (centimetre per second)
      //fact = 1.5*sin(PI*yy/12.0); // cm/ms (centimetre per millisecond)
      fact = 100.0*sin(PI*yy/12.0);
      //SolnData.var1Dot[ind]    =  fact*zz;
      //SolnData.var1Dot[ind+1]  =  0.0;
      //SolnData.var1Dot[ind+2]  = -fact*xx;
    }


    int ii;
    for(ii=0; ii<DirichletBCs.size(); ii++)
    {
      nn  = (int) (DirichletBCs[ii][0]);
      dof = (int) (DirichletBCs[ii][1]);

      SolnData.var1Dot(nn*ndof+dof) = 0.0;
    }

    //printVector(SolnData.var1Dot);

    // adjust the velocity values for the midnoes
    VectorXd  velTemp(SolnData.var1Dot.rows());
    velTemp = SolnData.var1Dot;

    for(nn=0; nn<nNode; nn++)
    {
      if( midnodeData[nn][0] ) // if midnode
      {
        for(dd=0; dd<ndof; dd++)
        {
          xx = 0.25*velTemp(midnodeData[nn][1]*ndof+dd) + 0.25*velTemp(midnodeData[nn][2]*ndof+dd);

          SolnData.var1Dot[nn*ndof+dd] = 2.0*(velTemp(nn*ndof+dd) - xx);
        }
        //cout << velTemp(nn*ndof) << '\t' << SolnData.var1Dot(nn*ndof) << endl;
      }
    }

    //printVector(SolnData.var1Dot);

    return;
}


void ElectromechFEM::assignBoundaryConditions()
{
    //double loadfactor=1.0;
    double loadfactor = timeFunction[0].prop;
    //loadfactor=0.0;
    //if ( timeFunction[0].prop < 1.0)
      //loadfactor = 0.0;
    //else
      //loadfactor = timeFunction[0].prop - 1.0;

    int idd, aa, ee, ii, jj, n1, n2, nn, dof, ind;

    double  xx=0.0, yy=0.0, zz=0.0, theta, disp[3], specVal;

    VectorXd  solnTemp(nNode*ndof);

    double  tCur = myTime.cur;

    SolnData.var1applied.setZero();
    solnTemp.setZero();

    for(ii=0; ii<DirichletBCs.size(); ii++)
    {
      nn       =  (int) (DirichletBCs[ii][0]);
      dof      =  (int) (DirichletBCs[ii][1]);

      xx = GeomData.NodePosOrig[nn][0];
      yy = GeomData.NodePosOrig[nn][1];

      if(ndim == 3)
        zz = GeomData.NodePosOrig[nn][2];

      //DirichletBCs[ii][2] = analy.computeValue(dof, xx, yy, zz, tCur);

      specVal  =  DirichletBCs[ii][2]*loadfactor;

      solnTemp[nn*ndof+dof] =  specVal;

      SolnData.var1applied[nn*ndof+dof] = specVal;
    }

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

    // electric potential DOFs

    SolnData.var3applied.setZero();
    solnTemp.setZero();

    for(ii=0; ii<DirichletBCs_Epot.size(); ii++)
    {
      nn       =  (int) (DirichletBCs_Epot[ii][0]);

      xx = GeomData.NodePosOrig[nn][0];
      yy = GeomData.NodePosOrig[nn][1];

      if(ndim == 3)
        zz = GeomData.NodePosOrig[nn][2];

      specVal  =  DirichletBCs_Epot[ii][2]*loadfactor;

      solnTemp[nn] =  specVal;

      SolnData.var3applied[nn] = specVal;
    }

    for(ii=0; ii<DirichletBCs_Epot.size(); ii++)
    {
      nn   = (int) (DirichletBCs_Epot[ii][0]);

      if( midnodeData[nn][0] )
      {
        xx = 0.25*solnTemp(midnodeData[nn][1]) + 0.25*solnTemp(midnodeData[nn][2]);

        SolnData.var3applied[nn] = 2.0*(solnTemp(nn) - xx);
      }

      SolnData.var3applied[nn] -= SolnData.var3[nn];
    }

    //printVector(SolnData.var3applied);

    // temperature DOFs

    SolnData.var4applied.setZero();
    solnTemp.setZero();

    for(ii=0; ii<DirichletBCs_Temp.size(); ii++)
    {
      nn       =  (int) (DirichletBCs_Temp[ii][0]);

      xx = GeomData.NodePosOrig[nn][0];
      yy = GeomData.NodePosOrig[nn][1];

      if(ndim == 3)
        zz = GeomData.NodePosOrig[nn][2];

      specVal  =  DirichletBCs_Temp[ii][2]*loadfactor;

      solnTemp[nn] =  specVal;

      SolnData.var4applied[nn] = specVal;
    }

    for(ii=0; ii<DirichletBCs_Temp.size(); ii++)
    {
      nn   = (int) (DirichletBCs_Temp[ii][0]);

      if( midnodeData[nn][0] )
      {
        xx = 0.25*solnTemp(midnodeData[nn][1]) + 0.25*solnTemp(midnodeData[nn][2]);

        SolnData.var4applied[nn] = 2.0*(solnTemp(nn) - xx);
      }

      SolnData.var4applied[nn] -= SolnData.var4[nn];
    }

    //printVector(SolnData.var4applied);

    return;
}



void ElectromechFEM::setTimeParam()
{

  SolnData.setTimeParam();

  return;
}



void ElectromechFEM::timeUpdate()
{
    cout <<  " ElectromechFEM::timeUpdate() ... STARTED " <<  endl;
    firstIter = true;
    localStiffnessError = 0;
    iterCount = 1;
    filecount++;

    //update time parameters for the solution algorithm
    SolnData.timeUpdate();

    int  idd = SolnData.ElemProp[0]->id;
    if( (idd == 1005) || (idd == 1010) || (idd == 1011) || (idd == 1045) || (idd == 1050) || (idd == 1051) )
    {
      for(int ee=0;ee<nElem;ee++)
      {
        elems[ee]->presPrev = elems[ee]->pres;
        elems[ee]->JbarPrev = elems[ee]->Jbar;
      }
    }

    // set Dirichlet boundary conditions
    if(firstIter)
      assignBoundaryConditions();

    //solnIncrStep.setZero();
    lambIncrStep = myTime.dt;
    lamb = 0.0;

    if(intVarFlag)
    {
      for(int e=0; e<nElem; e++)
      {
        elems[e]->ivar.saveSolution();
      }
    }

    updateIterStep();

    cout <<  " ElectromechFEM::timeUpdate() ... FINISHED " <<  endl;

    return;
}



void ElectromechFEM::updateIterStep()
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

        // second-order based formulation
        GeomData.specValNew[bb][ii] = SolnData.var1Dot[ind];
        GeomData.specValCur[bb][ii] = SolnData.var1DotCur[ind];
      }
    }

  return;
}


void ElectromechFEM::updateGeometry()
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

        // second-order based formulation
        GeomData.specValNew[bb][ii] = SolnData.var1Dot[ind];
        GeomData.specValCur[bb][ii] = SolnData.var1DotCur[ind];
      }
    }

    return;
}


void ElectromechFEM::reset()
{
  SolnData.reset();

  int  idd = SolnData.ElemProp[0]->id;
  if( (idd == 1005) || (idd == 1010) || (idd == 1011) || (idd == 1045) || (idd == 1050) || (idd == 1051) )
  {
    for(int ee=0;ee<nElem;ee++)
    {
      elems[ee]->pres = elems[ee]->presPrev;
      elems[ee]->Jbar = elems[ee]->JbarPrev;
    }
  }

  // copy internal variables intVar1 to intVar2
  if(intVarFlag)
    copyElemInternalVariables();

  return;
}


void ElectromechFEM::copyElemInternalVariables()
{
  // copy internal variables intVar1 to intVar2

  for(int e=0; e<nElem; e++)
  {
    elems[e]->ivar.reset();
  }

  return;
}



void ElectromechFEM::writeNodalData()
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

              //val = 0.0;
              //for(ii=0; ii<(OutputData.size()-3); ii++)
                //val += SolnData.force[OutputData[ee][3+ii]*ndof+dof];

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

              prgError(1,"ElectromechFEM::writeNodalData()","invalid value of 'type'!");
        break;
      }

      char        tmp[200];
      MyString    tmpStr;

      sprintf(tmp," \t %12.6E", val);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
    }

    return;
}





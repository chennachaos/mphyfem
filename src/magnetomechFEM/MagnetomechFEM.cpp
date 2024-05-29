
#include "MagnetomechFEM.h"
#include "DataBlockTemplate.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "NewElementMagnetoMech.h"
#include "NewMaterial.h"
#include "List.h"
#include <algorithm>
#include "utilitiesmaterial.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>


extern List<TimeFunction> timeFunction;
extern MyTime           myTime;

using namespace std;


MagnetomechFEM::MagnetomechFEM()
{

}


MagnetomechFEM::~MagnetomechFEM()
{

}



void MagnetomechFEM::printLogo()
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





void  MagnetomechFEM::readInputData(std::ifstream &Ifile, MyString &line)
{
/*
  char fct[] = "MagnetomechFEM::readInputData";

  MyString tmpl, *word;

  char tmp[30];

  int nw, i, j, k, n, nn, ii, bb;

  double fact;

  MyStringList   sTmp;
  List<Vector<int> > lviTmp;
  List<Vector<double> > lvdTmp;
  vector<double>  dblTemp;

  DataBlockTemplate t1, t2;

  
    magnetomechfem->key.addNew("element type", //0
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

  switch (domain[MAGNETOMECHFEM].key.whichBegins(line))
  {
    case  0: //cout << "     MagnetomechFEM: reading 'element type' ...\n\n";

            //SolnData.ElemProp.add(new PropertyItem(ELEMENTTYPE));
            //SolnData.ElemProp[SolnData.ElemProp.n-1].readInputData(Ifile,line,"input error in 'element type'!");

            break;

    case  1: //cout << "     MagnetomechFEM: reading 'material type' ...\n\n";

            //SolnData.MatlProp.add(new PropertyItem(MATERIAL));
            //SolnData.MatlProp[SolnData.MatlProp.n-1].readInputData(Ifile,line,"input error in 'material type'!");

            break;

    case  2: //cout << "     MagnetomechFEM: reading 'nodes' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1, fct, "invalid specification of 'nodes' !");

            nodePosData.resize(lvdTmp.n);

            for(i=0; i<lvdTmp.n; i++)
            {
              for(j=0;j<(lvdTmp[i].n-1);j++)
                nodePosData[i][j] = lvdTmp[i][j+1];
            }

            break;

    case  3: //cout << "     MagnetomechFEM: reading 'elements' ...\n\n";

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

    case  4: //cout << "     MagnetomechFEM: reading 'prescribed boundary conditions' ...\n\n";

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

              if(dblTemp[1] < ndm)
              {
                DirichletBCs.push_back(dblTemp);
              }
              else if(dblTemp[1] == ndm)
              {
                DirichletBCs_Mpot.push_back(dblTemp);
              }
              else if(dblTemp[1] == (ndm+2) )
              {
                DirichletBCs_Pres.push_back(dblTemp);
              }
              else
              {
                prgError(3, fct, "unknown DOF number in 'prescribed boundary conditions' !");
              }
            }

            break;

    case  5: //cout << "     MagnetomechFEM: reading 'nodal forces' ...\n\n";

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

    case  6: //cout << "     MagnetomechFEM: reading 'element face loads' ...\n\n";

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

    case  7: //cout << "     MagnetomechFEM: reading 'nodal data output' ...\n\n";

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

    case   8: //cout << "     MagnetomechFEM: reading 'control parameters' ...\n\n";

            if (!prgReadLnBrkSepListVectorDbl(Ifile,line,lvdTmp))
              prgError(1,fct,"invalid input in 'control parameters'!");

            if( lvdTmp[0].n < 3)
              cerr <<  " Error in (( MagnetomechFEM: reading 'control parameters' )) " << endl;

            tol      = lvdTmp[0][0];
            tis      = (int) lvdTmp[0][1];
            specRad = lvdTmp[0][2];

            if(tis < 100)
              IMPLICIT_SOLVER = true;
            else
              IMPLICIT_SOLVER = false;

            break;

    case  9: //cout << "     MagnetomechFEM: reading 'initial conditions' ...\n\n";

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

    case  10: cout << "     MagnetomechFEM: reading 'face load elements' ...\n\n";

             if (!prgReadLnBrkSepListVectorInt(Ifile,line,lviTmp))
               prgError(1,fct,"invalid input in 'face load elements'!");

            elemConnFaces.resize(lviTmp.n);

            for(i=0; i<lviTmp.n; i++)
            {
              if(lviTmp[i].n < 1)
                prgError(2, fct, "invalid number of 'face load elements' !");

              elemConnFaces[i].resize(lviTmp[i].n-1);

              for(j=0; j<lviTmp[i].n-1; j++)
                elemConnFaces[i][j] = lviTmp[i][1+j]-1;
            }

            cout << " aaaaaaaaaaaaa " << lviTmp.n << endl;

            break;


    case -1: // go and inherit from DOMAIN

            this->Domain::readInputData(Ifile,line);

            ndim  = ndm;
            ndof = ndf;

            break;
  }
*/
  return;
}




void MagnetomechFEM::readControlParameters(string& fname)
{
    cout << " Reading input data \n\n " << endl;

    std::ifstream  infile( fname );
    //std::ifstream  infile( string(infilename+".dat") );

    if(infile.fail())
    {
        cout << " Could not open the input nodes file " << endl;
        exit(1);
    }


    string line2;
    getline(infile,line2);


    MyString  line;

    line = "dimension";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "nodes";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "solid elements";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "prescribed boundary conditions";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "nodal forces";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "nodal data output";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "element type";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "material type";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    line = "control";

    readInputData(infile, line);

    ///////////////////////////////////////////////////
    ///////////////////////////////////////////////////

    return;
}







void MagnetomechFEM::prepareElemProp()
{
    if(SolnData.ElemProp.size() < 1)
      throw  runtime_error("\n\n Element property data missing! \n\n");


    for(int ii=0; ii<SolnData.ElemProp.size(); ii++)
    {
      string matkey(SolnData.ElemProp[ii]->name);

      SolnData.ElemProp[ii]->id = getElementID_MagnetoMech(matkey);

      cout << "SolnData.ElemProp[ii]->id = " << SolnData.ElemProp[ii]->id << endl;
    }

    return;
}





void MagnetomechFEM::prepareMatlProp()
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






void MagnetomechFEM::prepareInputData()
{
    printf("\n     MagnetomechFEM::prepareInputData()  .... STARTED ...\n");

    cout << " ndim  = " << ndim  << endl;
    cout << " nNode = " << nNode << endl;
    cout << " nElem = " << nElem << endl;
    cout << " ndof  = " << ndof  << endl;

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

    int ii, jj, kk, ee, aa, bb, cc, nn, ind, count, gp, r;

    npElem = elemConn[0].size()-3;

    elems = new ElementBase* [nElem];

    for(ee=0;ee<nElem;ee++)
    {
      elems[ee] = NewElementMagnetomechFEM(SolnData.ElemProp[elemConn[ee][0]]->id);

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

    elemsFaces = new ElementBase* [nElemFaces];
    for(ee=0; ee<nElemFaces; ee++)
    {
      //elemsFaces[ee] = NewElement(elemConnFaces[ee][0]);
      cout << "idd = " << SolnData.ElemProp[elemConnFaces[ee][0]]->id <<  endl;

      elemsFaces[ee] = NewElementMagnetomechFEM(SolnData.ElemProp[elemConnFaces[ee][0]]->id);

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

    // adjust boundary conditions for mid-nodes for quadratic triangle/tetrahedron element
    processForBernsteinElements();



    MIXED_ELEMENT = true;
    if( (idd == 2010) || (idd == 2060) )
      MIXED_ELEMENT = false;

    // MIXED_ELEMENT flag in the material object is always true since the formulation is the mixed formulation
    for(ii=0; ii<MatlDataList.size(); ii++)
      MatlDataList[ii]->MIXED_ELEMENT = true;


    COUPLED_PROBLEM = false;
    if( (idd == 2003) || (idd == 2004) || 
        (idd == 2054) || (idd == 2055) || (idd == 2056) )
        COUPLED_PROBLEM = true;


    dispDegree = presDegree = mpotDegree = tempDegree = -1;

    // set the polynomial degrees of different fields

    if( (idd == 2010) || (idd == 2060) )
    {
      dispDegree =  1;    presDegree =  0;
    }
    else if( (idd == 2001)  || (idd == 2002) || (idd == 2051) || (idd == 2052) || (idd == 2053) )
    {
      dispDegree =  2;    presDegree =  1;
    }
    else if( (idd == 2003)  || (idd == 2004) || (idd == 2054) || (idd == 2055) || (idd == 2056) )
    {
      dispDegree =  2;    mpotDegree =  2;    presDegree =  1;
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
    cout << " mpotDegree  = " << mpotDegree << endl;
    cout << " tempDegree  = " << tempDegree << endl;



    SolnData.TRULY_INCOMPRESSIBLE = false;
    for(ii=0; ii<MatlDataList.size(); ii++)
    {
      if(CompareDoubles(MatlDataList[ii]->getKinv(),0.0))
        SolnData.TRULY_INCOMPRESSIBLE = true;
    }

    printf("     MagnetomechFEM::prepareInputData()  .... FINISHED ...\n\n");

    return;
}





void MagnetomechFEM::processForBernsteinElements()
{
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

    return;
}




/*
void  MagnetomechFEM::prepareDataForPressure()
{
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

    return;
}
*/


//
void  MagnetomechFEM::prepareDataForPressure()
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
      cerr << " Wrong element type in 'MagnetomechFEM::prepareDataForPressure()' " << endl;
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




void  MagnetomechFEM::prepareDataForMagneticPotential()
{
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

    return;
}






void  MagnetomechFEM::prepareDataForTemperature()
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

    return;
}




bool MagnetomechFEM::converged()
{
  char fct[] = "MagnetomechFEM::converged";

  if (rNorm < tol && localStiffnessError == 0)
    return true;

  return false;
}




bool MagnetomechFEM::diverging(double factor)
{
  if (rNormPrev > -0.1 && (rNorm / rNormPrev) > factor) return true;

  //if (localStiffnessError != 0) return true;

  if (prgNAN(rNorm)) return true;

  return false;
}


void MagnetomechFEM::printComputerTime(bool reset, int detailFlg)
{
  cout << "----------------------------------------------------\n" << endl;

  //COUT; printf("MagnetomechFEM::calcStiffnessRes:  %7.3f sec ->%5.1f %\n", ctimCalcStiffRes, ctimCalcStiffRes/ctimSinceLastCall*100.);

  //COUT; printf("MagnetomechFEM::factoriseSolvUpdt: %7.3f sec ->%5.1f %\n", ctimFactSolvUpdt, ctimFactSolvUpdt/ctimSinceLastCall*100.);

  if (reset)
  {
    ctimFactSolvUpdt = 0.;
    ctimCalcStiffRes = 0.;
  }

  return;
}



void MagnetomechFEM::setInitialConditions()
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

      //SolnData.var1Dot[ind+2]  = 100.0;
      //SolnData.var1Dot[ind+2]  = 100.0*xx/100.0;
      //SolnData.var1Dot[ind+1]  = zz*10;
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
    SolnData.saveSolution();

    return;
}


void MagnetomechFEM::assignBoundaryConditions()
{
    if(DEBUG) {cout << endl; cout << endl; cout <<  " MagnetomechFEM::assignBoundaryConditions() ... STARTED " <<  endl;}

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

    //printVector(SolnData.var1applied);

    // pressure DOFs

    SolnData.var2applied.setZero();


    // magnetic potential DOFs

    SolnData.var3applied.setZero();
    solnTemp.setZero();

    for(ii=0; ii<DirichletBCs_Mpot.size(); ii++)
    {
      nn       =  (int) (DirichletBCs_Mpot[ii][0]);

      specVal  =  DirichletBCs_Mpot[ii][2]*loadfactor;

      solnTemp[nn] =  specVal;

      SolnData.var3applied[nn] = specVal;
    }

    //cout << "abcd " << endl;
    //printVector(SolnData.var3);
    //printVector(SolnData.var3applied);

    for(ii=0; ii<DirichletBCs_Mpot.size(); ii++)
    {
      nn   = (int) (DirichletBCs_Mpot[ii][0]);

      if( midnodeData[nn][0] )
      {
        xx = 0.25*solnTemp(midnodeData[nn][1]) + 0.25*solnTemp(midnodeData[nn][2]);

        SolnData.var3applied[nn] = 2.0*(solnTemp(nn) - xx);
      }

      //cout << "ii =  " << ii << '\t' << nn << endl;

      SolnData.var3applied[nn] -= SolnData.var3[nn];
    }

    //printVector(SolnData.var3applied);

    if(mpotDOF == 0)
        SolnData.var3 += SolnData.var3applied;

    if(DEBUG) {cout << endl; cout << endl; cout <<  " MagnetomechFEM::assignBoundaryConditions() ... ENDED " <<  endl;}

    return;
}



void MagnetomechFEM::setTimeParam()
{

  SolnData.setTimeParam();

  return;
}



void MagnetomechFEM::timeUpdate()
{
    if(DEBUG) {cout << endl; cout << endl; cout <<  " MagnetomechFEM::timeUpdate() ... STARTED " <<  endl;}

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

    if(DEBUG) {cout <<  " MagnetomechFEM::timeUpdate() ... FINISHED " <<  endl; cout << endl; cout << endl;}

    return;
}



void MagnetomechFEM::updateIterStep()
{
    if(DEBUG) {cout << "     MagnetomechFEM::updateIterStep ...STARTED \n\n";}

    SolnData.updateIterStep();

    //printVector(SolnData.var1DotDot);

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

    if(DEBUG) {cout << "     MagnetomechFEM::updateIterStep ...ENDED \n\n";}

    return;
}


void MagnetomechFEM::updateGeometry()
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



void MagnetomechFEM::saveSolution()
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





void MagnetomechFEM::reset()
{
    SolnData.reset();

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

    return;
}



void MagnetomechFEM::copyElemInternalVariables()
{
    // copy internal variables intVarPrev to intVar

    for(int ee=0; ee<nElem; ee++)
    {
      elems[ee]->ivar.reset();
    }

    return;
}



void MagnetomechFEM::writeNodalData()
{
    cout << " MagnetomechFEM::writeNodalData() " << endl;

    int ii, jj, bb, ee, type, nn, n1, n2, n3, dof;
    double val;

    char        tmp[200];
    MyString    tmpStr;

    sprintf(tmp," \t %12.6f \t %12.6f", myTime.cur, timeFunction[0].prop);

    tmpStr.append(tmp);

    prgWriteToTFile(tmpStr);


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

              prgError(1,"MagnetomechFEM::writeNodalData()","invalid value of 'type'!");
        break;
      }

      char        tmp[200];
      MyString    tmpStr;

      sprintf(tmp," \t %12.6E", val);

      tmpStr.append(tmp);

      prgWriteToTFile(tmpStr);
    }

    sprintf(tmp," \n");
    tmpStr.free();
    tmpStr.append(tmp);

    prgWriteToTFile(tmpStr);


    return;
}






#include "BaseFEM.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "NewMaterial.h"
#include <algorithm>
#include "utilitiesmaterial.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>

extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime               myTime;


using namespace std;



BaseFEM::BaseFEM()
{
    nElem = nNode = npElem = ndof = 0;

    dispDOF = presDOF = mpotDOF = tempDOF = totalDOF = 0;
    numMaterials = 0;

    IMPLICIT_SOLVER = true;
    MIXED_ELEMENT = false;
    ARC_LENGTH = false;
    DEBUG = false;

    conv_tol = -2.0;

    solverEigen  = NULL;
    solverPetsc  = NULL;
    elems  = NULL;
    elemsFaces = NULL;
    //MatlData = NULL;

    localStiffnessError = 0;
    fileCount = 0;
    IterNum = 1;

    loadStepConverged = 0;

    convergedFlagPrev = convergedFlag = false;

    loadFactor = 0.0; loadFactorPrev = 0.0; loadFactorPrev2 = 0.0;
    loadfactorVec.push_back(loadFactor);

    timeFinal = 1.0;
    itersMax  = 10;
    stepsMax  = 1000000;
}



BaseFEM::~BaseFEM()
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

  if(solverPetsc != NULL)
    delete solverPetsc;
  solverPetsc = NULL;

  for(vector<MaterialBase*>::iterator pObj = MatlDataList.begin(); pObj != MatlDataList.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  MatlDataList.clear(); // Purge the contents so no one tries to delete them again
}


void  BaseFEM::deallocatePetscObjects()
{
  if(solverPetsc != NULL)
    solverPetsc->free();

  return;
}



void BaseFEM::readInput(string& fname)
{
    cout << " Entering <BaseFEM::readInput> " << endl;

    ifstream  infile(fname);

    if(infile.fail())
    {
       cout << " Could not open input file " << endl;
       exit(-1);
    }

    vector<string>  stringlist;

    boost::trim(fname);
    //cout << line << endl;
    boost::algorithm::split(stringlist, fname, boost::is_any_of("."), boost::token_compress_on);

    std::string extn1=stringlist[stringlist.size()-1];

    if(extn1 == "msh")
      readInputGMSH(fname);
    else
      readInputMine(fname);

    return;
}



void BaseFEM::readInputGMSH(string& fname)
{
    cout << " BaseFEM::readInputGMSH " << endl;

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

    //map_elem_npelem = {1:2, 2:3, 3:4, 4:4, 5:8, 6:6, 7:5, 8:3, 9:6, 10:9, 11:10, 12:27, 13:18, 14:14, 15:1};

    while(getline(infile,line))
    {
      //std::cout << line << std::endl;

      boost::trim(line);

      // Physical names --- for BCs etc
      //
      if(line == "$PhysicalNames")
      {
        getline(infile,line);
        cout << line << endl;

        boost::trim(line);

        nPhysicalNames = stoi(line);
        cout << nPhysicalNames << endl;

        PhysicalNames.resize(nPhysicalNames);
        PhysicalNamesDims.resize(nPhysicalNames);
        PhysicalNamesTags.resize(nPhysicalNames);

        for(ii=0; ii<nPhysicalNames; ++ii)
        {
          getline(infile,line);
          cout << line << endl;

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
        cout << line << endl;

        ndim = *max_element(PhysicalNamesDims.begin(), PhysicalNamesDims.end());

        for(ii=0; ii<nPhysicalNames; ++ii)
        {
            if(PhysicalNamesDims[ii] == ndim)
              PhysicalNamesTagsDomains.push_back(PhysicalNamesTags[ii]);
            else
              PhysicalNamesTagsBoundaries.push_back(PhysicalNamesTags[ii]);
        }

        cout << " PhysicalNamesTagsDomains ... " << endl;
        printVector(PhysicalNamesTagsDomains);
        cout << endl;
        cout << " PhysicalNamesTagsBoundaries ... " << endl;
        printVector(PhysicalNamesTagsBoundaries);
        cout << endl;

        // prepare boundary patches

        numBoundaryPathes = 0;
        for(ii=0; ii<nPhysicalNames; ++ii)
        {
            if(PhysicalNamesDims[ii] < ndim)
                numBoundaryPathes++;
        }

        //BoundaryPatches.resize(numBoundaryPathes);

        count = 0;
        for(ii=0; ii<nPhysicalNames; ++ii)
        {
            if(PhysicalNamesDims[ii] < ndim)
            {
                BoundaryPatches[count].setLabel(PhysicalNames[ii]);
                BoundaryPatches[count].setTag(PhysicalNamesTags[ii]);
                count++;
            }
        }
      }

      // nodes
      //
      if(line == "$Nodes")
      {
        getline(infile,line);
        //cout << line << endl;
        boost::trim(line);

        nNode = stoi(line);
        cout << " nNode = " << nNode << endl;

        nodePosData.resize(nNode);

        for(ii=0; ii<nNode; ++ii)
        {
          getline(infile,line);
          //cout << line << endl;
          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

          nodePosData[ii][0] = stof(stringlist[1]);
          nodePosData[ii][1] = stof(stringlist[2]);
          nodePosData[ii][2] = stof(stringlist[3]);
        }

        getline(infile,line);
        cout << line << endl;
      }

      if(line == "$Elements")
      {
        getline(infile,line);
        cout << line << endl;
        boost::trim(line);
        boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

        count = stoi(stringlist[0]);
        cout << " count = " << count << endl;

        elemConnTemp.resize(count);
        elemConn.resize(count);

        // elm-number elm-type number-of-tags < tag > … node-number-list
        //
        nElem = 0;
        for(ee=0; ee<count; ++ee)
        {
          getline(infile,line);
          //cout << line << endl;
          //strip_spaces(line, stringlist);

          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

          ind = stringlist.size();
          cout << "ind = " << ind << endl;

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
              elemConn[nElem++] = vecIntTemp;
           }
          else
          {
              cout << "tag = " << tag << endl;
              for(j=0; j<numBoundaryPathes; ++j)
              {
                  if( tag == BoundaryPatches[j].getTag() )
                  {
                      index = j;
                      break;
                  }
              }
              cout << " index = " << index << endl;
              BoundaryPatches[index].addElement(vecIntTemp);
          }
        }
        
        getline(infile,line);
        cout << line << endl;
      }
    }


    //
    // Process data
    //
        for(ii=0; ii<numBoundaryPathes; ++ii)
        {
            BoundaryPatches[ii].processData();
        }

    cout << " nElem = " << nElem << '\t' << elemConn.size() << endl;

    while(elemConn.size() > nElem)
        elemConn.pop_back();

    cout << " nElem = " << nElem << '\t' << elemConn.size() << endl;
    
    for(ii=0; ii<nElem; ++ii)
      printVector(elemConn[ii]);

    cout << endl;
    cout << endl;

    for(ii=0; ii<numBoundaryPathes; ++ii)
    {
        BoundaryPatches[ii].printData();
    }
    


/*
        BoundaryNodes.resize(nPhysicalNames);

        for(ii=0; ii<nPhysicalNames; ++ii)
        {
          getline(infile,line);
          cout << line << endl;

          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

          boost::trim(stringlist[0]);
          PhysicalNamesDims[ii] = stoi(stringlist[0]);

          boost::trim(stringlist[1]);
          PhysicalNamesTags[ii] = stoi(stringlist[1]);

          boost::trim(stringlist[2]);
          PhysicalNames[ii] = stringlist[2].substr(1, stringlist[2].length()-2);
        }

        ndim = std::max_element(PhysicalNamesDims.begin(), PhysicalNamesDims.end());


        elemConn.resize(ind);

        for(ee=0; ee<ind; ++ee)
        {
          getline(infile,line);
          //cout << line << endl;
          //strip_spaces(line, stringlist);

          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);

          count = stringlist.size()-5;
          //cout << "count = " << count << endl;

          elemConnTemp[ee].resize(count);

          for(int j=0; j<count; j++)
            elemConnTemp[ee][j] = stoi(stringlist[j+5]) - 1;
        }

        
        npElem = 0;
        nElem  = 0;
        int elcount=1;

*/
    
    
    return;
}





void BaseFEM::readInputMine(string& fname)
{
    cout << " BaseFEM::readInput " << endl;

    ifstream  infile(fname);

    if(infile.fail())
    {
       cout << " Could not open input file " << endl;
       exit(-1);
    }


    string line;
    vector<string>  stringlist;

    while(getline(infile,line))
    {
        //if (line == "") continue;

        //std::cout << line << std::endl;

        //cout << " size       = " << line.size() << endl;

        boost::trim(line);
 
        //cout << " size       = " << line.size() << endl;
        //cout << " first word = " << line[0] << endl;

        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            //std::cout << line << std::endl;

            if(line.compare(string("Dimension")) == 0)
            {
                cout << "Dimension" << endl;
                readDimension(infile, line);
            }
            else if(line.compare(string("Nodes")) == 0)
            {
                cout << "Nodes" << endl;
                readNodes(infile, line);
            }
            else if(line.compare(string("Elements")) == 0)
            {
                cout << "Elements" << endl;
                readSolidElements(infile, line);
            }
            else if(line.compare(string("Prescribed Boundary Conditions")) == 0)
            {
                cout << "Prescribed Boundary Conditions" << endl;
                this->readPrescribedBCs(infile, line);
            }
            else if(line.compare(string("Face Load Elements")) == 0)
            {
                cout << "Face Load Elements" << endl;
                readSurfaceElements(infile, line);
            }
            else if(line.compare(string("Nodal Forces")) == 0)
            {
                cout << "Nodal Forces" << endl;
                readNodalForces(infile, line);
            }
            else if(line.compare(string("Element Type")) == 0)
            {
                cout << "Element Type" << endl;
                readElementProps(infile, line);
            }
            else if(line.compare(string("Material Type")) == 0)
            {
                cout << "Material Type" << endl;
                readMaterialProps(infile, line);
            }
            else if(line.compare(string("Time Function")) == 0)
            {
                cout << "Time Function" << endl;
                readTimeFunctions(infile, line);
            }
            else if(line.compare(string("Solver")) == 0)
            {
                cout << "Solver" << endl;
                readSolverDetails(infile, line);
            }
            else if(line.compare(string("Nodal Data Output")) == 0)
            {
                cout << "Output" << endl;
                readNodalDataForOutput(infile, line);
            }
            else if(line.compare(string("Output")) == 0)
            {
                cout << "Output" << endl;
                readOutputDetails(infile, line);

                //printVector(outputlist_nodal);
                //printVector(outputlist_elemental);
            }
            else if(line.compare(string("Target Shape")) == 0)
            {
                cout << "Target Shape" << endl;
                readTargetShape(infile, line);
            }
            else
            {
                cout << "key =  " <<  line << endl;
                throw runtime_error("Key not found in BaseFEM::readInput ...");
                //return -1;
            }

            //boost::split(stringlist, line, boost::is_any_of(" "));

            // trim any leading and trailing white spaces in the sub-strings.
            //for(auto& str: stringlist)  boost::trim(str);

            /*
            unordered_map<string,int>   map_keys = {
                        {"Element",       1},
                        {"Material",      2},
                        {"Time Function", 3},
                        {"Solver",        4},
                        {"Output",        5},
                        {"dummy1",        6},
                        {"dummy2",        7}
                        };

            unordered_map<string,int>::const_iterator got = map_keys.infiled(line);

            if( got == map_keys.end() )
            {
                throw runtime_error("Key not found in StandardFEM::readControlParameters ...");
                //return -1;
            }

            cout << " got->first  = " << got->first  << endl;
            cout << " got->second = " << got->second << endl;
            */

      }
    }

    infile.close();

    cout << " Input file is successfully read " << endl;

    return;
}


void BaseFEM::readDimension(ifstream& infile, string& line)
{
    // read the value of dimension
    infile >> ndim ;

    cout << "ndim = " << ndim << endl;

    assert( (ndim > 0) && (ndim < 4) );

    ndof = ndim;

    return;
}



void BaseFEM::readNodes(ifstream& infile, string& line)
{
    // read the number of nodes
    infile >> nNode ;

    cout << "nNode = " << nNode << endl;

    nodePosData.resize(nNode);

    int  ind;
    if(ndim == 1)
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosData[ii][0];
        }
    }
    else if(ndim == 2)
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosData[ii][0] >>  nodePosData[ii][1];
        }
    }
    else
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosData[ii][0] >>  nodePosData[ii][1] >>  nodePosData[ii][2];
        }
    }

    //cout.precision(17);
    //for(int ii=0; ii<nNode; ii++)
      //cout << ii << '\t' << setprecision(12) << nodePosData[ii][0]  << '\t' << setprecision(12) << nodePosData[ii][1]  << '\t' << setprecision(12) << nodePosData[ii][2] << endl;
    //cout << endl;    cout << endl;

    return;
}



void  BaseFEM::readTargetShape(ifstream& infile, string& line)
{
    cout << "nNode = " << nNode << endl;

    nodePosDataTarget.resize(nNode);

    int  ind;
    if(ndim == 1)
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosDataTarget[ii][0];
        }
    }
    else if(ndim == 2)
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosDataTarget[ii][0] >>  nodePosDataTarget[ii][1];
        }
    }
    else
    {
        for(int ii=0; ii<nNode; ii++)
        {
            infile  >>  ind >>  nodePosDataTarget[ii][0] >>  nodePosDataTarget[ii][1] >>  nodePosDataTarget[ii][2];
        }
    }

    return;
}


void  BaseFEM::readSolidElements(ifstream& infile, string& line)
{
    // read the number of elements
    getline(infile,line);    boost::trim(line);

    nElem = stoi(line);

    cout << "nElem = " << nElem << endl;

    vector<string>  stringlist;

    elemConn.resize(nElem);

    for(int ee=0; ee<nElem; ee++)
    {
        getline(infile,line);
        boost::trim(line);
        //cout << line << endl;
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        //stringlist.erase(std::remove(stringlist.begin(), stringlist.end(), " "), stringlist.end());


        int ind = stringlist.size();
        //cout << "ind = " << ind << endl;

        elemConn[ee].resize(ind-1);

        for(int j=1; j<ind; j++)
            elemConn[ee][j-1] = stoi(stringlist[j]) - 1;

        //printVector(elemConn[ee]);
    }

    return;
}



void BaseFEM::readSurfaceElements(ifstream& infile, string& line)
{
    // read the number of elements
    getline(infile,line);    boost::trim(line);

    nElemFaces = stoi(line);

    cout << "nElemFaces = " << nElemFaces << endl;

    vector<string>  stringlist;

    elemConnFaces.resize(nElemFaces);

    for(int ee=0; ee<nElemFaces; ee++)
    {
        getline(infile,line);
        boost::trim(line);
        //cout << line << endl;
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        //stringlist.erase(std::remove(stringlist.begin(), stringlist.end(), " "), stringlist.end());


        int ind = stringlist.size();
        //cout << "ind = " << ind << endl;

        elemConnFaces[ee].resize(ind-1);

        for(int j=1; j<ind; j++)
            elemConnFaces[ee][j-1] = stoi(stringlist[j]) - 1;

        //printVector(elemConnFaces[ee]);
    }

    return;
}



void BaseFEM::readPrescribedBCs(ifstream& infile, string& line)
{
    // read the number of prescribed BCs
    getline(infile,line);    boost::trim(line);

    nDBC = stoi(line);

    cout << "nDBC = " << nDBC << endl;

    vector<string>  stringlist;
    vector<double>  vecTemp(3);

    for(int ee=0; ee<nDBC; ee++)
    {
        getline(infile,line);
        //cout << line << endl;
        boost::trim(line);
        //cout << line << endl;
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        //stringlist.erase(std::remove(stringlist.begin(), stringlist.end(), " "), stringlist.end());


        int ind = stringlist.size();
        //cout << "ind = " << ind << endl;

        int dof = stoi(stringlist[1]);

        //cout << "dof = " << dof << endl;

        vecTemp[0] = stoi(stringlist[0])-1;
        vecTemp[1] = dof-1;
        vecTemp[2] = stod(stringlist[2]);

        if(dof <= ndim)
        {
            DirichletBCs.push_back(vecTemp);
        }
        else if( dof == (ndim+1) )
        {
            DirichletBCs_Pres.push_back(vecTemp);
        }
        else if( dof == (ndim+2) )
        {
            DirichletBCs_Temp.push_back(vecTemp);
        }
        else if( dof == (ndim+3) )
        {
            DirichletBCs_Epot.push_back(vecTemp);
        }
        else if( dof == (ndim+4) )
        {
            DirichletBCs_Mpot.push_back(vecTemp);
        }
        else
        {
            cerr << "unknown DOF number in 'prescribed boundary conditions' !" << endl;
        }

        //printVector(DirichletBCs[ee]);
    }

    return;
}



void BaseFEM::readNodalForces(ifstream& infile, string& line)
{
    // read the number of nodal forces
    getline(infile,line);    boost::trim(line);

    int  nFBCs = stoi(line);

    cout << "nFBCs = " << nFBCs << endl;

    vector<string>  stringlist;
    vector<double>  vecTemp(3);

    for(int ee=0; ee<nFBCs; ee++)
    {
        getline(infile,line);        boost::trim(line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        //stringlist.erase(std::remove(stringlist.begin(), stringlist.end(), " "), stringlist.end());


        int ind = stringlist.size();
        //cout << "ind = " << ind << endl;

        vecTemp[0] = stod(stringlist[0]);
        vecTemp[1] = stod(stringlist[1]);
        vecTemp[2] = stod(stringlist[2]);

        nodeForcesData.push_back(vecTemp);

        //printVector(DirichletBCs[ee]);
    }

    return;
}



void BaseFEM::readElementProps(ifstream& infile, string& line)
{
/*
    nElemTypes = SolnData.ElemProp.size();

    getline(infile,line);    boost::trim(line);

    vector<string>  stringlist;
    boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

    int ind = stringlist.size();
    cout << "ind = " << ind << endl;

    SolnData.ElemProp.push_back(new PropertyItem);

    SolnData.ElemProp[nElemTypes]->setID(stoi(stringlist[0]));

    SolnData.ElemProp[nElemTypes]->setName(stringlist[1]);

    ind -= 2;
    vector<double>  vecTemp(ind);
    for(int ii=0; ii<ind; ii++)
        vecTemp[ii] = stod(stringlist[ii+2]);

    SolnData.ElemProp[nElemTypes]->setData(vecTemp);

    nElemTypes++;
*/
    cout << "nElemTypes = " << nElemTypes << endl;

    return;
}




void BaseFEM::readMaterialProps(ifstream& infile, string& line)
{
    // read {
    getline(infile,line);    boost::trim(line);

    // read id
    getline(infile,line);    boost::trim(line);

    vector<string>  stringlist;
    boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);

    int  id = stoi(stringlist[1]);

    // read name
    getline(infile,line);    boost::trim(line);
    boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);


    std::string matkey(stringlist[1]);

    int  idd = getMaterialID(matkey);

    cout << " idd = " << idd << endl;

    if( (idd == 101) || (idd == 102) || (idd == 103) || (idd == 2006) || (idd == 2007) )
      intVarFlag = true;

    MatlDataList.push_back( NewMaterial(idd) );

    MatlDataList[numMaterials]->setID(id);

    MatlDataList[numMaterials]->SolnData = &(SolnData);

    //MatlDataList[numMaterials]->readData(infile, line);

    numMaterials++;

    cout << "numMaterials = " << numMaterials << endl;

    return;
}


void BaseFEM::readNodalDataForOutput(ifstream& infile, string& line)
{
    // read the number of rows
    getline(infile,line);    boost::trim(line);

    int  nOut = stoi(line);

    vector<string>  stringlist;
    vector<double>  vecTemp(3);

    for(int ee=0; ee<nOut; ee++)
    {
        getline(infile,line);        boost::trim(line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        //stringlist.erase(std::remove(stringlist.begin(), stringlist.end(), " "), stringlist.end());

        int ind = stringlist.size();
        //cout << "ind = " << ind << endl;

        vecTemp[0] = stod(stringlist[0]);
        vecTemp[1] = stod(stringlist[1]);
        vecTemp[2] = stod(stringlist[2]);

        OutputData.push_back(vecTemp);

        //printVector(DirichletBCs[ee]);
    }

    return;
}


void BaseFEM::readOutputDetails(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist[0] == "Frequency")
            {
                outputFrequency = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "Nodal")
            {
                outputlist_nodal.push_back(stringlist[1]);
            }
            else if(stringlist[0] == "Elemental")
            {
                outputlist_elemental.push_back(stringlist[1]);
            }
            else
            {
                throw runtime_error("Option not available in BaseFEM::readOutputDetails");
            }
        }
    }


    return;
}


void BaseFEM::readSolverDetails(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    unordered_map<string,int>   map_keys = {
                        {"library",       1},
                        {"type",          2},
                        {"dt",            3},
                        {"iterations",    4},
                        {"tolerance",     5},
                        {"tis",           6},
                        {"specrad",       7},
                        {"tf",            8},
                        {"maxsteps",      9},
                        {"masslumping",  10},
                        {"debug",        11}
                        };


    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            unordered_map<string,int>::const_iterator got = map_keys.find(stringlist[0]);

            if( got == map_keys.end() )
            {
                throw runtime_error("Key not found in BaseFEM::readSolverDetails ...");
                //return -1;
            }

            cout << " got->first  = " << got->first  << endl;
            cout << " got->second = " << got->second << endl;

            switch(got->second)
            {
                case 1:                                     // library

                    if( (stringlist[1] == "EIGEN") || (stringlist[1] == "eigen") )
                    {
                      SOLVER_LIB_TYPE = SOLVER_LIB_EIGEN;
                    }
                    else if( (stringlist[1] == "PARDISO") || (stringlist[1] == "pardiso") )
                    {
                      SOLVER_LIB_TYPE = SOLVER_LIB_EIGEN_PARDISO;
                    }
                    else if( (stringlist[1] == "PETSC") || (stringlist[1] == "petsc") )
                    {
                      SOLVER_LIB_TYPE = SOLVER_LIB_PETSC;
                    }
                    else
                    {
                        throw runtime_error("Library option not available in BaseFEM::readSolverDetails");
                    }

                break;

                case 2:                                     // algorithm type - Newton-Raphson or Arc-length

                    if(stringlist[1] == "standard")
                    {
                      //NONLINEAR_SOLVER_TYPE = SOLVER_TYPE_NEWTONRAPHSON;
                    }
                    else
                    {
                        throw runtime_error("Solver type not available in BaseFEM::readSolverDetails");
                    }

                break;

                case 3:                                     // time step

                    line = stringlist[1];
                    boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
                    for(auto& str: stringlist)  boost::trim(str);

                    cout << "time step " << stringlist.size() << endl;
                    //printVector(stringlist);

                    if(stringlist.empty())
                    {
                        throw runtime_error("no data for time step in BaseFEM::readSolverDetails");
                    }

                    if(stringlist.size() == 1)
                    {
                      myTime.set( stod(stringlist[0]), stod(stringlist[0]), stod(stringlist[0]) );
                    }
                    else if(stringlist.size() == 2)
                    {
                      myTime.set( stod(stringlist[0]), stod(stringlist[1]), stod(stringlist[0]) );
                    }
                    else
                    {
                      myTime.set( stod(stringlist[0]), stod(stringlist[1]), stod(stringlist[2]) );
                    }

                break;

                case 4:                                     // maximum iterations

                    itersMax = stoi(stringlist[1]);

                    if(itersMax < 0)
                    {
                        throw runtime_error("Maximum number of iterations cannot be negative in BaseFEM::readSolverDetails");
                    }

                break;

                case 5:                                     // convergence tolerance

                    conv_tol = stod(stringlist[1]);

                    if(conv_tol < 0.0)
                    {
                        throw runtime_error("Tolerance cannot be negative in BaseFEM::readSolverDetails");
                    }

                break;

                case 6:                                     // algorithm type - Newton-Raphson or Arc-length

                    //SolnData.setTimeIncrementType(stringlist[1]);

                break;

                case 7:                                     // algorithm type - Newton-Raphson or Arc-length

                    SolnData.setSpectralRadius(stod(stringlist[1]));

                break;

                case 8:                                     // final time

                    timeFinal = stod(stringlist[1]);

                break;

                case 9:                                     // maximum number of steps

                    stepsMax = stoi(stringlist[1]);

                break;

                case 10:                                    // mass lumping

                    SolnData.MassLumping = ( stoi(stringlist[1]) == 1);

                break;

                case 11:                                    // debugging

                    DEBUG = ( stoi(stringlist[1]) == 1);

                break;

                default:

                break;

            }
        }
    }

    return;
}


void BaseFEM::readTimeFunctions(ifstream& infile, string& line)
{
    // read {
    getline(infile,line);    boost::trim(line);

    //prgReadTimeFunctions(infile);

/*
    vector<string>  stringlist;

    while(getline(infile,line))
    {
        boost::trim(line);
        std::cout << line << std::endl;

        if(line[0] == '}') break;
    }
*/
    return;
}






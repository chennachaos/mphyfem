
#include "femBase.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include "NewLagrangeElement.h"
#include "NewMaterial.h"
#include <algorithm>
#include "myIDmaps.h"
#include <boost/algorithm/string.hpp>
#include <unordered_map>

extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern MyTime           myTime;
extern bool  debug;


using namespace std;



femBase::femBase()
{
    if(debug) cout << " femBase::femBase() ... constructor " << endl;

    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    ndof = 0;
    nElem_global = 0;
    nNode_global = 0;

    dispDOF = presDOF = mpotDOF = tempDOF = 0;
    numMaterials = 0;
    numElemTypes = 0;
    numDomains = 0;
    numBoundaryPatches = 0;

    firstIteration = true;
    IMPLICIT_SOLVER = true;
    MIXED_ELEMENT = false;
    ARC_LENGTH = false;
    COUPLED_PROBLEM = false;
    MIXED_ELEMENT_P0 = false;

    conv_tol = -1.0;

    solverEigen  = NULL;
    solverPetsc  = NULL;
    elems  = NULL;
    //MatlData = NULL;

    localStiffnessError = 0;
    filecount = 0;
    IterNum = 1;

    timeFinal = 1.0;
    iterationsMax = 10;
    stepsMax = 1000000;

    loadStepConverged = 0;
    ntotdofs_global = ntotdofs_local = 0;

    convergedFlagPrev = convergedFlag = false;

    loadFactor = loadFactorPrev = loadFactorPrev2 = 0.0;
    loadFactorVec.push_back(loadFactor);

    dispDegree = -1;
    presDegree = -1;
    mpotDegree = -1;

    outputfrequency = 1;
    outputfreq_vtk  = 1;

    dirname = get_current_dir_name();
    vector<string> stringlist;
    boost::algorithm::split(stringlist, dirname, boost::is_any_of("/"), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);

    dirname = stringlist[stringlist.size()-1];

    if(debug) cout << " femBase::femBase() ... constructor " << endl;
}



femBase::~femBase()
{
  if(elems != NULL)
  {
    for(int ii=0;ii<nElem_global;ii++)
      delete elems[ii];

    delete [] elems;
    elems = NULL;
  }

  for(vector<MaterialBase*>::iterator pObj = MatlDataList.begin(); pObj != MatlDataList.end(); ++pObj)
  {
    delete *pObj; // Note that this is deleting what pObj points to, which is a pointer
  }
  MatlDataList.clear(); // Purge the contents so no one tries to delete them again

  fout_nodaldata.close();
}





int  femBase::deallocatePetscObjects()
{
  solverPetsc->free();

  return 0;
}




void femBase::readInputGMSH(string& fname)
{
    cout << " femBase::readInputGMSH " << endl;

    ifstream infile(fname);

    if(infile.fail())
    {
       cout << " Could not open the input mesh file" << endl;
       exit(1);
    }

    string line;
    vector<string>  stringlist;

    int  ii, jj, j, ee, ind, index, count, tag, npElem;
    vector<string>        PhysicalNames;
    vector<int>           vecIntTemp, PhysicalNamesDims, PhysicalNamesTags, PhysicalNamesTagsDomains, PhysicalNamesTagsBoundaries;
    vector<vector<int> >  elemConnTemp;

    vector<int>  nodeNums;

    myPoint  pt;

    //read Mesh format (3 lines)
    getline(infile,line);
    getline(infile,line);
    getline(infile,line);

    //read $PhysicalNames
    getline(infile,line);

    //read number (of physical names)
    getline(infile,line);    boost::trim(line);

    //boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
    //for(auto& str: stringlist)  boost::trim(str);

    numPhysicalNames = stoi(line);
    cout << "numPhysicalNames = " << numPhysicalNames << endl;

    PhysicalNames.resize(numPhysicalNames);
    PhysicalNamesDims.resize(numPhysicalNames);
    PhysicalNamesTags.resize(numPhysicalNames);

    for(ii=0; ii<numPhysicalNames; ++ii)
    {
        getline(infile,line);
        cout << line << endl;

        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        PhysicalNamesDims[ii]  = stoi(stringlist[0]);

        PhysicalNamesTags[ii]  = stoi(stringlist[1]);

        // remove leading and trailing " from the name
        PhysicalNames[ii] = stringlist[2].substr(1,stringlist[2].length()-2);
    }

    //read $EndPhysicalNames
    getline(infile,line);
    cout << line << endl;

    ndim = *max_element(PhysicalNamesDims.begin(), PhysicalNamesDims.end());
    if(ndim == 1) ndim=2;

    for(ii=0; ii<numPhysicalNames; ++ii)
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

    numDomains = std::count(PhysicalNamesDims.begin(), PhysicalNamesDims.end(), ndim);
    numBoundaryPatches = PhysicalNamesDims.size() - numDomains;

    cout << " ndim                  = " << ndim << endl;
    cout << " numDomains            = " << numDomains << endl;
    cout << " numBoundaryPatches    = " << numBoundaryPatches << endl;

    // prepare boundary patches and domains
    numBoundaryPatches = 0;
    numDomains = 0;
    for(ii=0; ii<numPhysicalNames; ++ii)
    {
        if(PhysicalNamesDims[ii] < ndim)
        {
            BoundaryPatches.push_back(make_unique<BoundaryPatch>(PhysicalNames[ii], PhysicalNamesTags[ii], ndim));
            numBoundaryPatches++;
        }
        else
        {
            Domains.push_back(make_unique<Domain>(PhysicalNames[ii], PhysicalNamesTags[ii], ndim));
            numDomains++;
        }
    }
    cout << " numDomains         = " << numDomains << endl;
    cout << " numBoundaryPatches = " << numBoundaryPatches << endl;



    //read $Nodes
    getline(infile,line);
    cout << line << endl;

    // read number (of nodes)
    getline(infile,line);    boost::trim(line);
    cout << line << endl;

    nNode_global = stoi(line);
    cout << " nNode_global = " << nNode_global << endl;


    vector<myPoint>  nodecoords(nNode_global);

    for(ii=0; ii<nNode_global; ++ii)
    {
        getline(infile,line);
        cout << line << endl;

        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        pt[0] = stod(stringlist[1]);
        pt[1] = stod(stringlist[2]);
        pt[2] = stod(stringlist[3]);

        nodecoords[ii] = pt;
    }

    // read $EndNodes
    getline(infile,line);
    cout << line << endl;

    GeomData.setDimension(ndim);
    GeomData.setNodalPositions(nodecoords);


    // read $Elements
    getline(infile,line);
    cout << line << endl;

    // read number (of elements)
    getline(infile,line);    boost::trim(line);
    cout << line << endl;

    count = stoi(line);
    cout << " count = " << count << endl;

    elemConnTemp.resize(count);
    elemConn.resize(count);


    // elm-number elm-type number-of-tags < tag > … node-number-list
    nElem_global = 0;
    for(ee=0; ee<count; ee++)
    {
        getline(infile,line);
        cout << line << endl;

        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        ind = stringlist.size();

        vecIntTemp.resize(ind);
        for(ii=0; ii<ind; ii++)
          vecIntTemp[ii] = stoi(stringlist[ii]);

        elemConnTemp[ee] = vecIntTemp;

        //elm-number elm-type number-of-tags < tag > … node-number-list
        tag = elemConnTemp[ee][3];
        cout << " tag = " << tag << '\t' << ndim << endl;

        if( std::count(PhysicalNamesTagsDomains.begin(), PhysicalNamesTagsDomains.end(), tag) )
        {
            elemConn[nElem_global++] = vecIntTemp;
        }
        else
        {
            npElem = ind-5;
            nodeNums.resize(npElem);
            for(jj=0; jj<npElem; jj++)
              nodeNums[jj] = vecIntTemp[5+jj]-1;

            for(j=0; j<numBoundaryPatches; ++j)
            {
                if( tag == BoundaryPatches[j]->getTag() )
                {
                    index = j;
                    break;
                }
            }
            BoundaryPatches[index]->addElement(nodeNums);
        }
    }

    // read $EndElements
    getline(infile,line);
    cout << line << endl;

    infile.close();


    //
    // Process data
    //
    for(auto& bpt : BoundaryPatches)
    {
        bpt->processData();
        bpt->printData();
    }

    //elemConn.erase(elemConn.begin()+nElem_global, elemConn.end());    // for deletion in range
    while(elemConn.size() > nElem_global)
        elemConn.pop_back();


    string  fnameout = inputfilename+"-nodaloutput.dat";
    fout_nodaldata.open(fnameout, std::ofstream::out | std::ofstream::trunc);

    if(fout_nodaldata.fail())
    {
       cout << " Could not open the Output file for nodal data" << endl;
       exit(1);
    }

    //fout.setf(ios::fixed);
    //fout.setf(ios::showpoint);
    //fout.precision(8);

    return;
}





void femBase::readInput(string& configfname)
{
    cout << " femBase::readInput " << endl;

    //string  configfname = "./inputs/config";
    ifstream  infile(configfname);

    if(infile.fail())
    {
       cout << " Could not open input file " << endl;
       exit(-1);
    }


    string line;
    vector<string>  stringlist;

    while(getline(infile,line))
    {
        cout << line << endl;

        boost::trim(line);
 
        if( isActiveLine(line) )
        {
            //std::cout << line << std::endl;

            if(line.compare(string("Files")) == 0)
            {
                cout << "files" << endl;
                readFiles(infile, line);
            }
            else if(line.compare(string("Model Type")) == 0)
            {
                cout << "Model Type" << endl;
                readModelType(infile, line);
            }
            else if(line.compare(string("Domains")) == 0)
            {
                cout << "Domains" << endl;
                readDomainsData(infile, line);
            }
            else if(line.compare(string("Element")) == 0)
            {
                cout << "Element" << endl;
                readElementProps(infile, line);
            }
            else if(line.compare(string("Material")) == 0)
            {
                cout << "Material" << endl;
                readMaterialProps(infile, line);
            }
            else if(line.compare(string("Body Force")) == 0)
            {
                cout << "Body Force" << endl;
                readBodyForce(infile, line);
            }
            else if(line.compare(string("Time Functions")) == 0)
            {
                cout << "Time Functions" << endl;
                readTimeFunctions(infile, line);
            }
            else if(line.compare(string("Solver")) == 0)
            {
                cout << "Solver" << endl;
                readSolverDetails(infile, line);
            }
            else if(line.compare(string("Boundary Conditions")) == 0)
            {
                cout << "Boundary Conditions" << endl;
                this->readPrescribedBCs(infile, line);
            }
            else if(line.compare(string("Tractions")) == 0)
            {
                cout << "Tractions" << endl;
                readTractions(infile, line);
            }
            else if(line.compare(string("Nodal Forces")) == 0)
            {
                cout << "Nodal Forces" << endl;
                readNodalForces(infile, line);
            }
            else if(line.compare(string("Output")) == 0)
            {
                cout << "Output" << endl;
                readOutputDetails(infile, line);
            }
            else if(line.compare(string("Nodal Data Output")) == 0)
            {
                cout << "Nodal Output" << endl;
                readNodalDataOutputDetails(infile, line);
            }
            else
            {
                cout << "key =  " <<  line << endl;
                throw runtime_error("Key not found in femBase::readInput ...");
                //return -1;
            }
      }
    }

    infile.close();

    cout << " Control parameters are successfully read " << endl;

    return;
}



void femBase::readFiles(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);


    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);
        cout << line << endl;

        if(line[0] == '}') break;

        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            cout << stringlist[0] << '\t' << stringlist[1] << endl;

            if(stringlist[0] == "mesh")
            {
                inputfilename = stringlist[1];
                string gmshfname =  "./inputs/"+inputfilename+".msh";
                readInputGMSH(gmshfname);
            }
            else if(stringlist[0] == "solution")
            {
                restartfilename = stringlist[1];
            }
            else
            {
                cerr << "femBase::readFiles ... Invalid file type" << endl;
            }
        }
    }

    return;
}




void femBase::readModelType(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);


    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);
        cout << line << endl;

        if(line[0] == '}') break;

        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            cout << stringlist[0] << '\t' << stringlist[1] << endl;

            if(stringlist[0] == "dimension")
            {
                ndim = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "behaviour")
            {
                //readSolution(stringlist[1]);
            }
            else if(stringlist[0] == "thickness")
            {
                //readSolution(stringlist[1]);
            }
            else
            {
                cerr << "femBase::readModelType ... Invalid file type" << endl;
            }
        }
    }

    return;
}



void femBase::readDomainsData(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);


    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);
        cout << line << endl;

        if(line[0] == '}') break;

        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            cout << stringlist[0] << endl;

            int index=-1;
            for(auto& dom : Domains)
            {
              cout << " domain label = " << dom->getLabel() << endl;
              if( stringlist[0].compare( dom->getLabel() ) == 0 )
              {
                dom->readData(infile, line);
                index++;
              }
            }

            cout << " index= " << index << endl;
            if(index == -1)
            {
                throw runtime_error("femBase::readDomainsData ... invalid domain label");
            }
        }
    }

    return;
}



void femBase::readBodyForce(ifstream& infile, string& line)
{
    vector<string>  stringlist;
    myPoint  bforce;

    // read {
    getline(infile,line);    boost::trim(line);


    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);
        cout << line << endl;

        if(line[0] == '}') break;

        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            cout << stringlist[0] << '\t' << stringlist[1] << endl;

            if(stringlist[0] == "value")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
                for(auto& str: stringlist)  boost::trim(str);

                bforce[0] = stod(stringlist[0]);
                bforce[1] = stod(stringlist[1]);
                bforce[2] = stod(stringlist[2]);

                GeomData.setBodyForce(bforce);
            }
            else if( (stringlist[0] == "timefunction") || (stringlist[0] == "TimeFunction") )
            {
                GeomData.setBodyForceTimeFunction( stoi(stringlist[1])-1 );
            }
            else
            {
                throw runtime_error("femBase::readBodyForce ... Option not available");
            }
        }
    }

    return;
}





void femBase::readPrescribedBCs(ifstream& infile, string& line)
{
    if(debug) cout << " femBase::readPrescribedBCs " << endl;

    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}')
        {
            getline(infile,line);
            break;
        }


        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            cout << stringlist[0] << endl;

            numBoundaryConditions = boundaryConitions.size();

            boundaryConitions.push_back(make_unique<BoundaryCondition>());

            boundaryConitions[numBoundaryConditions]->label = stringlist[0];

            for(auto& bpt : BoundaryPatches)
            {
              if( stringlist[0] == bpt->getLabel() )
              {
                boundaryConitions[numBoundaryConditions]->bpt = bpt.get();
                break;
              }
            }

            // read {
            getline(infile,line);    boost::trim(line);

            while( infile && (line != "}") )
            {
                getline(infile,line);    boost::trim(line);

                cout << " line " << endl;
                cout << line << endl;

                if(line[0] == '}')
                {
                    getline(infile,line);
                    break;
                }

                if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
                {
                    boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
                    for(auto& str: stringlist)  boost::trim(str);

                    cout << "Inside patch ... " << stringlist[0] << '\t' << stringlist[1] << endl;

                    if(stringlist[0] == "type")
                    {
                        boundaryConitions[numBoundaryConditions]->BCType = stringlist[1];
                    }
                    else if(stringlist[0] == "dof")
                    {
                        boundaryConitions[numBoundaryConditions]->dof_specified_string = stringlist[1];

                        boundaryConitions[numBoundaryConditions]->dof_specified_int = getDOFfromString(stringlist[1]);

                        //set the default time function
                        boundaryConitions[numBoundaryConditions]->timeFunctionNum   = -1;
                    }
                    else if(stringlist[0] == "value")
                    {
                        boundaryConitions[numBoundaryConditions]->expression  = stringlist[1];
                    }
                    else if(stringlist[0] == "timefunction")
                    {
                        boundaryConitions[numBoundaryConditions]->timeFunctionNum   = stoi(stringlist[1])-1;
                    }
                    else
                    {
                        throw runtime_error("Option not available in femSolidmechanics::readBoundaryConditions");
                    }
                }//if
            } //while
        }//if
    }//while

    return;
}





void femBase::readTractions(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}')
        {
            getline(infile,line);
            break;
        }


        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            cout << stringlist[0] << endl;

            numTractionConditions = tractionConitions.size();

            tractionConitions.push_back(make_unique<TractionCondition>());

            tractionConitions[numTractionConditions]->label = stringlist[0];

            for(auto& bpt : BoundaryPatches)
            {
              if( stringlist[0] == bpt->getLabel() )
              {
                tractionConitions[numTractionConditions]->bpt = bpt.get();
                break;
              }
            }

            // read {
            getline(infile,line);    boost::trim(line);

            while( infile && (line != "}") )
            {
                getline(infile,line);    boost::trim(line);

                cout << " line " << endl;
                cout << line << endl;

                if(line[0] == '}')
                {
                    getline(infile,line);
                    break;
                }

                if( isActiveLine(line) )
                {
                    boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
                    for(auto& str: stringlist)  boost::trim(str);

                    cout << "Inside patch ... " << stringlist[0] << '\t' << stringlist[1] << endl;

                    if(stringlist[0] == "type")
                    {
                        tractionConitions[numTractionConditions]->BCType = stringlist[1];
                    }
                    else if(stringlist[0] == "value")
                    {
                        line =  stringlist[1];

                        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
                        for(auto& str: stringlist)  boost::trim(str);

                        assert(stringlist.size() >= ndim);

                        for(int ii=0; ii<ndim; ii++)
                          tractionConitions[numTractionConditions]->specValues[ii] = stod(stringlist[ii]);
                    }
                    else if(stringlist[0] == "timefunction")
                    {
                        tractionConitions[numTractionConditions]->timeFunctionNum   = stoi(stringlist[1])-1;
                    }
                    else
                    {
                        throw runtime_error("Option not available in femSolidmechanics::readTractionConditions");
                    }
                }//if
            } //while
        }//if
    }//while

    return;
}






void femBase::readNodalForces(ifstream& infile, string& line)
{
    vector<string>  stringlist;
    vector<double>  vecTemp(3);

    // read {
    getline(infile,line);    boost::trim(line);

    // read timefunction
    getline(infile,line);    boost::trim(line);
    boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);

    nodalForces.timeFunctionNum = stoi(stringlist[1])-1;

    // read count
    getline(infile,line);    boost::trim(line);
    boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);

    nodalForces.count = stoi(stringlist[1]);

    // read "data"
    for(int ee=0; ee<nodalForces.count; ee++)
    {
        getline(infile,line);        boost::trim(line);
        boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
        for(auto& str: stringlist)  boost::trim(str);

        assert(stringlist.size() == 3);
        //cout << "ind = " << ind << endl;

        vecTemp[0] = stod(stringlist[0])-1;
        vecTemp[1] = stod(stringlist[1])-1;
        vecTemp[2] = stod(stringlist[2]);

        nodalForces.data.push_back(vecTemp);
    }

    // read }
    getline(infile,line);    boost::trim(line);

    return;
}



void femBase::readElementProps(ifstream& infile, string& line)
{
    cout <<  " femBase::readElementProps ... " << endl;

    ElementTypeDataList.push_back( new ElementTypeData );

    ElementTypeDataList[numElemTypes]->readData(infile, line);

    numElemTypes++;

    cout << "numElemTypes = " << numElemTypes << endl;

    cout <<  " femBase::readElementProps ... " << endl;

    return;
}





void femBase::readMaterialProps(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    // read id
    getline(infile,line);    boost::trim(line);

    boost::algorithm::split(stringlist, line, boost::is_any_of(": "), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);

    int  id = stoi(stringlist[1]);

    // read name
    getline(infile,line);    boost::trim(line);
    boost::algorithm::split(stringlist, line, boost::is_any_of(": "), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);


    std::string matkey(stringlist[1]);

    int  matkeyid = getMaterialID(matkey);

    cout << " matkeyid = " << matkeyid << endl;

    vector<int>  matkeywithIVlist = {101, 102, 103, 2006, 2007};

    if( count(matkeywithIVlist.begin(), matkeywithIVlist.end(), matkeyid) )
      intVarFlag = true;

    MatlDataList.push_back( NewMaterial(matkeyid) );

    MatlDataList[numMaterials]->setID(id);

    MatlDataList[numMaterials]->SolnData = &(SolnData);

    MatlDataList[numMaterials]->readInput(infile, line);

    numMaterials++;

    cout << "numMaterials = " << numMaterials << endl;

    return;
}




void femBase::readSolverDetails(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist.size() == 1)
            {
                cerr << "\n\n\n ERROR: Input error in Solver block \n\n\n" << endl;
                exit(-2);
            }

            //cout << stringlist[0] << '\t' << stringlist[1] << endl;

            if(stringlist[0] == "library")
            {
                if( (stringlist[1] == "EIGEN") || (stringlist[1] == "eigen") )
                {
                    SOLVER_LIB_TYPE = SOLVER_LIB_EIGEN;
                }
                else if( (stringlist[1] == "PARDISO") || (stringlist[1] == "pardiso") )
                {
                    SOLVER_LIB_TYPE = SOLVER_LIB_EIGEN_PARDISO;
                }
                else
                {
                    throw runtime_error("Library option not available in femBase::readSolverDetails");
                }
            }
            else if(stringlist[0] == "solvertype") // algorithm type - Newton-Raphson or Arc-length
            {
                if(stringlist[1] == "newton")
                {
                    NONLINEAR_SOLVER_TYPE = SOLVER_TYPE_NEWTONRAPHSON;
                }
                else if(stringlist[1] == "arclength")
                {
                    NONLINEAR_SOLVER_TYPE = SOLVER_TYPE_ARCLENGTH;
                }
                else
                {
                    throw runtime_error("Solver type not available in femBase::readSolverDetails");
                }
            }
            else if(stringlist[0] == "timeStep")  // time step
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
                for(auto& str: stringlist)  boost::trim(str);

                cout << "time step " << stringlist.size() << endl;

                if(stringlist.empty())
                {
                    throw runtime_error("femBase::readSolverDetails ... no data for time step");
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
            }
            else if(stringlist[0] == "maximumIterations")
            {
                iterationsMax = stoi(stringlist[1]);

                if(iterationsMax < 0)
                {
                    throw runtime_error("femBase::readSolverDetails ... Maximum number of iterations cannot be negative");
                }
            }
            else if(stringlist[0] == "tolerance")
            {
                conv_tol = stod(stringlist[1]);

                if(conv_tol < 0.0)
                {
                    throw runtime_error("Tolerance cannot be negative in femBase::readSolverDetails");
                }
            }
            else if(stringlist[0] == "timescheme")
            {
                SolnData.setTimeIncrementScheme( stringlist[1] );
            }
            else if(stringlist[0] == "spectralRadius")
            {
                SolnData.setSpectralRadius(stod(stringlist[1]));
            }
            else if(stringlist[0] == "finalTime")
            {
                timeFinal = stod(stringlist[1]);
            }
            else if(stringlist[0] == "maximumSteps")
            {
                stepsMax = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "outputFrequency")
            {
                outputfreq_vtk = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "masslumping")
            {
                SolnData.MassLumping = ( stoi(stringlist[1]) == 1);
            }
            else if(stringlist[0] == "debug")
            {
                cout << " debug = " << stringlist[0] << endl;
                debug = ( stoi(stringlist[1]) == 1);
            }
        }
    }

    return;
}


void femBase::readTimeFunctions(ifstream& infile, string& line)
{
    vector<string>  stringlist;
    vector<double>  doublelist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            int id = stoi(stringlist[0]);

            if(id == 0)
            {
                cout << " Time function id '0' is not addmissible. It should be 1 or higher. " << endl;
                throw runtime_error(" Error in femSolids::readTimeFunctions...  ");
            }

            if(stringlist.size() < 11)
            {
                cout << " Number of temrs in Time function should be at least 11 " << endl;
                throw runtime_error(" Error in femSolids::readTimeFunctions...  ");
            }

            doublelist.resize(stringlist.size()-1);
            for(int i=0; i<stringlist.size()-1; ++i)
            {
                doublelist[i] = stod(stringlist[i+1]);
            }

            id -= 1;
            if(id == timeFunctions.size())
            {
                timeFunctions.push_back(make_unique<TimeFunction>());
                timeFunctions[id]->setID(id);
            }
            timeFunctions[id]->addTimeBlock(doublelist);
        }
    }

    for(auto& tmf : timeFunctions)
    {
        tmf->printSelf();
        tmf->update();
    }

    return;
}






void femBase::readOutputDetails(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist[0] == "Frequency")
            {
                outputfrequency = stoi(stringlist[1]);
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
                throw runtime_error("femBase::readOutputDetails ... Option not available");
            }
        }
    }

    return;
}



void  femBase::readNodalDataOutputDetails(ifstream& infile, string& line)
{
    vector<string>  stringlist;
    vector<double>  dbllist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( isActiveLine(line) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

            dbllist.clear();
            for(auto& str: stringlist)
            {
                boost::trim(str);
                dbllist.push_back(stod(str));
            }

            NodalDataOutput.push_back(dbllist);
        }
    }

    return;
}




void femBase::readInitialConditions(ifstream& infile, string& line)
{
    vector<string>  stringlist;
    string label;

    // read {
    getline(infile,line);    boost::trim(line);

    initialConitions.resize(ndof);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            label = stringlist[0];

            int dof = getDOFfromString(stringlist[0]);

            if( dof >= ndof)
            {
                throw runtime_error("femBase::readInitialConditions ... Specified DOF higher than #DOF ...");
            }

            initialConitions[dof] = stringlist[1];
        }
    }

    return;
}




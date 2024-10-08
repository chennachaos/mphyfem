
#include <algorithm>
#include <chrono>
#include "femINSmixed.h"
#include "KimMoinFlow.h"
#include "elementutilitiescfd.h"
#include "BernsteinElem2DINSTria6Node.h"
#include "BernsteinElem2DINSTriaP2bP1dc.h"
#include "BernsteinElem2DINSQuad9Node.h"
#include "BernsteinElem3DINSTetra10Node.h"
#include "FunctionsSolver.h"
#include "util.h"
#include "ImmersedIntegrationElement.h"
#include "ImmersedSolid.h"
#include "mymapkeys.h"


#include <boost/algorithm/string.hpp>
#include <unordered_map>
#include <unistd.h>

#include "MyTime.h"
#include "TimeFunction.h"


extern vector<unique_ptr<TimeFunction> > timeFunctions;
extern  MyTime                 myTime;
extern  bool  debug;


using namespace std;


int ElementBase::elemcount = 0;


femINSmixed::femINSmixed()
{
    ndof = 0; nElem_global = 0; nNode_global = 0; fileCount = 0;
    totalDOF_Velo = 0; totalDOF_Pres = 0; totalDOF_Lambda = 0; totalDOF_Solid=0;

    dispDOF  = 0;  presDOF = 0;  lambdaDOF = 0;  totalDOF = 0;

    nImmersedElems = 0;

    nDBC_Velo = 0;
    nDBC_Pres = 0;

    outputFreq = 1;
    
    GRID_CHANGED = false; IB_MOVED = false;
    debug = true;

    SCHEME_TYPE = SCHEME_TYPE_IMPLICIT;

    timeIntegrationScheme = "STEADY";
    spectralRadius = 0.0;
    conv_tol = 1.0e-6;

    td.resize(100);

    AlgoType = 2;

    elems = nullptr;


    dirname = get_current_dir_name();
    vector<string> stringlist;
    boost::algorithm::split(stringlist, dirname, boost::is_any_of("/"), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);

    dirname = stringlist[stringlist.size()-1];
}





femINSmixed::~femINSmixed()
{
    if(elems != NULL)
    {
      for(int ii=0;ii<nElem_global;++ii)
        delete elems[ii];

      delete [] elems;
      elems = NULL;
    }

    phase = -1; error = 0;

    pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &totalDOF, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

}






int  femINSmixed::readInputGMSH(string& fname)
{
    cout << " femINSmixed::readInputGMSH " << endl;

    ifstream  infile(fname);

    if(infile.fail())
    {
       cout << " Could not open input file " << endl;
       exit(-1);
    }

    string line;
    vector<string>  stringlist;

    int  nPhysicalNames, ii, jj, j, ee, ind, index, count, tag, npElem;
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

        ndim = *max_element(PhysicalNamesDims.begin(), PhysicalNamesDims.end());

        for(ii=0; ii<nPhysicalNames; ++ii)
        {
            if(PhysicalNamesDims[ii] == ndim)
              PhysicalNamesTagsDomains.push_back(PhysicalNamesTags[ii]);
            else
              PhysicalNamesTagsBoundaries.push_back(PhysicalNamesTags[ii]);
        }

        //cout << " PhysicalNamesTagsDomains ... " << endl;
        //printVector(PhysicalNamesTagsDomains);
        //cout << endl;
        //cout << " PhysicalNamesTagsBoundaries ... " << endl;
        //printVector(PhysicalNamesTagsBoundaries);
        //cout << endl;


        // prepare boundary patches and domains
        numBoundaryPatches = 0;
        numDomains = 0;
        for(ii=0; ii<nPhysicalNames; ++ii)
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
        cout << " numBoundaryPatches = " << numBoundaryPatches << endl;
        cout << " numDomains         = " << numDomains << endl;
      }

      if(ndim == 1)  ndim = 2;
      ndof = ndim;

      // nodes
      //
      if(line == "$Nodes")
      {
        getline(infile,line);
        boost::trim(line);

        nNode_global = stoi(line);
        cout << " nNode_global = " << nNode_global << endl;

        nodeCoords.resize(nNode_global);

        for(ii=0; ii<nNode_global; ++ii)
        {
          getline(infile,line);
          //cout << line << endl;
          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

          nodeCoords[ii][0] = stod(stringlist[1]);
          nodeCoords[ii][1] = stod(stringlist[2]);
          nodeCoords[ii][2] = stod(stringlist[3]);
        }

        getline(infile,line);
      }

      if(line == "$Elements")
      {
        getline(infile,line);
        cout << line << endl;
        boost::trim(line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

        count = stoi(stringlist[0]);
        cout << " count = " << count << endl;

        elemConnTemp.resize(count);
        elemConn.resize(count);

        // elm-number elm-type number-of-tags < tag > … node-number-list
        //
        nElem_global = 0;
        for(ee=0; ee<count; ++ee)
        {
          getline(infile,line);
          //cout << line << endl;
          //strip_spaces(line, stringlist);

          boost::trim(line);
          boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);

          ind = stringlist.size();
          //cout << "ind = " << ind << endl;

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
              elemConn[nElem_global++] = elemConnTemp[ee];
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
        cout << line << endl;
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

    return 0;
}





int femINSmixed::readConfiguration(string& fname)
{
    cout << " femINSmixed::readConfiguration " << endl;

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
        boost::trim(line);
 
        //cout << " size       = " << line.size() << endl;
        //cout << " first word = " << line[0] << endl;

        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            //std::cout << line << std::endl;

            if(line.compare(string("Domains")) == 0)
            {
                cout << "Domains" << endl;
                // read the Domains
                readDomains(infile, line);
            }
            else if(line.compare(string("Fluid Properties")) == 0)
            {
                cout << "Fluid Properties" << endl;
                readFluidProperties(infile, line);
                //readFluidProps(infile, line);
            }
            else if(line.compare(string("Body Force")) == 0)
            {
                cout << "Body Force" << endl;
                readBodyForce(infile, line);
            }
            else if(line.compare(string("Boundary Conditions")) == 0)
            {
                cout << "Boundary Conditions" << endl;
                // read the number of prescribed BCs
                readBoundaryConditions(infile, line);
            }
            else if(line.compare(string("Tractions")) == 0)
            {
                cout << "Tractions" << endl;
                // read the number of prescribed BCs
                readTractionConditions(infile, line);
            }
            else if(line.compare(string("Material")) == 0)
            {
                cout << "Material" << endl;
                //readMaterialProps(infile, line);
            }
            else if(line.compare(string("Element Properties")) == 0)
            {
                cout << "Element Properties" << endl;
                readElementProps(infile, line);
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
            else if(line.compare(string("Initial Conditions")) == 0)
            {
                cout << "Initial Conditions" << endl;
                readInitialConditions(infile, line);
            }
            else if(line.compare(string("Nodal Forces")) == 0)
            {
                cout << "Nodal Forces" << endl;
                readNodalForces(infile, line);
            }
            else if(line.compare(string("Output")) == 0)
            {
                cout << "Output" << endl;
                //readOutputDetails(infile, line);
            }
            else if(line.compare(string("Patch Output")) == 0)
            {
                cout << "\n\n Reading Patch Output \n\n" << endl;
                readOutputDetailsPatch(infile, line);
            }
            else
            {
                cout << "key =  " <<  line << endl;
                throw runtime_error("Key not found in femINSmixed::readConfiguration ...");
                //return -1;
            }
      }
    }

    infile.close();

    cout << " Configuration file is successfully read " << endl;

    return 0;
}






int femINSmixed::readFluidProperties(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if( (stringlist[0] == "density") || (stringlist[0] == "rho") )
            {
                fluidProperties[0] = stod(stringlist[1]);
            }
            else if( (stringlist[0] == "viscosity") || (stringlist[0] == "mu") )
            {
                fluidProperties[1] = stod(stringlist[1]);
            }
            else
            {
                throw runtime_error("Option not available in femINSmixed::readFluidProperties");
            }
        }
    }


    return 0;
}







int femINSmixed::readPrescribedBCs(ifstream& infile, string& line)
{
    // read the number of prescribed BCs
    getline(infile,line);    boost::trim(line);

    nDBC = stoi(line);

    cout << "nDBCs = " << nDBC << endl;

    vector<string>  stringlist;
    vector<double>  vecTemp(3);

    nDBC_Pres = 0.0;
    nDBC_Velo = 0.0;

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

        vecTemp[0] = stoi(stringlist[0])-1;
        vecTemp[1] = stoi(stringlist[1])-1;
        vecTemp[2] = stod(stringlist[2]);

        if( vecTemp[1] == ndim )
        {
          DirichletBCsPres.push_back(vecTemp);
          nDBC_Pres++;
        }
        else
        {
          DirichletBCsVelo.push_back(vecTemp);
          nDBC_Velo++;
        }
        //printVector(DirichletBCs[ee]);
    }

    cout << " nDBC_Velo      = " << '\t' << nDBC_Velo << endl;
    cout << " nDBC_Pres      = " << '\t' << nDBC_Pres << endl;

    return 0;
}



int femINSmixed::readNodalForces(ifstream& infile, string& line)
{
    // read the number of nodal forces
    getline(infile,line);    boost::trim(line);

    nFBC = stoi(line);

    cout << "nFBC = " << nFBC << endl;

    vector<string>  stringlist;
    vector<double>  vecTemp(3);

    for(int ee=0; ee<nFBC; ee++)
    {
        getline(infile,line);        boost::trim(line);
        boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
        //stringlist.erase(std::remove(stringlist.begin(), stringlist.end(), " "), stringlist.end());


        int ind = stringlist.size();
        //cout << "ind = " << ind << endl;

        vecTemp[0] = stod(stringlist[0])-1;
        vecTemp[1] = stod(stringlist[1])-1;
        vecTemp[2] = stod(stringlist[2]);

        nodeForcesData.push_back(vecTemp);

        //printVector(DirichletBCs[ee]);
    }

    return 0;
}



int femINSmixed::readElementProps(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist[0] == "type")
            {
                ELEMENT_TYPE = stringlist[1];
            }
            else
            {
                throw runtime_error("Option not available in femINSmixed::readElementProps");
            }
        }
    }

    return 0;
}




int femINSmixed::readFluidProps(ifstream& infile, string& line)
{
    // read {
    getline(infile,line);    boost::trim(line);

    vector<string>  stringlist;
    // read name
    getline(infile,line);    boost::trim(line);
    boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);


    //std::string matkey(stringlist[1]);
    //int  idd = getMaterialID(matkey);
    //cout << " idd = " << idd << endl;

    //FluidMatlDataList.push_back( NewMaterial(idd) );
    FluidMatlDataList.push_back( (FluidMaterialBase*) new FluidMaterialBase );

    FluidMatlDataList[numMaterials]->setID(0);

    //FluidMatlDataList[numMaterials]->SolnData = &(SolnData);

    FluidMatlDataList[numMaterials]->readInput(infile, line);

    numMaterials++;

    cout << "numMaterials = " << numMaterials << endl;

    return 0;
}





int femINSmixed::readOutputDetails(ifstream& infile, string& line)
{
/*
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
                throw runtime_error("Option not available in femINSmixed::readOutputDetails");
            }
        }
    }
*/

    return 0;
}







int  femINSmixed::readBodyForce(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        if(line[0] == '}') break;

        cout << line << endl;

        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            if(stringlist[0] == "value")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
                for(auto& str: stringlist)  boost::trim(str);

                bodyForce[0] = stod(stringlist[0]);
                bodyForce[1] = stod(stringlist[1]);
                bodyForce[2] = stod(stringlist[2]);
            }
            else if( (stringlist[0] == "timefunction") || (stringlist[0] == "TimeFunction") )
            {
                bodyForceTimeFunction = stoi(stringlist[1])-1;
            }
            else
            {
                throw runtime_error("Option not available in femINSmixed::readBodyForce");
            }
        }
    }


    return 0;
}





int  femINSmixed::readSolverDetails(ifstream& infile, string& line)
{
    vector<string>  stringlist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            cout << stringlist[0] << '\t' << stringlist[1] << endl;

            if(stringlist[0] == "schemetype")
            {
                if(stoi(stringlist[1]) == 0)
                  SCHEME_TYPE = SCHEME_TYPE_IMPLICIT;
                else if(stoi(stringlist[1]) == 1)
                  SCHEME_TYPE = SCHEME_TYPE_SEMIIMPLICIT;
                else if(stoi(stringlist[1]) == 2)
                  SCHEME_TYPE = SCHEME_TYPE_EXPLICIT;
            }
            else if(stringlist[0] == "timescheme")
            {
                timeIntegrationScheme = stringlist[1];
            }
            else if(stringlist[0] == "spectralRadius")
            {
                spectralRadius = stod(stringlist[1]);
            }
            else if(stringlist[0] == "finalTime")
            {
                timeFinal = stod(stringlist[1]);
            }
            else if(stringlist[0] == "timeStep")
            {
                line = stringlist[1];
                boost::algorithm::split(stringlist, line, boost::is_any_of("\t "), boost::token_compress_on);
                for(auto& str: stringlist)  boost::trim(str);

                if(stringlist.size() == 1)
                {
                    myTime.set( stod(stringlist[0]), stod(stringlist[0]), stod(stringlist[0]) );
                }
                else if(stringlist.size() == 1)
                {
                    myTime.set( stod(stringlist[0]), stod(stringlist[1]), stod(stringlist[0]) );
                }
                else
                {
                    myTime.set( stod(stringlist[0]), stod(stringlist[1]), stod(stringlist[2]) );
                }
            }
            else if(stringlist[0] == "maximumSteps")
            {
                stepsMax = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "maximumIterations")
            {
                iterationsMax = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "tolerance")
            {
                conv_tol = stod(stringlist[1]);
            }
            else if(stringlist[0] == "debug")
            {
                debug = ( stoi(stringlist[1]) == 1);
            }
            else if(stringlist[0] == "outputFrequency")
            {
                outputFreq = stoi(stringlist[1]);
            }
            else if(stringlist[0] == "CFL")
            {
                CFL = stod(stringlist[1]);
            }
            else
            {
                throw runtime_error("Option not available in femINSmixed::readSolverDetails");
            }
        }
    }


    return 0;
}





int femINSmixed::readTimeFunctions(ifstream& infile, string& line)
{
    // read {
    //getline(infile,line);    boost::trim(line);
    //prgReadTimeFunctions(infile);

    vector<string>  stringlist;
    vector<double>  doublelist;

    // read {
    getline(infile,line);    boost::trim(line);

    while( infile && (line != "}") )
    {
        getline(infile,line);    boost::trim(line);

        cout << line << endl;

        if(line[0] == '}') break;


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
        {
            boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            int id = stoi(stringlist[0]);

            if(id == 0)
            {
                cout << " Time function id '0' is not addmissible. It should be 1 or higher. " << endl;
                throw runtime_error(" Error in femSolidmechanics::readTimeFunctions...  ");
            }

            if(stringlist.size() < 11)
            {
                cout << " Number of temrs in Time function should be at least 11 " << endl;
                throw runtime_error(" Error in femSolidmechanics::readTimeFunctions...  ");
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

    return 0;
}





int femINSmixed::readDomains(ifstream& infile, string& line)
{
/*
    vector<string>  stringlist;

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

            string  label = stringlist[0];

            int index=-1, ii;
            for(ii=0; ii<numDomains; ++ii)
            {
                if( label == Domains[ii]->getLabel() )
                {
                    index = ii;
                    break;
                }
            }

            if(index == -1)
            {
                throw runtime_error("Domain label in femINSmixed::readDomains");
            }

            Domains[index]->readData(infile, line);
        }
    }
*/
    return 0;
}




int  femINSmixed::readBoundaryConditions(ifstream& infile, string& line)
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


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
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

                    cout << "Inside patch ... " << stringlist[0] << endl;

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
                        throw runtime_error("Option not available in femINSmixed::readBoundaryConditions");
                    }
                }//if
            } //while
        }//if
    }//while

    return 0;
}








int  femINSmixed::readTractionConditions(ifstream& infile, string& line)
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


        if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
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

                if( (line.size() > 1) && !( (line[0] == '!') || (line[0] == '#') || (line[0] == '%') ) )
                {
                    boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
                    for(auto& str: stringlist)  boost::trim(str);

                    cout << "Inside patch ... " << stringlist[0] << endl;

                    if(stringlist[0] == "type")
                    {
                        tractionConitions[numTractionConditions]->BCType = stringlist[1];
                    }
                    else if(stringlist[0] == "value")
                    {
                        assert(stringlist.size() >= ndim);

                        for(int i=0; i<ndim; i++)
                          tractionConitions[numTractionConditions]->specValues[i] = stod(stringlist[1+i]);
                    }
                    else if(stringlist[0] == "timefunction")
                    {
                        tractionConitions[numTractionConditions]->timeFunctionNum   = stoi(stringlist[1])-1;
                    }
                    else
                    {
                        throw runtime_error("Option not available in femINSmixed::readTractionConditions");
                    }
                }//if
            } //while
        }//if
    }//while

    return 0;
}





int femINSmixed::readInitialConditions(ifstream& infile, string& line)
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
            boost::algorithm::split(stringlist, line, boost::is_any_of(":"), boost::token_compress_on);
            for(auto& str: stringlist)  boost::trim(str);

            label = stringlist[0];

            int dof = getDOFfromString(stringlist[0]);

            cout << "dof = " << dof << '\t' << ndof << endl;

            if( dof >= ndof)
            {
                throw runtime_error("Specified DOF higher than #DOF in femINSmixed::readInitialConditions ...");
            }

            initialConitions[dof] = stringlist[1];
        }
    }

    return 0;
}







int  femINSmixed::readOutputDetailsPatch(ifstream& infile, string& line)
{
    vector<string>  stringlist;

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

            int ind=1;
            for(auto& bpt : BoundaryPatches)
            {
               if( stringlist[0] == bpt->getLabel() )
               {
                   bpt->setOutputFlag();
                   break;
               }
               ind++;
            }

            if(ind > numBoundaryPatches)
            {
                throw runtime_error("Patch not available in femINSmixed::readOutputDetailsPatch");
            }
        }
    }

    return 0;
}





int femINSmixed::prepareInputData()
{
    printf("\n     femINSmixed::prepareInputData()  .... STARTED ...\n");

    int ii, jj, kk, ee, nn, ind, n1, n2, dof;

    assert(ndim > 0 && ndim < 4);

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();
    
    npElemVelo = elemConn[0].size()-5;

    cout << "ELEMENT_TYPE = " << '\t' << ELEMENT_TYPE << endl;

    if( (ELEMENT_TYPE == "P2P1") || (ELEMENT_TYPE == "p2p1") )
    {
      if(ndim == 2)
      {
          npElemVelo = 6;
          npElemPres = 3;
      }
      else
      {
          npElemVelo = 10;
          npElemPres = 4;
      }
    }
    else if( (ELEMENT_TYPE == "P2bP1") || (ELEMENT_TYPE == "p2bp1") )
    {
      if(ndim == 2)
      {
          npElemVelo = 7;
          npElemPres = 3;
      }
      else
      {
          npElemVelo = 15;
          npElemPres = 4;
      }
    }
    else if( (ELEMENT_TYPE == "q2q1") || (ELEMENT_TYPE == "q2q1") )
    {
      if(ndim == 2)
      {
          npElemVelo = 9;
          npElemPres = 4;
      }
      else
      {
          npElemVelo = 27;
          npElemPres = 8;
      }
    }
    else
    {
        throw  runtime_error("Invalid ELEMENT_TYPE in femINSmixed::prepareInputData");
    }

    cout << " npElemVelo        =  " << npElemVelo << endl;
    cout << " npElemPres        =  " << npElemPres << endl;


    nodeCoordsCur.resize(nNode_global);
    for(int ii=0; ii<nNode_global; ii++)
        nodeCoordsCur[ii] = nodeCoords[ii];

    nNode_Velo = nNode_global;
    nNode_Pres = nNode_global;

    if( (npElemVelo == 7) )
      nNode_Pres = nElem_global*npElemPres;


    ///////////////////////////////////////////////////////////////////
    //
    ///////////////////////////////////////////////////////////////////

    // create elements and prepare element data
    elems = new ElementBase* [nElem_global];

    vector<int>  nodeNums, nodeNumsPres;

    for(ee=0;ee<nElem_global;++ee)
    {
      //elems[ee] = NewLagrangeElement(SolnData.ElemProp[elemConn[ee][0]].id);

      if(npElemVelo == 6)
        elems[ee] = new BernsteinElem2DINSTria6Node;
      else if(npElemVelo == 7)
        elems[ee] = new BernsteinElem2DINSTriaP2bP1dc;
      else if(npElemVelo == 9)
        elems[ee] = new BernsteinElem2DINSQuad9Node;
      else if(npElemVelo == 10)
        elems[ee] = new BernsteinElem3DINSTetra10Node;


      nodeNums.resize(npElemVelo);
      for(ii=0; ii<npElemVelo; ii++)
        nodeNums[ii] = elemConn[ee][5+ii]-1;

      elems[ee]->elenum   = ee;
      elems[ee]->nodeNums = nodeNums;

      //printVector(elems[ee]->nodeNums);

      //elems[ee]->FluidMatlData = FluidMatlDataList[0];

      elems[ee]->prepareElemData(nodeCoords);
    }

    cout << " elements are created and prepared " << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // find mid nodes and modify the nodal coordinates and Dirichlet BCs
    //
    ///////////////////////////////////////////////////////////////////

    // create the arrays for proce midnodes
    // index 0 --- 1 - midnode, 0 - not midnode
    // index 1 --- connecting node 1
    // index 2 --- connecting node 2

    midNodeData.resize(nNode_global);
    for(ii=0; ii<nNode_global; ++ii)
    {
      midNodeData[ii].resize(3);

      // set the default value to 0 ('not a mid-node')
      midNodeData[ii][0] = 0;  midNodeData[ii][1] = 0;  midNodeData[ii][2] = 0;
    }

    if(npElemVelo == 6)
    {
      for(ee=0; ee<nElem_global; ++ee)
      {
        nodeNums = elems[ee]->nodeNums;

        ii = nodeNums[3];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[0];
        midNodeData[ii][2] = nodeNums[1];

        ii = nodeNums[4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[1];
        midNodeData[ii][2] = nodeNums[2];

        ii = nodeNums[5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[2];
        midNodeData[ii][2] = nodeNums[0];

        nodeNumsPres = elems[ee]->nodeNumsPres;

        pressure_nodes.push_back(nodeNumsPres[0]);
        pressure_nodes.push_back(nodeNumsPres[1]);
        pressure_nodes.push_back(nodeNumsPres[2]);
      }
    }
    else if(npElemVelo == 7)
    {
      for(ee=0; ee<nElem_global; ++ee)
      {
        nodeNums = elems[ee]->nodeNums;

        ii = nodeNums[3];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[0];
        midNodeData[ii][2] = nodeNums[1];

        ii = nodeNums[4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[1];
        midNodeData[ii][2] = nodeNums[2];

        ii = nodeNums[5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[2];
        midNodeData[ii][2] = nodeNums[0];

        nodeNumsPres = elems[ee]->nodeNumsPres;

        pressure_nodes.push_back(nodeNumsPres[0]);
        pressure_nodes.push_back(nodeNumsPres[1]);
        pressure_nodes.push_back(nodeNumsPres[2]);
      }
    }
    else if(npElemVelo == 9)
    {
      for(ee=0; ee<nElem_global; ++ee)
      {
        nodeNums = elems[ee]->nodeNums;

        ii = nodeNums[4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[0];
        midNodeData[ii][2] = nodeNums[1];

        ii = nodeNums[5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[1];
        midNodeData[ii][2] = nodeNums[2];

        ii = nodeNums[6];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[2];
        midNodeData[ii][2] = nodeNums[3];

        ii = nodeNums[7];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[3];
        midNodeData[ii][2] = nodeNums[0];

        nodeNumsPres = elems[ee]->nodeNumsPres;

        pressure_nodes.push_back(nodeNumsPres[0]);
        pressure_nodes.push_back(nodeNumsPres[1]);
        pressure_nodes.push_back(nodeNumsPres[2]);
        pressure_nodes.push_back(nodeNumsPres[3]);
      }
    }
    else if(npElemVelo == 10)
    {
      for(ee=0; ee<nElem_global; ++ee)
      {
        nodeNums = elems[ee]->nodeNums;

        ii = nodeNums[4];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[0];
        midNodeData[ii][2] = nodeNums[1];

        ii = nodeNums[5];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[1];
        midNodeData[ii][2] = nodeNums[2];

        ii = nodeNums[6];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[0];
        midNodeData[ii][2] = nodeNums[2];

        ii = nodeNums[7];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[3];
        midNodeData[ii][2] = nodeNums[0];

        ii = nodeNums[8];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[1];
        midNodeData[ii][2] = nodeNums[3];

        ii = nodeNums[9];
        midNodeData[ii][0] = 1;
        midNodeData[ii][1] = nodeNums[2];
        midNodeData[ii][2] = nodeNums[3];


        nodeNumsPres = elems[ee]->nodeNumsPres;

        pressure_nodes.push_back(nodeNumsPres[0]);
        pressure_nodes.push_back(nodeNumsPres[1]);
        pressure_nodes.push_back(nodeNumsPres[2]);
        pressure_nodes.push_back(nodeNumsPres[3]);
      }
    }

    findUnique(pressure_nodes);

    nNode_Pres = pressure_nodes.size();

    pressure_nodes_map.resize(nNode_global,-1);
    for(ii=0; ii<nNode_Pres; ++ii)
    {
      pressure_nodes_map[pressure_nodes[ii]] = ii;
    }

    cout << " nNode_Pres = " << nNode_Pres << endl;

    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    ind = nNode_Velo*ndim;

    velo.resize(ind);
    velo.setZero();

    veloPrev  = velo;
    veloCur   = velo;
    veloPrev2 = velo;
    veloApplied = velo;

    veloDot     = velo;
    veloDotPrev = veloDot;
    veloDotCur  = veloDot;

    if(npElemVelo == 7)
      pres.resize(nNode_Pres);
    else
      pres.resize(nNode_Velo);

    pres.setZero();

    presPrev    = pres;
    presPrev2   = pres;
    presCur     = pres;
    presApplied = pres;
    presDot     = velo;
    presDotPrev = pres;

    double  xx, yy, zz, fact;

    cout << " ========== " << endl;

    // loop over the nodes and adjust nodal coordinates
    //
    for(nn=0; nn<nNode_global; nn++)
    {
        if( midNodeData[nn][0] )
        {
          n1 = midNodeData[nn][1];
          n2 = midNodeData[nn][2];

          xx = 0.25*nodeCoords[n1][0] + 0.25*nodeCoords[n2][0];
          yy = 0.25*nodeCoords[n1][1] + 0.25*nodeCoords[n2][1];

          nodeCoords[nn][0] = 2.0*(nodeCoords[nn][0] - xx);
          nodeCoords[nn][1] = 2.0*(nodeCoords[nn][1] - yy);

          if(ndim == 3)
          {
            zz = 0.25*nodeCoords[n1][2] + 0.25*nodeCoords[n2][2];

            nodeCoords[nn][2] = 2.0*(nodeCoords[nn][2] - zz);
          }
        }
    }
    //

    printf("     femINSmixed::prepareInputData()  .... FINISHED ...\n\n");

    return 0;
}



int  femINSmixed::readImmersedSolids(string& fname)
{
    cout << " Reading Immersed elements \n\n " << endl;

    ImmersedSolid  *imsolid = new ImmersedSolid;

    imsolid->readInput(fname);
    
    ImmersedBodyObjects.push_back(imsolid);    


    cout << " no. of immersed Nodes              =  " << ImmersedBodyObjects[0]->getNumberOfNodes() << endl;
    cout << " no. of immersed Elements           =  " << ImmersedBodyObjects[0]->getNumberOfElements() << endl;

    return 0;
}




int  femINSmixed::prepareImmersedSolids()
{
    if(ImmersedBodyObjects.size() == 0)
      return 0;

    printf("\n     femINSmixed::prepareImmersedSolids()  .... STARTED ...\n");

    ///////////////////////////////////////////////////
    //
    // prepare the data for immersed points
    // 
    ///////////////////////////////////////////////////

    int iie, bb, gp, elnum, ind1, ind2;
    myPoint  geom;

    ImmersedIntegrationElement *lme;
    vector<int> nodeNums;

    for(bb=0;bb<ImmersedBodyObjects.size();bb++)
    {
      for(iie=0; iie<ImmersedBodyObjects[bb]->getNumberOfElements(); iie++)
      {
        lme = ImmersedBodyObjects[bb]->ImmIntgElems[iie];

        //cout << bb << '\t' << iie << '\t' << lme->isActive() << '\t' << lme->gausspoints.size() << endl;

        nodeNums = lme->nodeNums;

        if( lme->isActive() )
        {
          for(int ii=0; ii<nodeNums.size(); ii++)
          {
            ind1 = ii*ndim;
            ind2 = nodeNums[ii]*ndim;

            for(int jj=0; jj<ndim; jj++)
              //lme->forAssyVec[ind1+jj] = totalDOF_Lambda++;
              lme->forAssyVec[ind1+jj] = ind2+jj;
          }

          for(gp=0; gp<lme->gausspoints.size(); gp++)
          {
            //lme->computePointAtGP(gp, geom);
            ImmersedBodyObjects[bb]->computePointAtGP(iie, lme->gausspoints[gp], geom);

            //cout << " finding Background element " << endl;

            elnum = findElementNumber(geom);

            lme->elemNums[gp] = elnum;

            //printf("xx = %12.6f, yy = %12.6f, zz = %12.6f, dvol = %5d, \n", geom[0], geom[1], geom[2], elnum);

            //lme->elems[gp] = elems[elnum];
          }
        }
      }
    }

    totalDOF_Lambda = ImmersedBodyObjects[0]->getNumberOfNodes()*ndim;

    lambdas.resize(totalDOF_Lambda);
    lambdasPrev = lambdas;
    lambdasIncr = lambdas;
    lambdasCur  = lambdas;

    totalDOF_Solid = 0;

    printf("\n     femINSmixed::prepareImmersedSolids()  .... FINISHED ...\n");

    return 0;
}





int femINSmixed::printInfo()
{
/*
  printf("\n Background Fluid Grid information \n");
    printf("\n ----------------------------------------------------------------------------- \n");
    printf("\t   Problem   Dimension        =  %5d\n\n", DIM);
    printf("\t   Polynomial Degree          =  %5d\t%5d\t%5d\n\n", degree[0], degree[1], degree[2]);
    printf("\t   Number of Elements         =  %5d\t%5d\t%5d\n\n", nelem[0], nelem[1], nelem[2]);
    printf("\t   Origin of the grid         =  %12.6f\t%12.6f\t%12.6f\n\n", origin[0], origin[1], origin[2]);
    printf("\t   Grid Length                =  %12.6f\t%12.6f\t%12.6f\n\n", gridLEN[0], gridLEN[1], gridLEN[2]);
    printf("\t   MAXIMUM REFINEMENT LEVEL   =  %5d\n\n", MAX_LEVEL);
    printf("\t   DOF at each Control Point  =  %5d\n\n", ndof);
    printf("\t   Total number of DOF        =  %5d\n\n", totalDOF);
    printf("\n ----------------------------------------------------------------------------- \n");
*/
    return 0;
}







int femINSmixed::setSpecifiedDOFs_Velocity(vector<vector<bool> >&  NodeDofType)
{
  //if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femINSmixed::setSpecifiedDOFs_Velocity() ... STARTED \n");
  if(debug)  printf("femINSmixed::setSpecifiedDOFs_Velocity() ... STARTED \n");

  vector<int>  nodeNums;
  int ii, jj, nn, ind, dof;

  dofs_specified_velo.clear();


    // set specified DOFs
    for(auto&  bc : boundaryConitions)
    {
        if(debug)
        {
        cout << " Patch                 = " << bc->label << endl;
        cout << " BCType                = " << bc->BCType << endl;
        cout << " dof_specified_string  = " << bc->dof_specified_string << endl;
        cout << " dof_specified_int     = " << bc->dof_specified_int << endl;
        cout << " expression            = " << bc->expression << endl;
        cout << " timeFunctionNum       = " << bc->timeFunctionNum << endl;
        cout << endl;  cout << endl;
        }
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
                dofs_specified_velo.push_back(nn*ndof+dof);
            }
        }
        else if(bc->BCType == "wall")
        {
            for(ii=0; ii<nodeNums.size(); ++ii)
            {
              nn = nodeNums[ii];
              ind = nn*ndof;
              for(dof=0; dof<ndof; ++dof)
              {
                NodeDofType[nn][dof] = true;
                dofs_specified_velo.push_back(ind+dof);
              }
            }
        }
        else
        {
            throw runtime_error("Patch type not available in femINSmixed::setSpecifiedDOFs_Velocity");
        }
    }

    findUnique(dofs_specified_velo);
    //printVector(dofs_specified_velo);

    //if(debug)  PetscPrintf(PETSC_COMM_WORLD, "femINSmixed::setSpecifiedDOFs_Velocity() ... ENDED \n");
    if(debug)  printf("femINSmixed::setSpecifiedDOFs_Velocity() ... ENDED \n");

    return 0;
}












int femINSmixed::setBoundaryConditions()
{
    //if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femINSmixed::setBoundaryConditions() ... STARTED \n");}
    if(debug) {printf("\n\n femINSmixed::setBoundaryConditions() ... STARTED \n");}

    int  ii, nn, dof, index=-1, timeFuncNum, ind;
    double xc, yc, zc, value, loadFactor;
    vector<int> nodeNums;

    veloApplied.setZero();
    presApplied.setZero();

    // set boundary conditions
    for(auto&  bc : boundaryConitions)
    {
        //cout << " Patch                 = " << bc->label << endl;
        //cout << " BCType                = " << bc->BCType << endl;
        //cout << " dof_specified_string  = " << bc->dof_specified_string << endl;
        //cout << " dof_specified_int     = " << bc->dof_specified_int << endl;
        //cout << " expression            = " << bc->expression << endl;
        //cout << " timeFunctionNum       = " << bc->timeFunctionNum << endl;
        //cout << endl;  cout << endl;

        // nodeNums contains new node numbers which are used for the solution
        // nodeCoordsOrig contains coordinates of nodes with numbers before domain decompositions

        nodeNums = bc->bpt->nodeNums;

        //printVector(nodeNums);

        if(bc->BCType == "specified")
        {
            timeFuncNum = bc->getTimeFunctionNumber();

            if(timeFuncNum == -1)
              loadFactor = 1.0;
            else
              loadFactor = timeFunctions[timeFuncNum]->getValue();

            //PetscPrintf(MPI_COMM_WORLD, " bc->label = %10s ...  bc->timeFunctionNum = %d ... load Factor = %12.6f \n", bc->label.c_str(), bc->timeFunctionNum,  loadFactor);
            //printf(" bc->label = %10s ...  bc->timeFunctionNum = %d ... load Factor = %12.6f \n", bc->label.c_str(), bc->timeFunctionNum,  loadFactor);

            myMathFunction  mathfun;
            mathfun.initialise(bc->expression);

            dof = bc->dof_specified_int;

            for(ii=0; ii<nodeNums.size(); ++ii)
            {
                //nn = node_map_get_old[nodeNums[i]];
                nn = nodeNums[ii];

                xc = nodeCoords[nn][0];
                yc = nodeCoords[nn][1];
                zc = nodeCoords[nn][2];

                value = mathfun.getValue(xc, yc, zc, myTime.cur) * loadFactor;

                //cout << xc << '\t' << yc << '\t' << zc << '\t' << loadFactor << '\t' << value << endl;

                veloApplied[nn*ndof+dof] = value;
            }
        }
        else if(bc->BCType == "wall")
        {
            for(ii=0; ii<nodeNums.size(); ++ii)
            {
              nn = nodeNums[ii];
              ind = nn*ndof;

              for(dof=0; dof<ndim; ++dof)
              {
                veloApplied[ind+dof] = 0.0;
              }
            }
        }
        else
        {
            throw runtime_error("Patch type not available in femSolidmechanics::setBoundaryConditions");
        }
    }

    //if(SCHEME_TYPE == IMPLICIT)
    //{
      veloApplied -= velo;
      presApplied -= pres;
    //}
    //printVector(veloApplied);

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

    //if(debug) {PetscPrintf(MPI_COMM_WORLD, "\n\n femINSmixed::setBoundaryConditions() ... ENDED \n");}
    if(debug) {printf("\n\n femINSmixed::setBoundaryConditions() ... ENDED \n");}

    return 0;
}




/*
int femINSmixed::assignBoundaryConditions()
{
    int ii, jj, nn, dof;

    for(ii=0; ii<nDBC_Velo; ++ii)
    {
        nn  = DirichletBCsVelo[ii][0];
        dof = DirichletBCsVelo[ii][1];

        jj = nn*ndim+dof;

        //cout << n1 << '\t' << n2 << '\t' << DirichletBCsVelo[ii][2] << '\t' << loadFactor << endl;

        veloApplied[jj] = DirichletBCsVelo[ii][2] * loadFactor - velo[jj];

        //veloApplied[jj] = analy.computeValue(n2, nodeCoords[n1][0], nodeCoords[n1][1], 0.0, timeNow) - velo[jj];
    }

    double fact;
    for(ii=0; ii<nDBC_Velo; ++ii)
    {
        nn  = DirichletBCsVelo[ii][0];
        dof = DirichletBCsVelo[ii][1];

        if( midNodeData[nn][0] )
        {
          fact = 0.25*veloApplied(midNodeData[nn][1]*ndim+dof) + 0.25*veloApplied(midNodeData[nn][2]*ndim+dof);

          veloApplied[nn*ndim+dof] = 2.0*(veloApplied(nn*ndim+dof) - fact);
        }
    }


    for(ii=0; ii<nDBC_Pres; ++ii)
    {
        jj = DirichletBCsPres[ii][0];

        presApplied[jj] = DirichletBCsPres[ii][2] * loadFactor - pres[jj];

        //presApplied[jj] = analy.computeValue(2, nodeCoords[jj][0], nodeCoords[jj][1], 0.0, timeNow) - pres[jj];
    }

    return 0;
}
*/

/*
int femINSmixed::applyBoundaryConditions()
{
    //cout <<  " applying boundary conditions .... " << endl;

    int ii, jj, n1, n2;
    double  fact1 = 1.0/(gamm1*dt), fact2 = (gamm1-1.0)/gamm1;

    for(ii=0; ii<nDBC_Velo; ++ii)
    {
        n1 = DirichletBCsVelo[ii][0];
        n2 = DirichletBCsVelo[ii][1];

        jj = n1*ndim+n2;

        velo[jj]    = veloApplied[jj] * loadFactor;
        veloDot[jj] = fact1 * (velo[jj]-veloPrev[jj]) + fact2 * veloDotPrev[jj];
    }

    for(ii=0; ii<nDBC_Pres; ++ii)
    {
        jj = DirichletBCsPres[ii][0];

        pres[jj]    = DirichletBCsPres[ii][2] * loadFactor;
        presDot[jj] = fact1 * (pres[jj]-presPrev[jj]) + fact2 * presDotPrev[jj];
    }

    return 0;
}
*/





int femINSmixed::addBoundaryConditions()
{
    //cout <<  " applying boundary conditions .... " << endl;

    // add specified Dirichlet boundary conditions if first iteration
        int ii, dof;
        for(ii=0; ii<dofs_specified_velo.size(); ii++)
        {
            dof = dofs_specified_velo[ii];
            velo[dof] += veloApplied[dof] ;
        }

        for(ii=0; ii<dofs_specified_pres.size(); ii++)
        {
            dof = dofs_specified_pres[ii];
            pres[dof] += presApplied[dof] ;
        }

    return 0;
}




int femINSmixed::setInitialConditions()
{
    double  xx=0.0, yy=0.0, zz=0.0, fact;
    //double  specVal;

    // adjust the velocity values for the midnoes
    //VectorXd  velTemp;

    /*
    for(nn=0; nn<nNode_global; nn++)
    {
      if( midNodeData[nn][0] ) // if midnode
      {
        for(dd=0; dd<ndof; dd++)
        {
          xx = 0.25*velTemp(midNodeData[nn][1]*ndof+dd) + 0.25*velTemp(midNodeData[nn][2]*ndof+dd);

          //SolnData.var1Dot[nn*ndof+dd] = 2.0*(velTemp(nn*ndof+dd) - xx);
        }
        //cout << velTemp(nn*ndof) << '\t' << SolnData.var1Dot(nn*ndof) << endl;
      }
    }
    */

    /*
    for(int ii=0; ii<nNode_Velo; ++ii)
    {
        xx = nodeCoords[ii][0];
        yy = nodeCoords[ii][1];
        zz = nodeCoords[ii][2];

        //veloPrev(ii*2) = 2.0*yy*(3.0-yy)/3.0;
        veloPrev(ii*ndim) = 1.0*yy;

        //veloPrev(ii*ndim) = 16.0*0.45*yy*zz*(0.41-yy)*(0.41-zz)/0.41/0.41/0.41/0.41;
    }
    velo = veloPrev;
    */

    return 0;
}






int femINSmixed::setTimeParam()
{
  //SolnData.setTimeParam();
  SetTimeParametersFluid(timeIntegrationScheme, spectralRadius, dt, td);

  return 0;
}



int femINSmixed::timeUpdate()
{
    //if(debug) {PetscPrintf(MPI_COMM_WORLD, " femSolids::timeUpdate() ... STARTED \n");}
    if(debug) {printf(" femSolids::timeUpdate() ... STARTED \n");}

    myTime.update();

    // set parameters for the time integration scheme.
    // need to be done every time step to account for adaptive time stepping
    SetTimeParametersFluid(timeIntegrationScheme, spectralRadius, myTime.dt, td);

    // update time functions
    for(auto& tmf : timeFunctions)
      tmf->update();

    saveSolution();

    // set Dirichlet boundary conditions
    setBoundaryConditions();
    //printVector(veloApplied);

    //cout << " aaaaaaaaaaaaaaa " << endl;

    forceCur   = td[1]*force   + (1.0-td[1])*forcePrev; // force

    return 0;
}



int femINSmixed::updateIterStep()
{
  //SolnData.updateIterStep();
  veloDot    = td[9]*velo + td[10]*veloPrev + td[11]*veloPrev2 + td[15]*veloDotPrev ;

  veloCur    = td[2]*velo    + (1.0-td[2])*veloPrev; // velocity
  presCur    = td[2]*pres    + (1.0-td[2])*presPrev; // pressure
  lambdasCur = td[2]*lambdas + (1.0-td[2])*lambdasPrev; // Lagrange multipliers
  veloDotCur = td[1]*veloDot + (1.0-td[1])*veloDotPrev;

  return 0;
}



bool femINSmixed::converged()
{
  if (rhsNorm < conv_tol)
    return true;

  return false;
}




int  femINSmixed::saveSolution()
{
    veloPrev2   = veloPrev;
    veloPrev    = velo;
    veloDotPrev = veloDot;
    presPrev    = pres;
    lambdasPrev = lambdas;
    forcePrev   = force;

    return 0;
}



int femINSmixed::reset()
{
    velo      = veloPrev;
    veloPrev  = veloPrev2;
    veloDot   = veloDotPrev;
    pres      = presPrev;
    lambdas   = lambdasPrev;
    force     = forcePrev;

    return 0;
}




int femINSmixed::writeNodalData()
{
  return 0;
}




int femINSmixed::writeOutputDataPatches()
{
    if(debug) printf("femINSmixed::writeOutputDataPatches() ... \n");

    vector<int> nodeNums;
    int  i, nn, dof, ind, size;
    double  totalForce[3], totalMoment[3];

    for(auto& bpt : BoundaryPatches)
    {
      if(bpt->getOutputFlag())
      {
          nodeNums = bpt->nodeNums;
          size = nodeNums.size();

          totalForce[0] = 0.0;
          totalForce[1] = 0.0;
          totalForce[2] = 0.0;

          totalMoment[0] = 0.0;
          totalMoment[1] = 0.0;
          totalMoment[2] = 0.0;

          if(ndim == 2)
          {
            for(i=0; i<size; ++i)
            {
              nn = nodeNums[i];
              totalForce[0] += reacVec[nn*ndof];
              totalForce[1] += reacVec[nn*ndof+1];
            }
          }
          else
          {
            for(i=0; i<size; ++i)
            {
              nn = nodeNums[i];
              totalForce[0] += reacVec[nn*ndof];
              totalForce[1] += reacVec[nn*ndof+1];
              totalForce[2] += reacVec[nn*ndof+2];
            }
          }

          bpt->forcedata << myTime.cur << '\t' << totalForce[0]  << '\t' << totalForce[1]  << '\t' << totalForce[2]
                                       << '\t' << totalMoment[0] << '\t' << totalMoment[1] << '\t' << totalMoment[2] << endl;
      }
    }

    if(debug) printf("femINSmixed::writeOutputDataPatches() ... \n");

/*
    PetscPrintf(MPI_COMM_WORLD, "femINSmixed::writeOutputDataPatches() ... \n");

    PetscScalar *arrayTempReac;
    Vec            vecseqReac;
    VecScatter     ctxReac;

    VecScatterCreateToZero(solverPetsc->reacVec, &ctxReac, &vecseqReac);

    VecAssemblyBegin(solverPetsc->reacVec);
    VecAssemblyEnd(solverPetsc->reacVec);

    VecScatterBegin(ctxReac, solverPetsc->reacVec, vecseqReac, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd(ctxReac,   solverPetsc->reacVec, vecseqReac, INSERT_VALUES, SCATTER_FORWARD);

    VecGetArray(vecseqReac, &arrayTempReac);


    vector<int> nodeNums;
    int  i, nn, dof, ind, size;
    double  totalForce[3], totalMoment[3];

    if(this_mpi_proc == 0)
    {
    for(auto& bpt : mesh->BoundaryPatches)
    {
      if(bpt->getOutputFlag())
      {
          nodeNums = bpt->nodeNums;
          size = nodeNums.size();

          totalForce[0] = 0.0;
          totalForce[1] = 0.0;
          totalForce[2] = 0.0;

          totalMoment[0] = 0.0;
          totalMoment[1] = 0.0;
          totalMoment[2] = 0.0;

          if(mesh->getNDIM() == 2)
          {
            for(i=0; i<size; ++i)
            {
              nn = nodeNums[i];
              totalForce[0] += arrayTempReac[nn*ndof];
              totalForce[1] += arrayTempReac[nn*ndof+1];
            }
          }
          else
          {
            for(i=0; i<size; ++i)
            {
              nn = nodeNums[i];
              totalForce[0] += arrayTempReac[nn*ndof];
              totalForce[1] += arrayTempReac[nn*ndof+1];
              totalForce[2] += arrayTempReac[nn*ndof+2];
            }
          }

          bpt->forcedata << myTime.cur << '\t' << totalForce[0]  << '\t' << totalForce[1]  << '\t' << totalForce[2]
                                       << '\t' << totalMoment[0] << '\t' << totalMoment[1] << '\t' << totalMoment[2] << endl;
      }
    }
    }

    VecRestoreArray(vecseqReac, &arrayTempReac);

    VecScatterDestroy(&ctxReac);
    VecDestroy(&vecseqReac);
    PetscFree(arrayTempReac);
*/
    return 0;
}





int femINSmixed::initialise_pardiso()
{
  soln.resize(totalDOF);
  soln.setZero();

  perm.resize(totalDOF);

  phase = 11; error = 0;
  int  mtxType = 11, idum;

  char *tmp;


  //mtxType =  1;  // real and structurally symmetric
  mtxType = 11; // real and unsymmetric


  SOLVER = 0;       // sparse direct solver
  //SOLVER = 1;       // multi-recursive iterative solver
  MTYPE  = mtxType; // matrix type
  MAXFCT = 1;       // maximum number of factorisations of same sparsity pattern
                    //      to be kept in memory
  MNUM   = 1;       // which factorisation to use
  NRHS   = 1;       // number of right hand sides
  MSGLVL = 0;       // output message level (1 -> print statistical information)

  IPARM[0] = 0;     // PARADISO will set IPARM to default values
  //IPARM[0] = 1;     // user input values to IPARM
  //IPARM[1] = 2;

  tmp = getenv("OMP_NUM_THREADS");

  if (tmp != NULL) 
  {
    sscanf(tmp,"%d", &idum);
  }
  else printf("set environment variable OMP_NUM_THREADS!");


  cout << "OMP_NUM_THREADS = " << idum << endl;

  IPARM[2] = max(1,idum);  // number of processors (no default value available)


  pardisoinit_(PT, &MTYPE, &SOLVER, IPARM, DPARM, &error);

  if (error != 0)
  {
    if (error == -10) printf("no license file found.");
    if (error == -11) printf("license is expired.");
    if (error == -12) printf("wrong username or hostname.");
  }

  int  *c1, *c2, ii;

  csr.resize(totalDOF+1);
  col.resize(matK.nonZeros());

  c1 = matK.outerIndexPtr();
  c2 = matK.innerIndexPtr();


  for(ii=0;ii<=totalDOF;ii++)
    csr[ii] = c1[ii] + 1;

  for(ii=0;ii<matK.nonZeros();ii++)
    col[ii] = c2[ii] + 1;

  array = matK.valuePtr();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &totalDOF, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &ddum, &ddum, &error, DPARM);

  if (error != 0)
  {
    cout << "PARDISO ERROR = " << error << "\n\n";
    printf("symbolic factorisation failed.");
  }

  //IPARM[4] = 0;  // user input permutation
  //IPARM[4] = 2;  // return the permutation
  IPARM[5] = 0; // do not overwrite RHS with solution
  IPARM[7] = 1; // max number of iterative refinement steps

  return 0;
}




int femINSmixed::factoriseAndSolve_pardiso()
{
  phase = 23; error = 0;

  soln.setZero();

  pardiso_(PT, &MAXFCT, &MNUM, &MTYPE, &phase,
           &totalDOF, array, &csr[0], &col[0], &perm[0], &NRHS,
           IPARM, &MSGLVL, &rhsVec[0], &soln[0], &error, DPARM);

  //cout << "error = " << error << endl;

//   printf("Peak memory [kB] phase 1       = %d \n", IPARM[14]);
//   printf("Permanent integer memory [kb]. = %d \n", IPARM[15]);
//   printf("Peak real memory [kB]          = %d \n", IPARM[16]);
//   printf("Number of nonzeros in LU.      = %d \n", IPARM[17]);
//   printf("Gflops for LU factorization.   = %d \n", IPARM[18]);

  return 0;
}

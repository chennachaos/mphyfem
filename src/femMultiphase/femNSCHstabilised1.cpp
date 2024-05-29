
#include <algorithm>
#include <chrono>
#include "femNSCHstabilised.h"
#include "KimMoinFlow.h"
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


femNSCHstabilised::femNSCHstabilised()
{
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    fileCount = 0;

    computerTimeAssembly = computerTimeSolver = computerTimePostprocess = computerTimePattern = computerTimeTimeLoop = 0.0;

    outputFreq = 1;
    conv_tol   = 1.0e-8;
    spectralRadius = 0.0;
    timeIntegrationScheme = "STEADY";

    infilename = get_current_dir_name();
    vector<string>  stringlist;
    boost::algorithm::split(stringlist, infilename, boost::is_any_of("/"), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);
    infilename = stringlist[stringlist.size()-1];
}


femNSCHstabilised::~femNSCHstabilised()
{
    //cout << " femNSCHstabilised::~femNSCHstabilised() " << this_mpi_proc << endl;

    //deallocatePetscObjects();

    //cout << " femNSCHstabilised::~femNSCHstabilised() " << this_mpi_proc << endl;
}




int  femNSCHstabilised::deallocatePetscObjects()
{
  NavierStokes.deallocatePetscObjects();
  CahnHilliard.deallocatePetscObjects();

  return 0;
}






int  femNSCHstabilised::readInputGMSH(string& fname)
{
    PetscPrintf(MPI_COMM_WORLD, " BaseFEM::readInputGMSH \n");
    
    mesh = make_shared<myMesh>();
    mesh->readInputGMSH(fname);

    ndim = mesh->getNDIM();

    NavierStokes.mesh = mesh;
    CahnHilliard.mesh = mesh;

    PetscPrintf(MPI_COMM_WORLD, " BaseFEM::readInputGMSH \n");

    return 0;
}




int  femNSCHstabilised::readConfiguration(string& fname, string& fnamefluid, string& fnamephase)
{
    NavierStokes.readConfiguration(fnamefluid);

    CahnHilliard.readConfiguration(fnamephase);

    readConfiguration(fname);

    return 0;
}




int femNSCHstabilised::readConfiguration(string& fname)
{
    cout << " femNSCHstabilised::readConfiguration " << endl;

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

            if(line.compare(string("Fluid Properties")) == 0)
            {
                cout << "Fluid Properties" << endl;
                readFluidProperties(infile, line);
            }
            else if(line.compare(string("Body Force")) == 0)
            {
                cout << "Body Force" << endl;
                // read the number of prescribed BCs
                //readBodyForce(infile, line);
            }
            else if(line.compare(string("Boundaries")) == 0)
            {
                cout << "Boundaries" << endl;
                // read the number of prescribed BCs
                readBoundariesData(infile, line);
            }
            else if(line.compare(string("Material Type")) == 0)
            {
                cout << "Material Type" << endl;
                //readMaterialProps(infile, line);
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
                //readInitialConditions(infile, line);
            }
            else if(line.compare(string("Output")) == 0)
            {
                cout << "Output" << endl;
                //readOutputDetails(infile, line);
            }
            else
            {
                cout << "key =  " <<  line << endl;
                throw runtime_error("Key not found in femNSCHstabilised::readConfiguration ...");
                //return -1;
            }
      }
    }

    infile.close();

    cout << " Configuration file is successfully read " << endl;

    return 0;
}









int  femNSCHstabilised::readBoundariesData(ifstream& infile, string& line)
{
    mesh->readBoundariesData(infile, line);

    return 0;
}





int  femNSCHstabilised::readSolverDetails(ifstream& infile, string& line)
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

            if(stringlist[0] == "timescheme")
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
                if(stringlist.size() == 2)
                {
                    myTime.set( stod(stringlist[1]), stod(stringlist[1]), stod(stringlist[1]) );
                }
                else if(stringlist.size() == 3)
                {
                    myTime.set( stod(stringlist[1]), stod(stringlist[2]), stod(stringlist[1]) );
                }
                else
                {
                    myTime.set( stod(stringlist[1]), stod(stringlist[2]), stod(stringlist[3]) );
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
            else
            {
                throw runtime_error("Option not available in femNSCHstabilised::readBodyForce");
            }
        }
    }


    return 0;
}







int femNSCHstabilised::readFluidProperties(ifstream& infile, string& line)
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

            cout << stringlist[0] << "\t" << stringlist[1] << endl;

            if(stringlist.size() < 3)
            {
                throw runtime_error("Insufficient values in femNSCHstabilised::readFluidProperties");
            }

            if(stringlist[0] == "density")
            {
                fluidDensity[0] = stod(stringlist[1]);
                fluidDensity[1] = stod(stringlist[2]);
            }
            else if(stringlist[0] == "viscosity")
            {
                fluidDynamicViscosity[0] = stod(stringlist[1]);
                fluidDynamicViscosity[1] = stod(stringlist[2]);
            }
            else
            {
                throw runtime_error("Option not available in femNSCHstabilised::readFluidProperties");
            }
        }
    }


    return 0;
}






int  femNSCHstabilised::readTimeFunctions(ifstream& infile, string& line)
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
                throw runtime_error(" Error in femCahnHilliard::readTimeFunctions...  ");
            }

            if(stringlist.size() < 11)
            {
                cout << " Number of temrs in Time function should be at least 11 " << endl;
                throw runtime_error(" Error in femCahnHilliard::readTimeFunctions...  ");
            }

            doublelist.resize(stringlist.size()-1);
            for(int i=0; i<stringlist.size()-1; ++i)
            {
                doublelist[i] = stod(stringlist[i+1]);
            }

            id -= 1;
            if(id == timeFunction.size())
            {
                timeFunction.push_back(make_unique<TimeFunction>());
                timeFunction[id]->setID(id);
            }
            timeFunction[id]->addTimeBlock(doublelist);
        }
    }

    for(auto& tmf : timeFunction)
    {
        tmf->printSelf();
        tmf->update();
    }

    return 0;
}






int femNSCHstabilised::prepareInputData()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n  femNSCHstabilised::prepareInputData()  .... STARTED ...\n\n");

    int ii, jj, kk, ee, nn, ind, n1, n2, dof;

    mesh->prepareInputData();
    
    NavierStokes.prepareInputData();
    CahnHilliard.prepareInputData();

    //elem_proc_id     = mesh->elem_proc_id;
    //node_proc_id     = mesh->node_proc_id;
    //node_map_get_old = mesh->node_map_get_old;
    //node_map_get_new = mesh->node_map_get_new;

    // 
    //
    // Set time integration schemes for the fluid and phase problems
    // ==================================================

    NavierStokes.timeIntegrationScheme = timeIntegrationScheme;
    NavierStokes.spectralRadius        = spectralRadius;
    NavierStokes.MULTIPHASEFLOW        = true;

    CahnHilliard.timeIntegrationScheme = timeIntegrationScheme;
    CahnHilliard.spectralRadius        = spectralRadius;
    CahnHilliard.MULTIPHASEFLOW        = true;

    //checkInputData();
    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    PetscPrintf(MPI_COMM_WORLD, "\n\n  femNSCHstabilised::prepareInputData()  .... FINISHED ...\n\n");

    return 0;
}





int femNSCHstabilised::setInitialConditions()
{
  CahnHilliard.setInitialConditions();

  NavierStokes.setInitialConditions();

  return 0;
}







int femNSCHstabilised::setTimeParam()
{
  //setTimeParam();

  return 0;
}




int femNSCHstabilised::timeUpdate()
{
  myTime.update();

  // set parameters for the time integration scheme.
  // need to be done every time step to account for adaptive time stepping
  SetTimeParametersFluid(timeIntegrationScheme, spectralRadius, myTime.dt, NavierStokes.td);
  SetTimeParametersFluid(timeIntegrationScheme, spectralRadius, myTime.dt, CahnHilliard.td);

  // update time functions
  for(auto& tmf : timeFunction)
    tmf->update();


  NavierStokes.setBoundaryConditions();
  CahnHilliard.setBoundaryConditions();

  return 0;
}



int femNSCHstabilised::updateIterStep()
{

  return 0;
}




int  femNSCHstabilised::reset()
{

  return 0;
}



int  femNSCHstabilised::saveSolution()
{

  return 0;
}




int femNSCHstabilised::writeOutputData()
{
    PetscPrintf(MPI_COMM_WORLD, "femNSCHstabilised::writeOutputData() ... \n");
/*
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






int femNSCHstabilised::printComputerTimes()
{
    PetscPrintf(MPI_COMM_WORLD, "\n====================================================================\n");
    PetscPrintf(MPI_COMM_WORLD, "\n       Computer time in seconds \n");
    PetscPrintf(MPI_COMM_WORLD, "\n====================================================================\n");

    PetscPrintf(MPI_COMM_WORLD, "\n\n     (1) Matrix Pattern                  = %12.6f \n", computerTimePattern );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (2) Matrix Assembly                 = %12.6f \n", computerTimeAssembly );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (3) Matrix Solver                   = %12.6f \n", computerTimeSolver );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (4) Post process                    = %12.6f \n", computerTimePostprocess );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (5) Total time in time loop         = %12.6f \n", computerTimeTimeLoop );
    PetscPrintf(MPI_COMM_WORLD, "\n\n     (6) Difference (5-2-3-4)            = %12.6f \n", (computerTimeTimeLoop-computerTimeAssembly-computerTimeSolver-computerTimePostprocess) );

    PetscPrintf(MPI_COMM_WORLD, "\n\n====================================================================\n\n");

    return 0;
}






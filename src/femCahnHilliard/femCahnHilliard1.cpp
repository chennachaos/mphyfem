
#include <algorithm>
#include <chrono>
#include "femCahnHilliard.h"
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


femCahnHilliard::femCahnHilliard()
{
    MPI_Comm_size(MPI_COMM_WORLD, &n_mpi_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &this_mpi_proc);

    ndof = 2; fileCount = 0;
    ntotdofs_local = ntotdofs_global = 0;
    numBoundaryConditions = 0;

    computerTimeAssembly = computerTimeSolver = computerTimePostprocess = computerTimePattern = computerTimeTimeLoop = 0.0;

    AlgoType = 2;
    bodyForceTimeFunction = -1;
    bodyForce[0] = 0.0;
    bodyForce[1] = 0.0;
    bodyForce[2] = 0.0;

    outputFreq = 1;
    conv_tol   = 1.0e-8;
    spectralRadius = 0.0;
    timeIntegrationScheme = "STEADY";

    SCHEME_TYPE = 2;
    
    MULTIPHASEFLOW = false;

    infilename = get_current_dir_name();
    vector<string>  stringlist;
    boost::algorithm::split(stringlist, infilename, boost::is_any_of("/"), boost::token_compress_on);
    for(auto& str: stringlist)  boost::trim(str);
    infilename = stringlist[stringlist.size()-1];

    td.resize(100);

    stabilisationFactors[0] = 1.0; //SUPG
    stabilisationFactors[1] = 1.0; //PSPG
    stabilisationFactors[2] = 0.0; //LSIC
}


femCahnHilliard::~femCahnHilliard()
{
    //cout << " femCahnHilliard::~femCahnHilliard() " << this_mpi_proc << endl;

    //deallocatePetscObjects();

    //cout << " femCahnHilliard::~femCahnHilliard() " << this_mpi_proc << endl;
}




int  femCahnHilliard::deallocatePetscObjects()
{
  solverPetsc->free();

  return 0;
}






int  femCahnHilliard::readInputGMSH(string& fname)
{
    PetscPrintf(MPI_COMM_WORLD, " BaseFEM::readInputGMSH \n");
    
    //mesh = make_unique<myMesh>();
    mesh = make_shared<myMesh>();

    mesh->readInputGMSH(fname);
    
    ndof = 2;

    PetscPrintf(MPI_COMM_WORLD, " BaseFEM::readInputGMSH \n");

    return 0;
}






int femCahnHilliard::readConfiguration(string& fname)
{
    cout << " femCahnHilliard::readConfiguration " << endl;

    ifstream  infile(fname);

    if(infile.fail())
    {
       cout << " Could not open input file " << endl;
       exit(-1);
    }

    ndim = mesh->getNDIM();
    ndof = 2;

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
                readBodyForce(infile, line);
            }
            else if(line.compare(string("Boundary Conditions")) == 0)
            {
                cout << "Boundary Conditions" << endl;
                // read the number of prescribed BCs
                readBoundaryConditions(infile, line);
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
                readInitialConditions(infile, line);
            }
            else if(line.compare(string("Patch Output")) == 0)
            {
                cout << "Patch Output" << endl;
                readOutputDetailsPatch(infile, line);
            }
            else
            {
                cout << "key =  " <<  line << endl;
                throw runtime_error("Key not found in femCahnHilliard::readConfiguration ...");
                //return -1;
            }
      }
    }

    infile.close();

    cout << " Configuration file is successfully read " << endl;

    return 0;
}







int  femCahnHilliard::readBodyForce(ifstream& infile, string& line)
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

            if(stringlist[0] == "value")
            {
                bodyForce[0] = stod(stringlist[1]);
                bodyForce[1] = stod(stringlist[2]);
                bodyForce[2] = stod(stringlist[3]);
            }
            else if(stringlist[0] == "timefunction")
            {
                bodyForceTimeFunction = stoi(stringlist[1])-1;
            }
            else
            {
                throw runtime_error("Option not available in femCahnHilliard::readBodyForce");
            }
        }
    }


    return 0;
}




int  femCahnHilliard::readBoundaryConditions(ifstream& infile, string& line)
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

            for(auto& bpt : mesh->BoundaryPatches)
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
                    boost::algorithm::split(stringlist, line, boost::is_any_of(" "), boost::token_compress_on);
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
                        throw runtime_error("Option not available in femCahnHilliard::readBoundaryConditions");
                    }
                }//if
            } //while
        }//if
    }//while

    return 0;
}





int  femCahnHilliard::readSolverDetails(ifstream& infile, string& line)
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
                throw runtime_error("Option not available in femCahnHilliard::readBodyForce");
            }
        }
    }


    return 0;
}



int  femCahnHilliard::readTimeFunctions(ifstream& infile, string& line)
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





int femCahnHilliard::readInitialConditions(ifstream& infile, string& line)
{
    vector<string>  stringlist;
    string label;

    // read {
    getline(infile,line);    boost::trim(line);
    
    initialConitions.resize(ndof);
    //for(vector<myMathFunction>::iterator itr = initialConitions.begin(); itr!=initialConitions.end(); ++itr)
      //itr->initialise("0.0");

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
                throw runtime_error("Specified DOF higher than #DOF in myMesh::readInitialConditions ...");
            }

            initialConitions[dof] = stringlist[1];
            //initialConitions[dof].initialise(stringlist[1]);
        }
    }

    return 0;
}









int femCahnHilliard::readFluidProperties(ifstream& infile, string& line)
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

            if(stringlist[0] == "density")
            {
                fluidProperties[0] = stod(stringlist[1]);
            }
            else if(stringlist[0] == "viscosity")
            {
                fluidProperties[1] = stod(stringlist[1]);
            }
            else
            {
                throw runtime_error("Option not available in femCahnHilliard::readFluidProperties");
            }
        }
    }


    return 0;
}







int  femCahnHilliard::readOutputDetailsPatch(ifstream& infile, string& line)
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
            for(auto& bpt : mesh->BoundaryPatches)
            {
               if( stringlist[0] == bpt->getLabel() )
               {
                   bpt->setOutputFlag();
                   break;
               }
               ind++;
            }
            
            cout << ind << '\t' << mesh->numBoundaryPatches << endl;

            if(ind > mesh->numBoundaryPatches)
            {
                throw runtime_error("Patch not available in femCahnHilliard::readOutputDetailsPatch");
            }
        }
    }

    return 0;
}




int femCahnHilliard::prepareInputData()
{
    PetscPrintf(MPI_COMM_WORLD, "\n\n  femCahnHilliard::prepareInputData()  .... STARTED ...\n\n");

    int ii, jj, kk, ee, nn, ind, n1, n2, dof;

    if(this_mpi_proc == 0)
    {
        cout << " Boundary patches ...... " << endl;
        cout << " -------------------------- " << endl;
        cout << endl;
        // add boundary Conditions
        //for(auto& bpt : mesh->BoundaryPatches)
          //bpt->printData();

        cout << " Boundary conditions ...... " << endl;
        cout << " -------------------------- " << endl;
        cout << endl;

        for(auto&  bc : boundaryConitions)
        {
            cout << " Patch                 = " << bc->label << endl;
            cout << " BCType                = " << bc->BCType << endl;
            cout << " dof_specified_string  = " << bc->dof_specified_string << endl;
            cout << " dof_specified_int     = " << bc->dof_specified_int << endl;
            cout << " expression            = " << bc->expression << endl;
            cout << " timeFunctionNum       = " << bc->timeFunctionNum << endl;
            cout << endl;
            cout << endl;
        }
    }


    if(!mesh->MESH_PREPARED)
      mesh->prepareInputData();
    
    elem_proc_id     = mesh->elem_proc_id;
    node_proc_id     = mesh->node_proc_id;
    node_map_get_old = mesh->node_map_get_old;
    node_map_get_new = mesh->node_map_get_new;

    // ==================================================
    //
    // Check the  consistency of input data
    //
    // ==================================================

    //checkInputData();
    ///////////////////////////////////////////////////////////////////
    //
    // set SolnData details
    //
    ///////////////////////////////////////////////////////////////////

    ndim = mesh->getNDIM();
    ind = mesh->getnNodeGlobal()*ndof;

    soln.resize(ind);
    soln.setZero();

    solnInit  = soln;
    solnPrev  = soln;
    solnPrev2 = soln;
    solnPrev3 = soln;
    solnPrev4 = soln;
    solnCur   = soln;

    solnDot     = soln;
    solnDotPrev = soln;
    solnDotCur  = soln;
    solnApplied = soln;
    
    
    PetscPrintf(MPI_COMM_WORLD, "\n\n  femCahnHilliard::prepareInputData()  .... FINISHED ...\n\n");

    return 0;
}




int femCahnHilliard::setSpecifiedDOFs(vector<vector<bool> >& NodeType)
{
    int  i, nn, dof;
    vector<int> nodeNums;

    dofs_specified.clear();

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
        //nodeNums = mesh->BoundaryPatches[index]->nodeNums;

        if(bc->BCType == "specified")
        {
            dof = bc->dof_specified_int;
            
            for(i=0; i<nodeNums.size(); ++i)
            {
                nn = nodeNums[i];

                NodeType[nn][dof] = true;
                dofs_specified.push_back(nn*ndof+dof);
            }
        }
        else if(bc->BCType == "wall")
        {
            for(i=0; i<nodeNums.size(); ++i)
            {
              for(dof=0; dof<mesh->getNDIM(); ++dof)
              {
                nn = nodeNums[i];

                NodeType[nn][dof] = true;
                dofs_specified.push_back(nn*ndof+dof);
              }
            }
        }
        else
        {
            throw runtime_error("Patch type not available in femCahnHilliard::setBoundaryConditions");
        }
    }
    findUnique(dofs_specified);
    //printVector(dofs_specified);

    return 0;
}




int femCahnHilliard::setBoundaryConditions()
{
    int  i, nn, dof, index=-1;
    double xc, yc, zc, value, timeFactor;
    vector<int> nodeNums;

    solnApplied.setZero();

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
        //cout << " Patch                 = " << bc->bpt->getLabel() << endl;

        // nodeNums contains new node numbers which are used for the solution
        // nodeCoordsOrig contains coordinates of nodes with numbers before domain decompositions

        nodeNums = bc->bpt->nodeNums;

        if(bc->BCType == "specified")
        {
            myMathFunction  mathfun;
            mathfun.initialise(bc->expression);

            if(bc->timeFunctionNum == -1)
              timeFactor = 1.0;
            else
              timeFactor = timeFunction[bc->timeFunctionNum]->getFactor();

            PetscPrintf(MPI_COMM_WORLD, " bc->label = %10s ...  bc->timeFunctionNum = %d ... time Factor = %12.6f \n", bc->label.c_str(), bc->timeFunctionNum,  timeFactor);

            dof = bc->dof_specified_int;
            
            for(i=0; i<nodeNums.size(); ++i)
            {
                nn = node_map_get_old[nodeNums[i]];

                xc = mesh->nodeCoordsOrig[nn][0];
                yc = mesh->nodeCoordsOrig[nn][1];
                zc = mesh->nodeCoordsOrig[nn][2];

                value = mathfun.getValue(xc, yc, zc) * timeFactor;

                //cout << xc << '\t' << yc << '\t' << zc << '\t' << timeFactor << '\t' << value << endl;

                solnApplied[nodeNums[i]*ndof+dof] = value;
            }
        }
        else if(bc->BCType == "wall")
        {
            for(i=0; i<nodeNums.size(); ++i)
            {
              for(dof=0; dof<mesh->getNDIM(); ++dof)
              {
                solnApplied[nodeNums[i]*ndof+dof] = 0.0;
              }
            }
        }
        else
        {
            throw runtime_error("Patch type not available in femCahnHilliard::setBoundaryConditions");
        }
    }

    solnApplied -= soln;

    //printVector(solnApplied);

    return 0;
}




int femCahnHilliard::setInitialConditions()
{
    PetscPrintf(MPI_COMM_WORLD, "femCahnHilliard::setInitialConditions() ... \n");

    double  x=0.0, y=0.0, z=0.0, eps = 0.02;
    string  expr;

    cout << solnPrev.rows() << endl;

    solnPrev.setZero();
    
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

        for(int ii=0; ii<mesh->getnNodeGlobal(); ++ii)
        {
          x = mesh->nodeCoordsOrig[ii][0];
          y = mesh->nodeCoordsOrig[ii][1];
          z = mesh->nodeCoordsOrig[ii][2];

          solnPrev(node_map_get_new[ii]*ndof + dof) = expression.value();
        }
      }
    }
    PetscPrintf(MPI_COMM_WORLD, "femCahnHilliard::setInitialConditions() ... \n");

    // add boundary Conditions
    //setBoundaryConditions();
    //solnPrev += solnApplied;

    soln = solnPrev;

    PetscPrintf(MPI_COMM_WORLD, "femCahnHilliard::setInitialConditions() ... \n");

    return 0;
}






int femCahnHilliard::setTimeParam()
{
  //setTimeParam();

  return 0;
}




int femCahnHilliard::timeUpdate()
{
  myTime.update();

  // set parameters for the time integration scheme.
  // need to be done every time step to account for adaptive time stepping
  SetTimeParametersFluid(timeIntegrationScheme, spectralRadius, myTime.dt, td);

  // update time functions
  for(auto& tmf : timeFunction)
    tmf->update();


  setBoundaryConditions();
  //if(this_mpi_proc == 0) printVector(SolnData.solnApplied);


  return 0;
}



int femCahnHilliard::updateIterStep()
{
  solnDot    = td[9]*soln + td[10]*solnPrev + td[11]*solnPrev2 + td[12]*solnPrev3 + td[13]*solnPrev4 + td[15]*solnDotPrev ;

  solnCur    = td[2]*soln    + (1.0-td[2])*solnPrev;
  solnDotCur = td[1]*solnDot + (1.0-td[1])*solnDotPrev;
  
  return 0;
}




int  femCahnHilliard::reset()
{
  solnPrev  = solnPrev2;
  solnPrev2 = solnPrev3;
  solnPrev3 = solnPrev4;

  solnDot   = solnDotPrev;

  return 0;
}



int  femCahnHilliard::saveSolution()
{
  solnPrev4  = solnPrev3;
  solnPrev3  = solnPrev2;
  solnPrev2  = solnPrev;
  solnPrev   = soln;

  solnDotPrev  = solnDot;

  return 0;
}





int femCahnHilliard::getPhi(VectorXd& phi)
{
    assert(phi.size() == mesh->getnNodeGlobal());

    for(int i=0; i<mesh->getnNodeGlobal(); i++)
    {
        phi[i] = soln[i*2];
    }
    
    return 0;
}





int femCahnHilliard::getEta(VectorXd& eta)
{
    assert(eta.size() == mesh->getnNodeGlobal());

    for(int i=0; i<mesh->getnNodeGlobal(); i++)
    {
        eta[i] = soln[i*2+1];
    }
    
    return 0;
}





bool  femCahnHilliard::converged()
{
  if(rhsNorm < conv_tol)
    return true;

  return false;
}




int femCahnHilliard::writeOutputData()
{
    PetscPrintf(MPI_COMM_WORLD, "femCahnHilliard::writeOutputData() ... \n");
    
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

    return 0;
}



int  femCahnHilliard::writeResult(string&  fname)
{
/*
    ofstream fout(fname);

    if(fout.fail())
    {
       cout << " Could not open the file in 'femCahnHilliard::writeResult' ... " << endl;
       exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(8);

    //fout << nNode_global <<endl;

    int  ii, nn, ind;
    if(ndim == 2)
    {
      for(int ii=0; ii<nNode_global; ii++)
      {
        ind = node_map_get_new[ii]*ndof;
        fout << soln[ind] << '\t' << soln[ind+1] << '\t' << soln[ind+2] << endl;
      }
    }
    else
    {
      for(int ii=0; ii<nNode_global; ii++)
      {
        ind = node_map_get_new[ii]*ndof;
        fout << soln[ind] << '\t' << soln[ind+1] << '\t' << soln[ind+2] << '\t' << soln[ind+3] << endl;
      }
    }

    fout.close();
*/
    return  0;
}






int  femCahnHilliard::checkResult(string&  fname)
{
/*
    std::ifstream  datfile(fname);

    if(datfile.fail())
    {
        cout << " Could not open the file in 'femCahnHilliard::checkResult' ... " << endl;
        exit(1);
    }

    string  line, stringVal, stringVec[10];
    double  tempDbl;

    int  nnodes, ii, ind;
    // read number of nodes
    datfile  >>  nnodes;
    if(nnodes != nNode_global)
      return  -1;

    double  val[ndof];
    bool  flag = true;
    if(ndim == 2)
    {
      for(ii=0; ii<nNode_global; ++ii)
      {
        datfile >> val[0] >> val[1] >> val[2];
        //cout << val[0] << '\t' << val[1] << '\t' << val[2] << endl;

        ind = node_map_get_new[ii]*ndof;
        val[0] -= soln[ind];
        val[1] -= soln[ind+1];
        val[2] -= soln[ind+2];

        if( sqrt(val[0]*val[0] + val[1]*val[1] + val[2]*val[2]) > 1.0e-4)
        {
          flag = false;
          break;
        }
      }
    }
    else
    {
      for(ii=0; ii<nNode_global; ++ii)
      {
        datfile >> val[0] >> val[1] >> val[2] >> val[3];
      }
    }

    datfile.close();
    if(!flag)
      return -2;
*/

    return  0;
}




int  femCahnHilliard::readResult(string&  fname)
{
    ofstream fout(fname);

    if(fout.fail())
    {
       cout << " Could not open the file in 'femCahnHilliard::readResult' ... " << endl;
       exit(1);
    }

    fout.setf(ios::fixed);
    fout.setf(ios::showpoint);
    fout.precision(8);


  return  0;
}








int femCahnHilliard::printComputerTimes()
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





int  femCahnHilliard::postProcess()
{
    if(this_mpi_proc != 0)
      return 0;

    PetscPrintf(MPI_COMM_WORLD, " \n femCahnHilliard::postProcess \n");

    //
    // setup and write vtk data
    //
    //////////////////////////////////////////////


    vtkSmartPointer<vtkUnstructuredGrid>     uGridVTK      =  vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints>               pointsVTK     =  vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkVertex>               vertexVTK     =  vtkSmartPointer<vtkVertex>::New();

    vtkSmartPointer<vtkTriangle>             triaVTK       =  vtkSmartPointer<vtkTriangle>::New();
    vtkSmartPointer<vtkQuad>                 quadVTK       =  vtkSmartPointer<vtkQuad>::New();
    vtkSmartPointer<vtkTetra>                tetraVTK      =  vtkSmartPointer<vtkTetra>::New();
    vtkSmartPointer<vtkHexahedron>           hexaVTK       =  vtkSmartPointer<vtkHexahedron>::New();

    vtkSmartPointer<vtkQuadraticTriangle>    tria6VTK      =  vtkSmartPointer<vtkQuadraticTriangle>::New();
    vtkSmartPointer<vtkBiQuadraticQuad>      quad9VTK      =  vtkSmartPointer<vtkBiQuadraticQuad>::New();
    vtkSmartPointer<vtkQuadraticTetra>       tetra10VTK    =  vtkSmartPointer<vtkQuadraticTetra>::New();
    vtkSmartPointer<vtkQuadraticHexahedron>  hexa27VTK     =  vtkSmartPointer<vtkQuadraticHexahedron>::New();

    vtkSmartPointer<vtkFloatArray>           vecVTK        =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK        =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           scaVTK2       =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           procIdVTKnode =  vtkSmartPointer<vtkFloatArray>::New();
    vtkSmartPointer<vtkFloatArray>           procIdVTKcell =  vtkSmartPointer<vtkFloatArray>::New();

    vtkSmartPointer<vtkXMLUnstructuredGridWriter>  writerUGridVTK =  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();


    int  ee, ii, jj, kk, nn, n1, n2, npElem;
    double  val, vec[3]={0.0,0.0,0.0};
    vector<int> nodeNums;

    scaVTK->SetName("phi");
    scaVTK->SetNumberOfTuples(mesh->nNode_global);

    scaVTK2->SetName("eta");
    scaVTK2->SetNumberOfTuples(mesh->nNode_global);

    procIdVTKcell->SetName("procId");
    procIdVTKcell->SetNumberOfTuples(mesh->nElem_global);

    procIdVTKnode->SetName("procId");
    procIdVTKnode->SetNumberOfTuples(mesh->nNode_global);

    vtkIdType pt[10];

    if(mesh->ndim == 2)
    {
      for(ii=0; ii<mesh->nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(mesh->nodeCoordsOrig[ii][0], mesh->nodeCoordsOrig[ii][1], 0.0);

        procIdVTKnode->InsertTuple1(ii, node_proc_id[ii]);

        nn = node_map_get_new[ii];
        kk = nn*ndof;

        scaVTK->SetTuple1(ii, soln[kk]);
        scaVTK2->SetTuple1(ii, soln[kk+1]);
      }

      for(ee=0; ee<mesh->nElem_global; ee++)
      {
        procIdVTKcell->SetTuple1(ee, elem_proc_id[ee]);

        nodeNums = mesh->elemConn[ee];
        npElem   = nodeNums.size();

        if(npElem == 3)
        {
          for(ii=0; ii<npElem; ii++)
            triaVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
        }
        else if(npElem == 6)
        {
          for(ii=0; ii<npElem; ii++)
            tria6VTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(tria6VTK->GetCellType(), tria6VTK->GetPointIds());
        }
        else if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            quadVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
        }
        else if(npElem == 9)
        {
          for(ii=0; ii<npElem; ii++)
            quad9VTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(quad9VTK->GetCellType(), quad9VTK->GetPointIds());
        }
      }
    }
    else // (ndim == 3)
    {
      for(ii=0; ii<mesh->nNode_global; ii++)
      {
        pt[0] = pointsVTK->InsertNextPoint(mesh->nodeCoordsOrig[ii][0], mesh->nodeCoordsOrig[ii][1], mesh->nodeCoordsOrig[ii][2]);

        procIdVTKnode->InsertTuple1(ii, node_proc_id[ii]);

        nn = node_map_get_new[ii];
        kk = nn*ndof;

        vec[0] = soln[kk];
        vec[1] = soln[kk+1];
        vec[2] = soln[kk+2];
        val    = soln[kk+3];

        vecVTK->InsertTuple(ii, vec);
        scaVTK->SetTuple1(ii, val);
      }

      for(ee=0; ee<mesh->nElem_global; ee++)
      {
        procIdVTKcell->SetTuple1(ee, elem_proc_id[ee]);

        nodeNums = mesh->elemConn[ee];
        npElem   = nodeNums.size();

        if(npElem == 4)
        {
          for(ii=0; ii<npElem; ii++)
            tetraVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(tetraVTK->GetCellType(), tetraVTK->GetPointIds());
        }
        else if(npElem == 8)
        {
          for(ii=0; ii<npElem; ii++)
            hexaVTK->GetPointIds()->SetId(ii, node_map_get_old[nodeNums[ii]] );

          uGridVTK->InsertNextCell(hexaVTK->GetCellType(), hexaVTK->GetPointIds());
        }
      }
    }

    uGridVTK->SetPoints(pointsVTK);
    uGridVTK->GetPointData()->SetScalars(scaVTK);
    uGridVTK->GetPointData()->AddArray(scaVTK2);
    uGridVTK->GetPointData()->AddArray(procIdVTKnode);

    uGridVTK->GetCellData()->SetScalars(procIdVTKcell);

    //Write the file.

    char VTKfilename[200];

    sprintf(VTKfilename,"%s%s%06d%s", infilename.c_str(), "-",fileCount,".vtu");

    writerUGridVTK->SetFileName(VTKfilename);

    writerUGridVTK->SetInputData(uGridVTK);
    writerUGridVTK->Write();

    fileCount++;

    PetscPrintf(MPI_COMM_WORLD, " \n femCahnHilliard::postProcess \n");

    return 0;
}









#include "BoundaryPatch.h"
#include "util.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "mymapkeys.h"

extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime               myTime;


BoundaryPatch::BoundaryPatch(string& label_, int tag_, int ndim_)
{
    label = label_;
    tag   = tag_;
    ndim  = ndim_;
    //ndof  = ndof_;

    outputflag = false;
}



BoundaryPatch::~BoundaryPatch()
{

}




void BoundaryPatch::addElement(vector<int> & elnodeNums)
{
    elemConn.push_back(elnodeNums);

    for(int i=0; i<elnodeNums.size(); ++i)
    nodeNums.push_back(elnodeNums[i]);

    return;
}



void BoundaryPatch::processData()
{
    nElem = elemConn.size();

    findUnique(nodeNums);
    nNode = nodeNums.size();

    return;
}






void  BoundaryPatch::setOutputFlag()
{
    outputflag = true;

    forcedata.open(string("forces-"+label+".dat"), ios::out | ios::trunc );

    if(forcedata.fail())
    {
        cout << " Could not open the output file for writing forces " << endl;
        exit(1);
    }

    forcedata.setf(ios::fixed);
    forcedata.setf(ios::showpoint);
    forcedata.precision(14);

    return;
}





int  BoundaryPatch::readData(ifstream& infile, string& line)
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

            if(stringlist[0] == "type")
            {
                BCType.push_back(stringlist[1]);
            }
            else if(stringlist[0] == "dof")
            {
                cout << stringlist[0] << "\t" << stringlist[1] << endl;

                dof_specified.push_back(getDOFfromString(stringlist[1]));

                //set the default time function
                timeFunctionNums.push_back(-1);
            }
            else if(stringlist[0] == "value")
            {
                //myMathFunction  mathfun;
                //mathfun.initialise(stringlist[1]);

                expressions.push_back(stringlist[1]);
                //valueexpressions.push_back(myMathFunction(stringlist[1]));
            }
            else if(stringlist[0] == "timefunction")
            {
                //timeFunctionNums.push_back(stoi(stringlist[1]));
                // update the time function
                timeFunctionNums[timeFunctionNums.size()-1] = stoi(stringlist[1])-1;
            }
            else if(stringlist[0] == "output")
            {
                outputflag = (stringlist[1] == "on" || stringlist[1] == "yes");
            }
            else
            {
                throw runtime_error("Option not available in BoundaryPatch::readData");
            }
        }
    }

    return 0;
}





void BoundaryPatch::updateNodeNumbers(vector<int>& node_map_get_new)
{
    int size = BCType.size(), el, i, nn, dof;

    for(i=0; i<nNode; ++i)
    {
        nn = nodeNums[i];
        nodeNums[i] = node_map_get_new[nn];
    }
    
    for(el=0; el<nElem; ++el)
    {
        for(i=0; i<elemConn[el].size(); ++i)
        {
            nn = elemConn[el][i];

            elemConn[el][i] = node_map_get_new[nn];
        }
    }

    return;
}





/*
void BoundaryPatch::setBoundaryConditions(int ndof, vector<int>& node_map_get_old, vector<myPoint>& node_coords, VectorXd& solnApplied)
{
    int size = BCType.size(), el, i, nn, dof;
    double xc, yc, zc, value, timeFactor;

    //cout << " BCType ... size = " << size  << endl;
    //cout << " nNode           = " << nNode << endl;
    //cout << " ndim            = " << ndim  << endl;
    //cout << " ndof            = " << ndof  << endl;

    for(el=0; el<size; ++el)
    {
        //cout << " BCType[el]    ... = " << BCType[el] << endl;

        if(BCType[el] == "specified")
        {
            myMathFunction  mathfun;
            mathfun.initialise(expressions[el]);

            timeFactor = timeFunction[timeFunctionNums[el]]->getFactor();

            //cout << el << '\t' << timeFunctionNums[el] << '\t' << timeFactor << endl;

            dof = dof_specified[el];
            
            // nodeNums contains new node numbers which are used for the solution
            // node_coords contains coordinates of nodes with numbers before domain decompositions

            for(i=0; i<nNode; ++i)
            {
                nn = node_map_get_old[nodeNums[i]];

                xc = node_coords[nn][0];
                yc = node_coords[nn][1];
                zc = node_coords[nn][2];

                value = mathfun.getValue(xc, yc, zc) * timeFactor;

                //cout << xc << '\t' << yc << '\t' << zc << '\t' << timeFactor << '\t' << value << endl;

                //DirichletBCs.push_back(make_tuple(nn, dof, value));
                solnApplied[nodeNums[i]*ndof+dof] = value;
            }
        }
        else if(BCType[el] == "wall")
        {
            for(i=0; i<nNode; ++i)
            {
              for(dof=0; dof<ndim; ++dof)
              {
                //DirichletBCs.push_back(make_tuple(nodeNums[i], dof, 0.0));
                solnApplied[nodeNums[i]*ndof+dof] = 0.0;
              }
            }
        }
        else if(BCType[el] == "traction")
        {
        }
        else
        {
            throw runtime_error("Patch type not available in BoundaryPatch::processBoundaryConditions");
        }
    }

    return;
}




void BoundaryPatch::setSpecifiedDOFs(int ndof, vector<int>& node_map_get_new, vector<vector<bool> >& NodeType, vector<int>& dofs_specified)
{
    int size = BCType.size(), el, i, nn, dof;

    for(el=0; el<size; ++el)
    {
        //cout << " BCType[el]    ... = " << BCType[el] << endl;
        if(BCType[el] == "specified")
        {
            dof = dof_specified[el];

            for(i=0; i<nNode; ++i)
            {
                //nn = node_map_get_new[nodeNums[i]];
                nn = nodeNums[i];

                NodeType[nn][dof] = true;
                dofs_specified.push_back(nn*ndof+dof);
            }
        }
        else if(BCType[el] == "wall")
        {
            for(i=0; i<nNode; ++i)
            {
              for(dof=0; dof<ndim; ++dof)
              {
                //nn = node_map_get_new[nodeNums[i]];
                nn = nodeNums[i];

                NodeType[nn][dof] = true;
                dofs_specified.push_back(nn*ndof+dof);
              }
            }
        }
        else if(BCType[el] == "traction")
        {
        }
        else
        {
            throw runtime_error("Patch type not available in BoundaryPatch::processBoundaryConditions");
        }
    }

    return;
}
*/



/*
void BoundaryPatch::processBoundaryConditions(vector<myPoint>& node_coords)
{
    int size = BCType.size(), el, i, nn, dof;
    
    cout << " BCType ... size = " << size  << endl;
    cout << " nNode           = " << nNode << endl;
    cout << " ndof            = " << ndof  << endl;

    for(el=0; el<size; ++el)
    {
        if(BCType[el] == "specified")
        {
            cout << " BCType[el]    ... = " << BCType[el] << endl;
            cout << " valuetype[el] ... = " << valuetype[el] << endl;

            dof = dof_specified[el];

            if(valuetype[el] == "constant")
            {

              for(i=0; i<nNode; ++i)
              {
                DirichletBCs.push_back(make_tuple(nodeNums[i], dof, fielddata[el][0]));
              }

            }
            else if(valuetype[el] == "parabolicX")
            {

              double Umagn = fielddata[el][0], x0 = fielddata[el][1], x1 = fielddata[el][2];
              double factor = Umagn * (6.0)/(x1-x0)/(x1-x0), value;

              for(i=0; i<nNode; ++i)
              {
                nn = nodeNums[i];

                value = factor * (node_coords[nn][0]-x0) * (x1-node_coords[nn][0]);
                DirichletBCs.push_back(make_tuple(nn, dof, value));
              }

            }
            else if(valuetype[el] == "parabolicY")
            {

              double Umagn = fielddata[el][0], y0 = fielddata[el][1], y1 = fielddata[el][2], y;
              double factor = Umagn * (6.0)/(y1-y0)/(y1-y0), value;

              cout << " valuetype[el] ... = " << Umagn << '\t' << y0 << '\t' << y1 << '\t' << factor << endl;

              for(i=0; i<nNode; ++i)
              {
                nn = nodeNums[i];
                y = node_coords[nn][1];
                value = factor * (y-y0) * (y1-y);

                //cout << " value ... = " << y << '\t' << value << endl;

                DirichletBCs.push_back(make_tuple(nn, dof, value));
              }

            }
            else if(valuetype[el] == "parabolicZ")
            {

              double Umagn = fielddata[el][0], z0 = fielddata[el][1], z1 = fielddata[el][2];
              double factor = Umagn * (6.0)/(z1-z0)/(z1-z0), value;

              for(i=0; i<nNode; ++i)
              {
                nn = nodeNums[i];

                value = factor * (node_coords[nn][2]-z0) * (z1-node_coords[nn][2]);
                DirichletBCs.push_back(make_tuple(nn, dof, value));
              }
           }
        }
        else if(BCType[el] == "wall")
        {
            for(i=0; i<nNode; ++i)
            {
              for(dof=0; dof<ndof; ++dof)
              {
                DirichletBCs.push_back(make_tuple(nodeNums[i], dof, 0.0));
              }
            }
        }
        else
        {
            throw runtime_error("Patch type not available in BoundaryPatch::processBoundaryConditions");
        }
    }

    return;
}
*/




void BoundaryPatch::printData()
{
    cout << "-----------------------------------" << endl;
    cout << "Patch tag   ... = " << tag  << endl;
    cout << "Patch label ... = " << label << endl;
    cout << "nNode           = " << nNode << endl;
    cout << "nElem           = " << nElem << endl;
    cout << endl;

    cout << " nElem = " << nElem << endl;
    cout << "Elements ..."  << endl;
    for(int ii=0; ii<nElem; ++ii)
      printVector(elemConn[ii]);

    cout << endl;
    cout << endl;

    cout << "Number of nodes = " << nNode << endl;
    cout << "Node numbers "  << endl;
      printVector(nodeNums);


    cout << "Boundary conditions" << endl;
    cout << endl;

    for(int el=0; el<BCType.size(); ++el)
    {
        cout << " BCType[el]        ... = " << BCType[el] << endl;
        cout << " dof_specified[el] ... = " << dof_specified[el] << endl;
        cout << " expressions[el]   ... = " << expressions[el] << endl;
    }

    cout << "-----------------------------------" << endl;
    cout << endl;
    cout << endl;

    return;
}




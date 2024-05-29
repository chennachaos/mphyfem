#include "BoundaryCondition.h"
#include "util.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "mymapkeys.h"

extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime               myTime;


BoundaryCondition::BoundaryCondition()
{
    dof_specified_int = 0;
    timeFunctionNum   = -1;
}



BoundaryCondition::~BoundaryCondition()
{

}

/*
void BoundaryCondition::setBoundaryConditions(int ndof, vector<int>& node_map_get_old, vector<myPoint>& node_coords, VectorXd& solnApplied)
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
            throw runtime_error("Patch type not available in BoundaryCondition::processBoundaryConditions");
        }
    }

    return;
}
*/


/*
void BoundaryCondition::setSpecifiedDOFs(int ndof, vector<int>& nodeNums, vector<vector<bool> >& NodeType, vector<int>& dofs_specified)
{
    int i, nn, dof;

    for(el=0; el<size; ++el)
    {
        //cout << " BCType[el]    ... = " << BCType[el] << endl;
        if(BCType[el] == "specified")
        {
            dof = dof_specified[el];

            for(i=0; i<nNode; ++i)
            {
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
            throw runtime_error("Patch type not available in BoundaryCondition::processBoundaryConditions");
        }
    }

    return;
}

*/



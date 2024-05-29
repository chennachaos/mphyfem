#ifndef incl_BoundaryCondition_h
#define incl_BoundaryCondition_h

#include <string>
#include <memory>
#include "BoundaryPatch.h"

using std::string;


class BoundaryCondition
{
public:
    string           label;
    string           BCType;               // specified or wall or traction
    string           dof_specified_string; // specified DOF
    int              dof_specified_int;    // specified DOF
    string           expression;          // expressions for specified boundary conditions
    int              timeFunctionNum;     // time function numbers
    //shared_ptr<BoundaryPatch>   bpt;
    BoundaryPatch*   bpt;


    BoundaryCondition();

    ~BoundaryCondition();

    int  getTimeFunctionNumber()
    { return timeFunctionNum; }

    //void setSpecifiedDOFs(int ndof, vector<int>& nodeNums, vector<vector<bool> >& NodeTypeOld, vector<int>& dofs_specified);
};

#endif


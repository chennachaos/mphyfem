#ifndef incl_TractionCondition_h
#define incl_TractionCondition_h

#include <string>
#include <memory>
#include "BoundaryPatch.h"

using std::string;


class TractionCondition
{
public:
    string           label;
    string           BCType;               // specified or wall or traction
    //string           expression;          // expressions for specified boundary conditions
    vector<double>   specValues;
    int              timeFunctionNum;     // time function numbers
    //shared_ptr<BoundaryPatch>   bpt;
    BoundaryPatch*   bpt;


    TractionCondition();

    ~TractionCondition();
};



#endif


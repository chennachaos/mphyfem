#ifndef incl_NodalForces_h
#define incl_NodalForces_h

#include <string>
#include <memory>

using std::string;


struct NodalForces
{
    int   timeFunctionNum=-1;     // time function number
    int   count=0;
    vector<vector<double> >  data;
};



#endif


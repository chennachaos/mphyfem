
#ifndef incl_Global_h
#define incl_Global_h


#include "MyTime.h"
#include "TimeFunction.h"
#include <vector>
#include <memory>

MyTime myTime;

vector<unique_ptr<TimeFunction> > timeFunctions;

bool debug = true;



#endif




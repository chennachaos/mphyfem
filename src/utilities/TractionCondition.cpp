#include "TractionCondition.h"
#include "util.h"
#include "MyTime.h"
#include "TimeFunction.h"
#include <unordered_map>
#include <boost/algorithm/string.hpp>
#include "mymapkeys.h"

extern   std::vector<unique_ptr<TimeFunction> > timeFunction;
extern MyTime               myTime;


TractionCondition::TractionCondition()
{
    timeFunctionNum  = -1;
    specValues.resize(3);
}



TractionCondition::~TractionCondition()
{

}

#include "InternalVariables.h"
#include "util.h"


InternalVariables::InternalVariables()
{
}



InternalVariables::~InternalVariables()
{

}



int  InternalVariables::initialise(int rows, int cols)
{
    var.resize(rows, cols); var.setZero();

    varPrev    = var;
    varDot     = var;
    varDotPrev = var;

    return 0;
}



int InternalVariables::update(VectorXd&  stre, MatrixXd&  Cmat)
{
    return 0;
}




int InternalVariables::reset()
{
    var    = varPrev;
    varDot = varDotPrev;

    return 0;
}




int InternalVariables::saveSolution()
{
    varPrev    = var;
    varDotPrev = varDot;

    return 0;
}





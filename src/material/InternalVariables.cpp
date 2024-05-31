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

    for(int jj=0; jj<cols; jj++)
    {
        var(0,jj) = 1.0;
        var(4,jj) = 1.0;
        var(8,jj) = 1.0;
    }

    varPrev    = var;
    varDot     = var*0.0;
    varDotPrev = var*0.0;

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





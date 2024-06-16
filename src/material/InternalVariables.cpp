#include "InternalVariables.h"
#include "util.h"


InternalVariables::InternalVariables()
{
}



InternalVariables::~InternalVariables()
{

}



int  InternalVariables::initialise(int nGP, int nivGP)
{
    int  size = nGP*nivGP;
    var.resize(size);
    varPrev.resize(size);
    varDot.resize(size);
    varDotPrev.resize(size);

    Matrix3d Id = Matrix3d::Identity();

    for(int ii=0; ii<var.size(); ii++)
    {
        var[ii] = Id;
        varPrev[ii] = Id;

        varDot[ii].setZero();
        varDotPrev[ii].setZero();
    }

    return 0;
}



int InternalVariables::update(VectorXd&  stre, MatrixXd&  Cmat)
{
    return 0;
}




int InternalVariables::reset()
{
    for(int ii=0; ii<var.size(); ii++)
    {
        var[ii]    = varPrev[ii];
        varDot[ii] = varDotPrev[ii];
    }

    return 0;
}




int InternalVariables::saveSolution()
{
    for(int ii=0; ii<var.size(); ii++)
    {
        varPrev[ii]    = var[ii];
        varDotPrev[ii] = varDot[ii];
    }

    return 0;
}





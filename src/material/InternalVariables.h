#ifndef incl_InternalVariables_h
#define incl_InternalVariables_h

#include <vector>
#include "headersEigen.h"

using std::vector;
using std::cout;
using std::endl;

using Eigen::VectorXd;
using Eigen::MatrixXd;



class InternalVariables
{
  public:

    //member variables
    vector<double>  matData;

    MatrixXd  var, varPrev, varDot, varDotPrev;

    //member functions

    InternalVariables();

    ~InternalVariables();

    int  initialise(int rows, int cols);

    int  update(VectorXd&  stre, MatrixXd&  Cmat);

    int  saveSolution();

    int  reset();
};




#endif

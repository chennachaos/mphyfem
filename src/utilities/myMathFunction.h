#ifndef incl_myMathFunction_h
#define incl_myMathFunction_h


//https://github.com/ArashPartow/exprtk/tree/master
//https://www.partow.net/programming/exprtk/index.html
#include "exprtk.hpp"

#include <vector>
#include <iostream>

using namespace std;



class myMathFunction
{
  public:

    //member variables
    double  x, y, z;

    exprtk::expression<double>   expression;

    string          expr;

    //member functions

    myMathFunction();

    ~myMathFunction();

    myMathFunction(const string& expr_);

    void initialise(const string& expr_);

    double getValue(double x_=0, double y_=0, double z_=0);

    string  getExpression()
    { return expr; }
};

#endif


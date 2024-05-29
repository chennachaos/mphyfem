
#include "myExpression.h"
#include <iostream>
#include "util.h"


using namespace std;


void myExpression::set(double dt_, double dtMin_, double dtMax_)
{
    dt    = dt_;
    dtMin = dtMin_;
    dtMax = dtMax_;

    stack.push_back(dt);

    dtPrev2 = dt;
    dtPrev  = dt;

    return;
}



void myExpression::update()
{
    prev2 = prev;
    prev  = cur;
    cur  += dt;

    return;
}



void myExpression::cut()
{
  cur -= dt;

  dtPrev2 = dtPrev;
  dtPrev  = dt;

  dt *= 0.5;

  if(dt < stack[0]*0.01)
  {
      throw runtime_error("Time step size is less than the minimum value");
  }

  stack.push_back(dt);

  return;
}



void myExpression::stck()
{ 
  cout << "time step stack" << endl;
  printVector(stack);

  int m = stack.size();

  if(m < 1)
    cerr << "myExpression::stck ... stack totally empty!" << endl;

  dtPrev2 = dtPrev;
  dtPrev  = dt;

  if(m == 1)
  {
    dt = stack[0];
    return;
  }

  if(m > 1)
  {
    dt = stack[m-1];
    stack.pop_back();

    return; 
  }

  return;
}



void myExpression::reset()
{
  dtPrev = -1.0;
  dt     = -1.0;
  cur    =  0.0;
  prev   = -1.0;
  prev2  = -2.0;
  write  =  1.0e9;
  dtOK   = false;
  stack.clear();

  return; 
}




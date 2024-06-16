
#include <iostream>

#include "TimeFunction.h"
#include "MyTime.h"
#include <assert.h>
#include <cmath>

extern MyTime myTime;


TimeFunction::TimeFunction()
{
  id = -1;
  nblocks = 0;
}



TimeFunction::~TimeFunction()
{
}


void  TimeFunction::addTimeBlock(vector<double> & ps_)
{
    assert(ps_.size() >= 10);

    t0.push_back(ps_[0]);
    tf.push_back(ps_[1]);

    int i, size = ps_.size()-2;
    std::vector<double>  doublelist(size);

    for(int i=0; i<size; ++i)
        doublelist[i] = ps_[i+2];

    ps.push_back(doublelist);

    nblocks = ps.size();

    return;
}



void TimeFunction::update()
{
  double t = myTime.cur, tmt0;
  std::vector<double>  p;

  tol = 1.e-12;

  factor = 0.;
  for(int i=0; i<nblocks; ++i)
  {
      if( (t < t0[i]-tol)  || (t > tf[i]-tol) )
        continue;

      p = ps[i];

      tmt0 = t-t0[i];

      factor += ( p[0] + p[1]*tmt0 + p[2]*sin(p[3]*tmt0+p[4]) + p[5]*cos(p[6]*tmt0+p[7]) );
  }

  return;
}



void TimeFunction::printSelf()
{
  for(int i=0; i<nblocks; i++)
  {
    std::cout << id << '\t' <<  t0[i] << '\t' << tf[i] << '\t'
    << ps[i][0] << '\t' << ps[i][1] << '\t' << ps[i][2] << '\t' << ps[i][3]
    << ps[i][4] << '\t' << ps[i][5] << '\t' << ps[i][6] << '\t' << ps[i][7] << std::endl;
  }
  std::cout << "\n\n\n" << std::endl;

  return;
}




double TimeFunction::evalFirstDerivative(double tcur)
{
    return 0.0;
}




double TimeFunction::evalSecondDerivative(double tcur)
{
    return 0.0;
}





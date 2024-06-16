
#ifndef incl_TimeFunction_h
#define incl_TimeFunction_h


#include  <vector>
using namespace std;

class TimeFunction
{
  private:
    int id, nblocks;

  public:

    double tol, factor;

    vector<double>  t0, tf;
    vector<vector<double> >  ps; // coefficients


    TimeFunction();

    ~TimeFunction();

    void setID(int id_)
    {
        id = id_; return;
    }

    void addTimeBlock(std::vector<double> & ps);

    void update();

    double evalFirstDerivative(double tcur=0.0);

    double evalSecondDerivative(double tcur=0.0);

    void printSelf();

    double  getValue()
    {
      return  factor;
    }

};



#endif


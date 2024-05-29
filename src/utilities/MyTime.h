
#ifndef incl_MyTime_h
#define incl_MyTime_h


#include <vector>


class MyTime
{ 
  public:
	
    MyTime() { reset(); return; }

    double cur, prev, prev2, dt, dtPrev, dtPrev2, dtMax, dtMin, dtDdtn, write;

    bool dtOK;

    std::vector<double> stack;

    void set(double dt_, double dtMin_, double dtMax_);

    void update();

    void cut();

    void stck();

    void reset();

};

#endif



#ifndef incl_myExpression_h
#define incl_myExpression_h


#include <vector>


class myExpression
{ 
  public:
	
    myExpression() { reset(); return; }

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



#ifndef incl_SolutionData_h
#define incl_SolutionData_h

#include "headersBasic.h"
#include "headersEigen.h"

using namespace std;

class  SolutionData
{
  private:
    int  size1, size2, size3, size4, size5;

  public:

    string  timeIntegrationScheme;
    bool firstIter, TRULY_INCOMPRESSIBLE, MassLumping;
    int  timeStepCount;
    double  spectralRadius;

    VectorXd  var1, var1Prev, var1Cur, var1Dot, var1DotPrev, var1DotCur, var1DotDot, var1DotDotPrev, var1DotDotCur, var1Incr, var1Extrap;
    VectorXd  var2, var2Prev, var2Cur, var2Incr;
    VectorXd  var3, var3Prev, var3Cur, var3Incr;
    VectorXd  var4, var4Prev, var4Cur, var4Incr, var4Dot, var4DotPrev, var4DotCur;
    VectorXd  var5, var5Prev, var5Cur, var5Incr, var5Dot, var5DotPrev, var5DotCur;
    VectorXd  var6, var6Prev, var6Cur, var6Incr, var6Dot, var6DotPrev, var6DotCur;
    VectorXd  dDot, dDotPrev;
    VectorXd  var1Init, var2Init, var3Init, var4Init, var5Init, var6Init;
    VectorXd  var1applied, var2applied, var3applied, var4applied, var5applied, var6applied;
    VectorXd  var1Prev2, var2Prev2, var3Prev2, var4Prev2, var5Prev2, var6Prev2;
    VectorXd  var1Prev3, var2Prev3, var3Prev3, var4Prev3, var5Prev3, var6Prev3;
    VectorXd  var1Prev4, var2Prev4, var3Prev4, var4Prev4, var5Prev4, var6Prev4;
    VectorXd  dispTarget;

    VectorXd  td, rhsVec, force, forcePrev, reac;
    VectorXd  soln, solnInit;

    //vector<PropertyItem*>  ElemProp, MatlProp;

    vector<int>  node_map_new_to_old;
    vector<int>  node_map_old_to_new;


    SolutionData();

    ~SolutionData(){}

    void setTimeIncrementType(string& ttt);

    void setSpectralRadius(double ttt);

    void  printSelf();

    void  initialise(int size1=0, int size2=0, int size3=0, int size4=0);

    void  setTimeParam();

    void  timeUpdate();

    void  updateIterStep();

    void  reset();

    void  saveSolution();
};





#endif



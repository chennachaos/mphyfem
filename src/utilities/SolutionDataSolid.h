
#ifndef incl_SolutionDataSolid_h
#define incl_SolutionDataSolid_h

#include "headersBasic.h"
#include "headersEigen.h"

using namespace std;

class  SolutionDataSolid
{
  private:
    int  size1, size2, size3, size4, size5;

  public:

    string  timeIntegrationScheme;
    bool firstIter, TRULY_INCOMPRESSIBLE, MassLumping;
    int  timeStepCount;
    double  spectralRadius;

    VectorXd  disp, dispPrev, dispCur, dispDot, dispDotPrev, dispDotCur, dispIncr, dispExtrap;
    VectorXd  velo, veloPrev, veloCur, veloIncr;
    VectorXd  acce, accePrev, acceCur, acceIncr;
    VectorXd  mpot, mpotPrev, mpotCur, mpotIncr;
    VectorXd  dispInit, veloInit, mpotInit;
    VectorXd  dispApplied, mpotApplied;
    VectorXd  dispPrev2, veloPrev2, mpotPrev2, dispPrev3;
    VectorXd  pres, presPrev, presPrev2, presCur, presIncr, presInit, presApplied;


    VectorXd  td, rhsVec, force, forcePrev, reac;

    vector<int>  node_map_new_to_old;
    vector<int>  node_map_old_to_new;


    SolutionDataSolid();

    ~SolutionDataSolid(){}

    void  setTimeIncrementScheme(string& ttt);

    string  getTimeIncrementScheme()
    {
      return  timeIntegrationScheme;
    }

    double  getSpectralRadius()
    {
      return  spectralRadius;
    }

    void  setSpectralRadius(double ttt);

    void  printSelf();

    void  initialise(int size1=0, int size2=0, int size3=0, int size4=0);

    void  setTimeParam();

    void  timeUpdate();

    void  updateIterStep();

    void  reset();

    void  saveSolution();
};





#endif



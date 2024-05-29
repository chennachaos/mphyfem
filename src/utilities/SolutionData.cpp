
#include "SolutionData.h"
#include "TimeFunction.h"
#include "MyTime.h"
#include "util.h"


extern MyTime myTime;
extern vector<unique_ptr<TimeFunction> > timeFunctions;



SolutionData::SolutionData()
{
    td.resize(100);

    size1 = size2 = size3 = size4 = size5 = 0;

    timeStepCount = 0;
    timeIntegrationScheme = "STATIC";
    spectralRadius = 0.0;

    MassLumping = false;

    return;
}



void SolutionData::setTimeIncrementType(string& ttt)
{
    timeIntegrationScheme = ttt;

    return;
}



void SolutionData::setSpectralRadius(double ttt)
{
    spectralRadius = ttt;

    if(spectralRadius < 0.0)
    {
        throw runtime_error("Spectral radius cannot be negative in SolutionData::setSpectralRadius");
    }

    return;
}



void SolutionData::initialise(int s1, int s2, int s3, int s4)
{
    size1 = s1;
    size2 = s2;
    size3 = s3;
    size4 = s4;

    var1.resize(size1);
    var1.setZero();

    var1Prev    = var1;
    var1Prev2   = var1;
    var1Prev3   = var1;
    var1Prev4   = var1;
    var1Cur     = var1;
    var1Incr    = var1;

    var1Dot     = var1;
    var1DotPrev = var1;
    var1DotCur  = var1;

    var1DotDot     = var1;
    var1DotDotPrev = var1;
    var1DotDotCur  = var1;

    var1Init    = var1;
    var1applied = var1;

    rhsVec   = var1;

    reac     = var1;

    dDot     = var1;
    dDotPrev = var1;


    var2.resize(size2);
    var2.setZero();

    var2Prev  = var2;
    var2Prev2 = var2;
    var2Prev3 = var2;
    var2Cur   = var2;
    var2Incr  = var2;
    var2Init  = var2;

    var3.resize(size3);
    var3.setZero();

    var3Prev  = var3;
    var3Cur   = var3;
    var3Incr  = var3;
    var3Init  = var3;

    var4.resize(size4);
    var4.setZero();

    var4Prev = var4;
    var4Cur  = var4;
    var4Incr = var4;

    var4Dot     = var4;
    var4DotPrev = var4;
    var4DotCur  = var4;

    var5.resize(size4);
    var5.setZero();

    var5Prev = var5;
    var5Cur  = var5;
    var5Incr = var5;

    var5Dot     = var5;
    var5DotPrev = var5;
    var5DotCur  = var5;

    return;
}



void  SolutionData::setTimeParam()
{
    // set the time integration parameters
    SetTimeParametersSolid(timeIntegrationScheme, spectralRadius, myTime.dt, td);

    return;
}



void  SolutionData::timeUpdate()
{
    // cout << "SolutionData::timeUpdate()" <<  endl;
    // store the variables

    timeStepCount++;

    // set the time integration parameters
    SetTimeParametersSolid(timeIntegrationScheme, spectralRadius, myTime.dt, td);


    //dtDdtn = 1.0; //myTime.dt/myTime.dtPrev;
    myTime.dtDdtn = myTime.dt/myTime.dtPrev;
    myTime.dtDdtn *= 0.5;

    //printf("time step = %14.10f \t %14.10f \t %14.10f \n", myTime.dt, myTime.dtPrev, myTime.dtDdtn);

    //myTime.dtDdtn = 0.0;

    //if(timeStepCount > 2)
    var1  = (1.0+myTime.dtDdtn)*var1Prev - myTime.dtDdtn*var1Prev2;
    //var2  = (1.0+myTime.dtDdtn)*var2Prev - myTime.dtDdtn*var2Prev2;
    //var3  = (1.0+myTime.dtDdtn)*var3Prev - myTime.dtDdtn*var3Prev2;
    //var4  = (1.0+myTime.dtDdtn)*var4Prev - myTime.dtDdtn*var4Prev2;
    //var5  = (1.0+myTime.dtDdtn)*var5Prev - myTime.dtDdtn*var5Prev2;

    //var1 = var1Prev;
    //var2 = var2Prev;
    //var3 = var3Prev;
    //var4 = var4Prev;
    //var5 = var5Prev;

    //else if(timeStepCount > 4)
      //var1 = 4.0*var1Prev - 6.0*var1Prev2 + 4.0*var1Prev3 - var1Prev4;

    return;
}



void SolutionData::updateIterStep()
{
  var1Dot       = td[10]*var1 + td[11]*var1Prev + td[12]*var1DotPrev + td[13]*var1DotDotPrev + td[14]*dDotPrev;
  var1DotDot    = td[15]*var1 + td[16]*var1Prev + td[17]*var1DotPrev + td[18]*var1DotDotPrev + td[19]*dDotPrev;
  dDot          = td[20]*var1 + td[21]*var1Prev + td[24]*dDotPrev;

  var1Cur       = td[2]*var1       + (1.0-td[2])*var1Prev;
  var1DotCur    = td[2]*var1Dot    + (1.0-td[2])*var1DotPrev;
  var1DotDotCur = td[1]*var1DotDot + (1.0-td[1])*var1DotDotPrev;

  var2Cur       = td[2]*var2    + (1.0-td[2])*var2Prev;     // pressure
  var3Cur       = td[2]*var3    + (1.0-td[2])*var3Prev;     // electric potential
  var4Cur       = td[2]*var4    + (1.0-td[2])*var4Prev;     // temperature
  var5Cur       = td[2]*var5    + (1.0-td[2])*var5Prev;     // chemical potential

  var4Dot       = td[20]*var4 + td[21]*var4Prev + td[24]*var4DotPrev;
  var5Dot       = td[20]*var5 + td[21]*var5Prev + td[24]*var5DotPrev;

  var4DotCur    = td[1]*var4Dot + (1.0-td[1])*var4DotPrev;  // temperature derivative

  var5DotCur    = td[1]*var5Dot + (1.0-td[1])*var5DotPrev;  // chemical potential derivative

  return;
}


void  SolutionData::reset()
{
  var1 = var1Prev;
  var2 = var2Prev;
  var3 = var3Prev;
  var4 = var4Prev;
  var5 = var5Prev;


  var1Dot = var1DotPrev ;
  var1DotDot = var1DotDotPrev;

  dDot = dDotPrev;

  var4Dot = var4DotPrev;
  var5Dot = var5DotPrev;

  force = forcePrev;

  return;
}


void SolutionData::saveSolution()
{
  // store the variables 

  var1Prev4      = var1Prev3;
  var1Prev3      = var1Prev2;
  var1Prev2      = var1Prev;
  var1Prev       = var1;

  var1DotPrev    = var1Dot;
  var1DotDotPrev = var1DotDot;
  dDotPrev       = dDot;
  forcePrev      = force;


  var2Prev4      = var2Prev3;
  var2Prev3      = var2Prev2;
  var2Prev2      = var2Prev;
  var2Prev       = var2;

  return;
}


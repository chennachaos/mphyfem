
#include "SolutionDataSolid.h"
#include "TimeFunction.h"
#include "MyTime.h"
#include "util.h"


extern MyTime myTime;
extern vector<unique_ptr<TimeFunction> > timeFunctions;



SolutionDataSolid::SolutionDataSolid()
{
    td.resize(100);

    size1 = size2 = size3 = size4 = size5 = 0;

    timeStepCount = 0;
    timeIntegrationScheme = "STEADY";
    spectralRadius = 0.0;

    MassLumping = false;

    return;
}



void SolutionDataSolid::setTimeIncrementScheme(string& ttt)
{
    timeIntegrationScheme = ttt;

    return;
}



void SolutionDataSolid::setSpectralRadius(double ttt)
{
    spectralRadius = ttt;

    if(spectralRadius < 0.0)
    {
        throw runtime_error("Spectral radius cannot be negative in SolutionDataSolid::setSpectralRadius");
    }

    return;
}



void SolutionDataSolid::initialise(int s1, int s2, int s3, int s4)
{
    size1 = s1;
    size2 = s2;
    size3 = s3;
    size4 = s4;

    disp.resize(size1);
    disp.setZero();

    dispPrev    = disp;
    dispPrev2   = disp;
    dispPrev3   = disp;
    dispCur     = disp;
    dispIncr    = disp;

    velo        = disp;
    veloPrev    = disp;
    veloCur     = disp;

    acce        = disp;
    accePrev    = disp;
    acceCur     = disp;

    dispInit    = disp;
    dispApplied = disp;

    dispDot     = disp;
    dispDotPrev = disp;

    rhsVec      = disp;
    reac        = disp;

    pres.resize(size2);
    pres.setZero();

    presPrev    = pres;
    presPrev2   = pres;
    presCur     = pres;
    presIncr    = pres;
    presInit    = pres;
    presApplied = pres;

    return;
}



void  SolutionDataSolid::setTimeParam()
{
    // set the time integration parameters
    SetTimeParametersSolid(timeIntegrationScheme, spectralRadius, myTime.dt, td);

    return;
}



void  SolutionDataSolid::timeUpdate()
{
    // cout << "SolutionDataSolid::timeUpdate()" <<  endl;
    // store the variables

    timeStepCount++;

    // set the time integration parameters
    SetTimeParametersSolid(timeIntegrationScheme, spectralRadius, myTime.dt, td);


    //dtDdtn = 1.0; //myTime.dt/myTime.dtPrev;
    myTime.dtDdtn = myTime.dt/myTime.dtPrev;

    //printf("time step = %14.10f \t %14.10f \t %14.10f \n", myTime.dt, myTime.dtPrev, myTime.dtDdtn);

    //myTime.dtDdtn = 0.0;

    //if(timeStepCount > 2)
    //disp  = (1.0+myTime.dtDdtn)*dispPrev - myTime.dtDdtn*dispPrev2;
    //pres  = (1.0+myTime.dtDdtn)*presPrev - myTime.dtDdtn*presPrev2;
    //var3  = (1.0+myTime.dtDdtn)*var3Prev - myTime.dtDdtn*var3Prev2;
    //var4  = (1.0+myTime.dtDdtn)*var4Prev - myTime.dtDdtn*var4Prev2;
    //var5  = (1.0+myTime.dtDdtn)*var5Prev - myTime.dtDdtn*var5Prev2;

    return;
}



void SolutionDataSolid::updateIterStep()
{
  velo      = td[10]*disp + td[11]*dispPrev + td[12]*veloPrev + td[13]*accePrev + td[14]*dispDotPrev;
  acce      = td[15]*disp + td[16]*dispPrev + td[17]*veloPrev + td[18]*accePrev + td[19]*dispDotPrev;
  dispDot   = td[20]*disp + td[21]*dispPrev + td[22]*veloPrev + td[23]*accePrev + td[24]*dispDotPrev;

  dispCur   = td[2]*disp  + (1.0-td[2])*dispPrev;
  veloCur   = td[2]*velo  + (1.0-td[2])*veloPrev;
  acceCur   = td[1]*acce  + (1.0-td[1])*accePrev;

  presCur   = td[2]*pres  + (1.0-td[2])*presPrev;

  return;
}


void  SolutionDataSolid::reset()
{
  disp    = dispPrev;
  velo    = veloPrev;
  acce    = accePrev;
  pres    = presPrev;
  dispDot = dispDotPrev;
  force   = forcePrev;

  return;
}


void SolutionDataSolid::saveSolution()
{
  // store the variables

  dispPrev2      = dispPrev;
  dispPrev       = disp;
  veloPrev2      = veloPrev;
  veloPrev       = velo;
  accePrev       = acce;
  dispDotPrev    = dispDot;
  forcePrev      = force;
  presPrev       = pres;

  return;
}


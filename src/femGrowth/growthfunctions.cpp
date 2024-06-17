
#include  "growthfunctions.h"
#include  "myConstants.h"

using namespace std;


int  growthfunction_isotropic(int  ftype, int  ndim,  double  growthFactor,  MatrixXd&  Fg)
{
    double  grp1 = 1.0+growthFactor;

    switch(ftype)
    {
        case 1:

            Fg(0,0) = grp1;    Fg(0,1) = 0.0;     Fg(0,2) = 0.0;
            Fg(1,0) = 0.0;     Fg(1,1) = 1.0;    Fg(1,2) = 0.0;
            Fg(2,0) = 0.0;     Fg(2,1) = 0.0;     Fg(2,2) = grp1;

            if(ndim == 2)  Fg(2,2) = 1.0;

        break;

        case 2:

            Fg(0,0) = grp1;    Fg(0,1) = 0.0;     Fg(0,2) = 0.0;
            Fg(1,0) = 0.0;     Fg(1,1) = grp1;    Fg(1,2) = 0.0;
            Fg(2,0) = 0.0;     Fg(2,1) = 0.0;     Fg(2,2) = grp1;

            if(ndim == 2)  Fg(2,2) = 1.0;

        break;

        case 3:

            Fg(0,0) = grp1;    Fg(0,1) = 0.0;     Fg(0,2) = 0.0;
            Fg(1,0) = 0.0;     Fg(1,1) = grp1;    Fg(1,2) = 0.0;
            Fg(2,0) = 0.0;     Fg(2,1) = 0.0;     Fg(2,2) = grp1;

            if(ndim == 2)  Fg(2,2) = 1.0;

        break;

        default:
            cout << "growthfunction_isotropic not defined for ftype = " << ftype << endl;
        break;
    }

    return 0;
}



int  growthfunction_anisotropic(int  ftype, int  ndim,  MatrixXd&  F)
{
    return 0;
}


int  growthfunction_beam_singlelayer(double  timeFact,  double* geom,  MatrixXd&  Fg)
{
    double  lambda10=1.0, lambda30=1.0, beamlength=1.0, lambda1, lambda3, lambda11;

    /*
    if(timeFact <= 3.00001)
    {
      lambda10 = 1.0;
      lambda11 = timeFact*PI*0.5;
      lambda30 = 1.0;
    }
    else
    {
      lambda10 = 1.0;
      lambda11 = PI*0.5;
      lambda30 = 1.0 + sin(0.5*PI*(timeFact-1.0));
    }

    Fg(0,0) = lambda10 + lambda11*geom[1]/beamlength;
    Fg(1,1) = lambda30;
    Fg(2,2) = 1.0;
    */
    int type = 1;
    double  val1, val2, X = geom[0], Y = geom[1], h = 0.02;

    Fg.setZero();

    switch(type)
    {
        case 1: // circular curve

            //Fg(0,0) = 1.0 + (0.5*PI + Y*PI - 1.0) *timeFact;
            Fg(0,0) = 1.0 + (PI*Y) *timeFact;
            Fg(1,1) = 1.0;
            Fg(2,2) = 1.0;

        break;

        case 2: // butterfly curve

            val1 = 0.5*PI*sqrt( cos(PI*geom[0])*cos(PI*geom[0]) + 4.0*cos(2.0*PI*geom[0])*cos(2.0*PI*geom[0]) );
            val2 = 2.0*PI*(3.0*sin(PI*geom[0])+sin(3.0*PI*geom[0]))/(5.0+cos(2.0*PI*geom[0])+4.0*cos(4.0*PI*geom[0]));

            Fg(0,0) = 1.0 + (val1 + Y*val2 - 1.0)*timeFact;
            Fg(1,1) = 1.0;
            Fg(2,2) = 1.0;

        break;

        case 3:  // elliptic curve

            val1 = 5.0+3.0*cos(2.0*PI*X);

            lambda10 = PI*sqrt(val1)/3.0/sqrt(2.0);
            lambda11 = 2.0*sqrt(2.0)*PI*PI/3.0/sqrt(val1);
            lambda11 += 32*sqrt(2)*h*h*PI*PI*PI*PI*(17+15*cos(2*PI*X))/27/pow(val1, 2.5);
            lambda11 += sqrt(2)*h*PI*PI*PI*(215+252*cos(2*PI*X)+45*cos(4*PI*X))/9/pow(val1, 2.5);

            Fg(0,0) = 1.0 + (lambda10 + Y*lambda11 - 1.0)*timeFact;
            Fg(1,1) = Fg(0,0);
            Fg(2,2) = 1.0;

        break;

        case 4:  // helical curve

            val1 = 26.0+25.0*X*(2.0+X);

            lambda10 = 0.5*sqrt(val1);
            lambda11 = 2.5*(27.0+25.0*X*(2.0+X))/sqrt(val1);
            //lambda11 += 25.0*h*( 753.0 + 125.0*X*(2.0+X)*(11.0+5.0*X*(2.0+X)) )/3.0/pow(val1, 1.5);
            //lambda11 += 500.0*h*h*(27.0+25.0*X*(2.0+X))*( 753.0 + 125.0*X*(2.0+X)*(11.0+5.0*X*(2.0+X)) )/9.0/pow(val1, 2.5);

            Fg(0,0) = 1.0 + (lambda10 + Y*lambda11 - 1.0)*timeFact;
            Fg(1,1) = Fg(0,0);
            Fg(2,2) = 1.0;

        break;

        default:

            cout << "This growth function not defined yet " << endl;

        break;
    }

    return 0;
}


int  growthfunction_beam_singlelayer3D(double  timeFact,  double* geom,  MatrixXd&  Fg)
{
    double  lambda10=1.0, lambda30=1.0, beamlength=1.0, lambda1, lambda3, lambda11;
    // growth tensors for thin circular rod example


    double  val1, val2, X = geom[0], Y = geom[1], Z = geom[2];

    Fg.setZero();

    /*
    helical shape 1
    Fg(0,0) = 1.0;
    Fg(1,1) = 1.0;
    Fg(2,2) = 1.0 - (0.5*PI*Y)*timeFact ;
    Fg(0,2) = timeFact*0.1;
    */

    //
    // circular shape
    Fg(0,0) = 1.0;
    Fg(1,1) = 1.0 ;
    Fg(2,2) = 1.0 - 0.5*PI*X*timeFact ;
    //

    /*
    //spiral shape
    Fg(0,0) = 1.0 ;
    Fg(1,1) = 1.0 ;
    Fg(2,2) = 1.0 - Z*(0.5*PI*X)*timeFact;
    */

    /*
    Fg(0,0) = 1.0;
    Fg(0,0) += (Y*PI*0.0+(Z-0.02)*2) *timeFact;
    Fg(1,1) = 1.0;
    Fg(2,2) = 1.0;
    */

    /*
    // helical curling
    //Fg(0,0) = 1.0 + (Y*PI+(Z-0.02)*2) *timeFact;
    Fg(0,0) = 1.0;
    Fg(1,1) = 1.0;
    Fg(2,2) = 1.0;
    Fg(2,2) += timeFact*Z*(X*cos(PI*Z)+Y*sin(PI*Z));
    //Fg(2,2) += timeFact*Z*(X*cos(PI*Z)+Y*sin(PI*Z));
    */

    /*
    // helical curling
    //Fg(0,0) = 1.0 + (Y*PI+(Z-0.02)*2) *timeFact;
    Fg(0,0) = 1.0;
    Fg(1,1) = 1.0;
    Fg(2,2) = 1.0;
    Fg(2,2) += timeFact*Z*(X*cos(PI*Z) );
    //Fg(2,2) += timeFact*Z*(X*cos(PI*Z)+X*sin(PI*Z)+Y*sin(PI*Z)+Y*cos(PI*Z));
    */

    return 0;
}



int  growthfunction_circularplate_singlelayer(double  timeFact, double* geom,  MatrixXd&  Fg)
{
    double  lambdaR=1.0, lambdaT=1.0, beamlength=1.0;
    // Axis is along the Z-direction
    double  radSq = geom[0]*geom[0] + geom[1]*geom[1], rad = sqrt(radSq);

/*
    if(timeFact <= 1.00001)
    {
      lambdaR = 1.0;
      lambdaT = 1.0+timeFact*0.3;
    }
    else
    {
      lambdaR = 1.0+(timeFact-1.0)*0.3;
      lambdaT = 1.3;
    }

    //lambdaR = 1.0+timeFact*0.1*exp(radSq-4.0/9.0);
    //lambdaT = 1.0+timeFact*0.1;
*/

    lambdaR = 1.0 + timeFact*(PI-1.0);
    lambdaT = 1.0 + timeFact*(6.0*sin(5.0*PI*rad/6.0)/5.0/rad - 1.0);

    //lambdaR = 1.0;
    //lambdaT = 1.0;
    //lambdaR = 1.0+timeFact*sin(PI*rad)*(0.04-geom[2]);
    lambdaR = 1.0+timeFact*(0.04-geom[2]);
    lambdaT = 1.0;

    Fg.setZero();
    Fg(0,0) = (geom[0]*geom[0]*lambdaR + geom[1]*geom[1]*lambdaT)/radSq;
    Fg(0,1) =  geom[0]*geom[1]*(lambdaR - lambdaT)/radSq;

    Fg(1,0) = Fg(0,1);
    Fg(1,1) = (geom[1]*geom[1]*lambdaR + geom[0]*geom[0]*lambdaT)/radSq;

    //Fg(2,2) = 1.0+timeFact*sqrt(radSq);
    Fg(2,2) = 1.0+timeFact*geom[2]*sqrt(radSq);

    return 0;
}


int  growthfunction_beam_multilayer(double  timeFact,  double* geom,  MatrixXd&  Fg)
{
    int  nlayers = 4;
    double  dlambda1, dlambda2, dlambda3, lambda1=1.0;
    double   h = 0.1; // total thickness
    double  dh = 0.1/nlayers; // thickness of each layer

    if(nlayers == 2)
    {
      dlambda1 = 0.0;
      dlambda2 = 2.0*PI/3.0;

      if(geom[1] < -0.025)
      {
        //lambda1 = 1.0;
        lambda1 = 1.0 + (53.0*PI/240.0 - 1.0) *timeFact;
      }
      else
      {
        //lambda1 = 1.0 + timeFact*dlambda2*dh;
        lambda1 = 1.0 + (61.0*PI/240.0 - 1.0) *timeFact;
      }
    }
    else if(nlayers == 3)
    {
      dlambda1 = 0.0;
      dlambda2 = 9.0*PI/16.0;
      dlambda3 = 9.0*PI/8.0;

      if(geom[1] < -0.05)
      {
        lambda1 = 1.0;
      }
      else if( (geom[1] > -0.05) && (geom[1] < -0.01666666666666667) )
      {
        lambda1 = 1.0 + timeFact*dlambda2*dh;
      }
      else
      {
        lambda1 = 1.0 + timeFact*dlambda3*dh;
      }
    }
    else if(nlayers == 4)
    {
      if(geom[1] < -0.0625)
      {
        lambda1 = 1.0 + (1.1+PI/80.0 - 1.0) *timeFact;
      }
      else if( (geom[1] > -0.0625) && (geom[1] < -0.0375) )
      {
        lambda1 = 1.0 + (1.1+3.0*PI/80.0 - 1.0) *timeFact;
      }
      else if( (geom[1] > -0.0375) && (geom[1] < -0.0125) )
      {
        lambda1 = 1.0 + (1.1+5.0*PI/80.0 - 1.0) *timeFact;
      }
      else
      {
        lambda1 = 1.0 + (1.1+7.0*PI/80.0 - 1.0) *timeFact;
      }
    }

    Fg.setZero();
    Fg(0,0) = lambda1;
    Fg(1,1) = 1.0;
    Fg(2,2) = 1.0;

    return 0;
}


int  growthfunction_rod(int ftype, int ndim, double timeFact, double* geom,  Eigen::MatrixXd& Fg)
{
    Fg(0,0) = 1.0;     Fg(0,1) = 0.0;     Fg(0,2) = 0.0;
    Fg(1,0) = 0.0;     Fg(1,1) = 1.0;     Fg(1,2) = 0.0;
    Fg(2,0) = 0.0;     Fg(2,1) = 0.0;     Fg(2,2) = 1.0;

    // curvature only
    //Fg(2,2) += timeFact*(geom[0]+geom[1]);
    // curvature and twist
    if(ftype == 1) // circle
    {
      Fg(2,2) -= timeFact*0.5*PI*geom[0];
    }
    if(ftype == 2) // spiral
    {
      Fg(2,2) -= timeFact*0.5*PI*geom[0]*geom[2];
    }
    if(ftype == 3) // helical
    {
      Fg(0,2) += timeFact*0.4;
      Fg(2,2) -= timeFact*2.0*PI*geom[1];
    }
    else if(ftype == 4)
    {
      double  phi = PI/3.0;
      Fg(2,2) += timeFact*geom[2]*(geom[0]*cos(2.0*PI*geom[2]/40.0)+geom[1]*sin(2.0*PI*geom[2]/40.0));
    }

    return 0;
}





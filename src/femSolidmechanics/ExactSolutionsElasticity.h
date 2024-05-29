
#ifndef  EXACT_SOLN_ELASTICITY_H
#define  EXACT_SOLN_ELASTICITY_H



class PlateWithHole
{
  private:
    double  T, R, nu, G, E, K;

  public:

    PlateWithHole(int sss, double G1, double nu1)
    {
      G  = G1;
      nu = nu1;
      
      E = 2.0*G*(1.0+nu);

      T = 10.0;
      R = 1.0;
      
      if(sss==1)// plane-stress
      {
        K = (3.0-nu)/(1.0+nu);
      }
      else //if(sss==2)// plane-strain 
      {
        K = 3.0-4.0*nu;
      }
    }
    
    ~PlateWithHole() {}
    
    double  dispX(double r, double t)
    {
      return   (T*R/8.0/G)*(r*(K+1.0)*cos(t)/R + 2.0*R*((1.0+K)*cos(t)+cos(3.0*t))/r - 2.0*R*R*R*cos(3.0*t)/r/r/r);
    }

    double  dispY(double r, double t)
    {
      return  (T*R/8.0/G)*(r*(K-3.0)*sin(t)/R + 2.0*R*((1.0-K)*sin(t)+sin(3.0*t))/r - 2.0*R*R*R*sin(3.0*t)/r/r/r);
    }
    
    double  forceX(double r, double t)
    {
       return 0.0;
    }
    double  forceY(double r, double t)
    {
       return 0.0;
    }

    double  stressXX(double r, double t)
    {
      return  (T - T*R*R*(1.5*cos(2.0*t)+cos(4.0*t))/r/r + T*1.5*R*R*R*R*cos(4.0*t)/r/r/r/r);
    }

    double  stressYY(double r, double t)
    {
      return  ( - T*R*R*(0.5*cos(2.0*t)-cos(4.0*t))/r/r - T*1.5*R*R*R*R*cos(4.0*t)/r/r/r/r);
    }

    double  stressXY(double r, double t)
    {
      return  (- T*R*R*(0.5*sin(2.0*t)+sin(4.0*t))/r/r + T*1.5*R*R*R*R*sin(4.0*t)/r/r/r/r);
    }

    void  stresses(double r, double t, double* val)
    {
      val[0] = T - T*R*R*(1.5*cos(2.0*t)+cos(4.0*t))/r/r + T*1.5*R*R*R*R*cos(4.0*t)/r/r/r/r ;

      val[1] = - T*R*R*(0.5*cos(2.0*t)-cos(4.0*t))/r/r - T*1.5*R*R*R*R*cos(4.0*t)/r/r/r/r ;

      val[2] = - T*R*R*(0.5*sin(2.0*t)+sin(4.0*t))/r/r + T*1.5*R*R*R*R*sin(4.0*t)/r/r/r/r ;

      return;
    }

    double  strainXX(double r, double t)
    {
      return  ((T*R/8.0/G)*((1.0+K)/R + 6.0*R*R*R*cos(4.0*t)/r/r/r/r - 2.0*R*(K*cos(2.0*t)+2.0*cos(4.0*t))/r/r ));
    }

    double  strainYY(double r, double t)
    {
      return  ((T*R/8.0/G)*((K-3.0)/R - 6.0*R*R*R*cos(4.0*t)/r/r/r/r - 2.0*R*(K*cos(2.0*t)-2.0*cos(4.0*t)-2.0*cos(2.0*t))/r/r ));
    }

    double  strainXY(double r, double t)
    {
      return  ((T*R/8.0/G)*( 6.0*R*R*R*sin(4.0*t)/r/r/r/r - 2.0*R*(sin(2.0*t)+2.0*sin(4.0*t))/r/r ));
    }

    void  strains(double r, double t, double* val)
    {
      val[0] =  ((T*R/8.0/G)*((1.0+K)/R + 6.0*R*R*R*cos(4.0*t)/r/r/r/r - 2.0*R*(K*cos(2.0*t)+2.0*cos(4.0*t))/r/r ));
      val[1] =  ((T*R/8.0/G)*((K-3.0)/R - 6.0*R*R*R*cos(4.0*t)/r/r/r/r - 2.0*R*(K*cos(2.0*t)-2.0*cos(4.0*t)-2.0*cos(2.0*t))/r/r ));
      val[2] =  ((T*R/8.0/G)*( 6.0*R*R*R*sin(4.0*t)/r/r/r/r - 2.0*R*(sin(2.0*t)+2.0*sin(4.0*t))/r/r ));

      return;
    }

};




class LaplaceNew
{
  public:

  LaplaceNew(){}
  
  ~LaplaceNew() {}
  
  double  value(double  x, double y)
  {
    return (x*x*x*(1.0-x*x*x)*y*y*y*(1.0-y*y*y));
  }

  double  force(double  x, double y)
  {
    return -6.0*( x*(1.0-5.0*x*x*x)*y*y*y*(1.0-y*y*y) + x*x*x*(1.0-x*x*x)*y*(1.0-5.0*y*y*y));
  }

  void  gradient(double  x, double y, double* grad)
  {
    grad[0] =  3.0*x*x*(1.0-2.0*x*x*x)*y*y*y*(1.0-y*y*y);
    grad[1] =  3.0*x*x*x*(1.0-x*x*x)*y*y*(1.0-2.0*y*y*y) ;
    return ;
  }

};



class  ElasticityEx1
{
private:
  double  mu;
public:
  
  ElasticityEx1(double mu1)
  {
    mu = mu1;
  }
  
  ~ElasticityEx1() {}
  
  double  dispX(double x, double y)
  {
    double fact = x*x*y*y*y*y * (x*x+y*y-16.0) * (x*x+y*y-1.0) ;

    fact *= (5*x*x*x*x + 18.0*x*x*y*y - 85.0*x*x + 13.0*y*y*y*y + 80.0 - 153.0*y*y);
    
    fact /= 1.0e6;

    return fact;
  }
  
  double  dispY(double x, double y)
  {
    double fact = -2.0*( x*y*y*y*y*y * (x*x+y*y-16.0) * (x*x+y*y-1.0) );

    fact *= (5*x*x*x*x - 51.0*x*x + 6.0*x*x*y*y - 17.0*y*y + 16.0 + y*y*y*y);
    
    fact /= 1.0e6;

    return fact;
  }

  double  forceX(double x, double y)
  {
    double fact = 91485.0*pow(x,4.0)*pow(y,2.0) + 2296.0*pow(x,6.0)*pow(y,4.0) + 2790.0*pow(x,4.0)*pow(y,6.0) ;
    
    fact += 7680.0*pow(x,2.0) + 645.0*pow(x,8.0)*pow(y,2.0) - 15470.0*pow(x,6.0)*pow(y,2.0) - 3808.0*pow(y,4.0) + 2889.0*pow(y,6.0);
    
    fact += 1280.0*pow(y,2.0) - 36414.0*pow(x,4.0)*pow(y,4.0) + 107856.0*pow(x,2.0)*pow(y,4.0) + 13.0*pow(y,10.0) - 374.0*pow(y,8.0);

    fact += 1122.0*pow(y,8.0)*pow(x,2.0) - 22338.0*pow(y,6.0)*pow(x,2.0) - 16320.0*pow(x,4.0) - 73440.0*pow(x,2.0)*pow(y,2.0);
    
    fact += 9630.0*pow(x,6.0) - 1020.0*pow(x,8.0) + 30.0*pow(x,10.0);
    
    fact *= -2.0*y*y*mu/1.0e6;

    return fact;
  }
  
  double  forceY(double x, double y)
  {
    double fact = 770.0*pow(y,4.0) + 1280.0 + 258.0*pow(x,6.0)*pow(y,2.0) + 492.0*pow(x,4.0)*pow(y,4.0) ;
    
    fact += -4641.0*pow(x,4.0)*pow(y,2.0) + 310.0*pow(y,6.0)*pow(x,2.0) - 5202.0*pow(x,2.0)*pow(y,4.0) + 25.0*pow(x,8.0) - 680.0*pow(x,6.0);
    
    fact += 4815.0*pow(x,4.0) + 51.0*pow(y,8.0) - 1241.0*pow(y,6.0) + 18297.0*pow(x,2.0)*pow(y,2.0) - 5440.0*pow(x,2.0) - 7334.0*pow(y,2.0);

    fact *= 8.0*y*y*y*mu/1.0e6;

    return fact;
  }

};





class  ElasticityEx2
{
private:
  double  G, K, L;
public:
  
  ElasticityEx2(double K1, double G1)
  {
    G = G1;
    K = K1;
    L = K - 2.0*G/3.0;
  }
  
  ~ElasticityEx2() {}
  
  double  dispX(double x, double y)
  {
    double fact = x*x*y*y*y*y * (x*x+y*y-16.0) * (x*x+y*y-1.0)/1.0e6 ;

    return fact;
  }
  
  double  dispY(double x, double y)
  {
    double fact = -2.0*( x*y*y*y*y*y * (x*x+y*y-16.0) * (x*x+y*y-1.0) )/1.0e6;

    return fact;
  }

  double  forceX(double x, double y)
  {
    double fact = 102.0*G*pow(x,4.0) - 96.0*G*pow(x,2.0) - 6.0*G*pow(x,6.0) + 48.0*G*pow(y,2.0) ;
    fact += - 85.0*G*pow(y,4.0) + 7.0*G*pow(y,6.0) + 64.0*L*pow(y,2.0) - 102.0*L*pow(y,4.0) ;
    fact += 8.0*L*pow(y,6.0) + 204.0*G*pow(x,2.0)*pow(y,2.0) - 10.0*G*pow(x,2.0)*pow(y,4.0);
    fact += - 35.0*G*pow(x,4.0)*pow(y,2.0) - 153.0*L*pow(x,2.0)*pow(y,2.0) ;
    fact += 30.0*L*pow(x,2.0)*pow(y,4.0) + 10.0*L*pow(x,4.0)*pow(y,2.0);

    fact *= y*y/500000.0;

    return fact;
  }
  
  double  forceY(double x, double y)
  {
    double  fact = 144.0*G + 64.0*L - 136.0*G*pow(x,2.0) + 7.0*G*pow(x,4.0) - 357.0*G*pow(y,2.0) + 37.0*G*pow(y,4.0);
    fact += -51.0*L*pow(x,2.0) + 2.0*L*pow(x,4.0) - 153.0*L*pow(y,2.0) + 16.0*L*pow(y,4.0);
    fact +=  41.0*G*pow(x,2.0)*pow(y,2.0) + 15.0*L*pow(x,2.0)*pow(y,2.0);

    fact *= x*pow(y,3.0)/125000.0;

    return fact;
  }
  
    double  stressXX(double x, double y)
    {
      return  (x*pow(y,4.0)*(32.0*G - 64.0*L - 68.0*G*x*x + 6.0*G*pow(x,4.0) - 34.0*G*y*y + 2*G*pow(y,4.0) + 51.0*L*x*x - 2.0*L*pow(x,4.0) + 102.0*L*y*y - 8.0*L*pow(y,4.0) + 8.0*G*x*x*y*y - 10*L*x*x*y*y))/500000.0;
    }

    double  stressYY(double x, double y)
    {
      return  -(x*pow(y,4.0)*(160.0*G + 64.0*L - 170.0*G*pow(x,2.0) + 10.0*G*pow(x,4.0) - 238.0*G*y*y + 18.0*G*pow(y,4.0) - 51.0*L*x*x + 2.0*L*pow(x,4.0) - 102.0*L*y*y + 8.0*L*pow(y,4.0) + 28.0*G*x*x*y*y + 10.0*L*x*x*y*y))/500000.0;
    }

    double  stressXY(double x, double y)
    {
      return  -(G*pow(y,3.0)*(- 2*pow(x,6.0) - pow(x,4.0)*pow(y,2.0) + 34.0*pow(x,4.0) + 2*pow(x,2.0)*pow(y,4.0) - 32*x*x + pow(y,6.0) - 17*pow(y,4.0) + 16*y*y))/500000.0;
    }

};




class  ElasticityEx3
{
private:
  double  G, K, L;
public:
  
  ElasticityEx3(double K1, double G1)
  {
    G = G1;
    K = K1;
    L = K - 2.0*G/3.0;
  }
  
  ~ElasticityEx3() {}
  
  double  dispX(double x, double y)
  {
    return  (x*x* y*y * (x-2.0)*(x-2.0) * (y-1.0) );
  }
  
  double  dispY(double x, double y)
  {
    return  ( x*x* y*y * (x-2.0)   * (y-1.0)*(y-1.0) );
  }

  double  forceX(double x, double y)
  {
    double  fact = 8.0*G*pow(x,2.0) - 8.0*G*pow(x,3.0) + 2.0*G*pow(x,4.0) + 16.0*G*pow(y,2.0) - 16.0*G*pow(y,3.0);
    fact  +=  8.0*L*pow(y,2.0) - 8.0*L*pow(y, 3.0) + 42.0*G*pow(x,2.0)*pow(y,2.0) - 36.0*G*pow(x,2.0)*pow(y,3.0);
    fact  +=  30.0*L*pow(x,2.0)*pow(y,2.0) - 24.0*L*pow(x,2.0)*pow(y,3.0) + 8.0*G*x*y + 8.0*L*x*y ;
    fact  +=  -72.0*G*x*pow(y,2.0) - 30.0*G*pow(x,2.0)*y + 64.0*G*x*pow(y,3.0) + 24.0*G*pow(x,3.0)*y ;
    fact  +=  -6.0*G*pow(x,4.0)*y - 48.0*L*x*pow(y,2.0) - 6.0*L*pow(x,2.0)*y + 40.0*L*x*pow(y,3.0);

    return fact;
  }
  
  double  forceY(double x, double y)
  {
    double fact = 8.0*G*pow(x,2.0) - 4.0*G*pow(x,3.0) + 4.0*G*pow(y,2.0) - 8.0*G*pow(y,3.0) + 4.0*G*pow(y,4.0);
    fact += 4.0*L*pow(x,2.0) - 2.0*L*pow(x,3.0) + 84.0*G*pow(x,2.0)*pow(y,2.0);
    fact += - 36.0*G*pow(x,3.0)*pow(y,2.0) + 60.0*L*pow(x,2.0)*pow(y,2.0) - 24.0*L*pow(x,3.0)*pow(y,2.0);
    fact += 16.0*G*x*y + 16.0*L*x*y - 30.0*G*x*pow(y,2.0) - 72.0*G*pow(x,2.0)*y + 12.0*G*x*pow(y,3.0);
    fact += 32.0*G*pow(x,3.0)*y - 6.0*G*x*pow(y,4.0) - 24.0*L*x*pow(y,2.0) - 48.0*L*pow(x,2.0)*y + 20.0*L*pow(x,3.0)*y;

    return fact;
  }

  double  strainXX(double x, double y)
  {
    return  ( 4.0*x*y*y*(x-1.0)*(x-2.0)*(y-1.0) );
  }
  
  double  strainYY(double x, double y)
  {
    return  ( 2.0*x*x*y*(2.0*y-1.0)*(x-2.0)*(y-1.0) );
  }

  double  strainXY(double x, double y)
  {
    return  ( 0.5*( x*x*(x-2.0)*(x-2.0)*(3.0*y*y-2.0*y) + (3.0*x*x-4.0*x)*y*y*(y-1.0)*(y-1.0) ) );
  }

  double  stressXX(double x, double y)
  {
    return  ( (2.0*G+L)*strainXX(x,y) + L*strainYY(x,y) );
  }
  
  double  stressYY(double x, double y)
  {
    return  ( (2.0*G+L)*strainYY(x,y) + L*strainXX(x,y) );
  }

  double  stressXY(double x, double y)
  {
    return  ( 2.0*G*strainXY(x, y) );
  }
  
};





class  ThickCylinder
{
  private:
    double  pin, ri, ro, G, E, nu, fact, a, b1, b2, A;

  public:

    ThickCylinder(int sss, double E1, double nu1)
    {
      //E = 200.0;
      //nu = 0.3;
      E = E1;
      nu = nu1;

      ri = 100.0;
      ro = 200.0;

      pin = 1.0;

      fact = pin*ri*ri/(ro*ro-ri*ri);
      A = fact;

      a = fact/E;

      b2 = (1.0+nu)*ro*ro;

      if(sss==1)// plane-stress
      {
        b1 = 1.0 - nu ;
      }
      else //if(sss==2)// plane-strain
      {
        b1 = 1.0 - nu - 2.0*nu*nu;
      }
    }

    ~ThickCylinder() {}

    double  forceX(double x, double y)
    { return 0.0; }

    double  forceY(double x, double y)
    { return 0.0; }

    double  forceR(double r, double t)
    { return 0.0; }

    double  dispR(double r, double t)
    {
      return   ( a*( b1*r + b2/r) );
    }

    double  dispX(double r, double t)
    {
      return   ( a*( b1*r + b2/r)*cos(t) );
    }

    double  dispY(double r, double t)
    {
      return   ( a*( b1*r + b2/r)*sin(t) );
    }

    //double  dispX(double  x, double y)
    //{
      //return   ( a*( b1 + b2/(x*x+y*y))*x );
    //}

    //double  dispY(double  x, double y)
    //{
      //return   ( a*( b1 + b2/(x*x+y*y))*y );
    //}

    double  stressRR(double r, double t)
    {
      return  (fact*(1.0-ro*ro/r/r));
    }

    double  stressTT(double r, double t)
    {
      return  (fact*(1.0+ro*ro/r/r));
    }

    double  stressRT(double r, double t)
    {
      return  0.0;
    }

    double  stressZZ(double r, double t)
    {
      return  nu*(stressRR(r, t)+stressTT(r, t));
    }

    double  pressure(double r, double t)
    {
      return  2*(1.0+nu)*A/3.0;
    }
    //double  stressXX(double x, double y)
    //{
      //return  fact*(1.0 - ro*ro*(x*x-y*y)/(x*x+y*y)/(x*x+y*y));
    //}

    //double  stressYY(double x, double y)
    //{
      //return  fact*(1.0 + ro*ro*(x*x-y*y)/(x*x+y*y)/(x*x+y*y));
    //}

    //double  stressXY(double x, double y)
    //{
      //return  fact*( - ro*ro*(2.0*x*y)/(x*x+y*y)/(x*x+y*y));
    //}

    double  stressXX(double r, double t)
    {
      return  fact*(1.0 - ro*ro*cos(2.0*t)/r/r);
    }

    double  stressYY(double r, double t)
    {
       return  fact*(1.0 + ro*ro*cos(2.0*t)/r/r);
    }

    double  stressXY(double r, double t)
    {
      return  fact*(- ro*ro*sin(2.0*t)/r/r);
    }

    //void  stresses(double r, double t, double* val)
    //{
      //val[0] = fact*(1.0 - ro*ro*cos(2.0*t)/r/r);
      //val[1] = fact*(1.0 + ro*ro*cos(2.0*t)/r/r);
      //val[2] = fact*(- ro*ro*sin(2.0*t)/r/r);
    //}

    void  stresses1(double r, double t, double* val)
    {
      double srr = stressRR(r, t), stt = stressTT(r, t), srt = 0.0;
      double ct = cos(t), st = sin(t);

      val[0] = srr*ct*ct + stt*st*st - 2*srt*st*ct; // xx
      val[1] = srr*st*st + stt*ct*ct + 2*srt*st*ct; // yy
      val[2] = (srr-stt)*st*ct + srt*(ct*ct-st*st); // xy
      val[3] = nu*(srr+stt); //zz
    }

    void  stresses2(double x, double y, double* val)
    {
      /*
      // ur = a*( b1*r + b2/r);

      double  du[3][3];
      double  a1 = a*b1;
      double  a2 = a*b2;
      double  r = sqrt(x*x+y*y);
      double  rP3 = r*r*r;
      double  rP5 = r*r*r*r*r;

      du[0][0] = -3.0*a2*x*x/rP5 + a1 + a2/rP3;
      du[0][1] = -3.0*a2*x*y/rP5;
      du[0][2] =  0.0;

      du[1][0] = -3.0*a2*y*x/rP5;
      du[1][1] = -3.0*a2*y*y/rP5 + a1 + a2/rP3;
      du[1][2] =  0.0;

      du[2][0] = 0.0;
      du[2][1] = 0.0;
      du[2][2] = 0.0;

      double volstrn = du[0][0]+du[1][1]+du[2][2];

      val[0] = 2.0*G*du[0][0] + L*volstrn ;
      val[1] = 2.0*G*du[1][1] + L*volstrn ;
      val[2] = 2.0*G*du[2][2] + L*volstrn ;
      val[3] = G*(du[0][1]+du[1][0]);
      val[4] = G*(du[1][2]+du[2][1]);
      val[5] = G*(du[2][0]+du[0][2]);
      */
    }


    double  strainXX(double r, double t)
    {
      return  ( a*(b1 - b2/r/r)*cos(t)*cos(t) + a*(b1 + b2/r/r)*sin(t)*sin(t) ) ;
    }

    double  strainYY(double r, double t)
    {
      return  ( a*(b1 - b2/r/r)*sin(t)*sin(t) + a*(b1 + b2/r/r)*cos(t)*cos(t) ) ;
    }

    double  strainXY(double r, double t)
    {
      return  ( -a*(b2/r/r)*sin(2.0*t) ) ;
    }

    void  strains(double r, double t, double* val)
    {
      val[0] =  a*(b1 - b2/r/r)*cos(t)*cos(t) + a*(b1 + b2/r/r)*sin(t)*sin(t) ;
      val[1] =  a*(b1 - b2/r/r)*sin(t)*sin(t) + a*(b1 + b2/r/r)*cos(t)*cos(t) ;
      val[2] = -a*(b2/r/r)*sin(2.0*t) ;

      return;
    }
};




class  ThickSphere
{
  private:
    double  pin, ri, ro, E, nu, G, K, L, A, B, fact, fact1, fact2, fact3, fact4;

  public:

    ThickSphere(double E1, double nu1)
    {
      E  = E1;
      nu = nu1;

      G = E/2.0/(1.0+nu);
      K = E/3.0/(1.0-2.0*nu);
      L = K - 2.0*G/3.0;

      ri = 100.0;
      ro = 200.0;

      pin = 1.0;

      A = pin*ri*ri*ri/(ro*ro*ro-ri*ri*ri);
      B = A*ro*ro*ro;

      fact4 = ri*ri*ri/(ro*ro*ro-ri*ri*ri);

      fact1 = pin*fact4/E;

      fact2 = (1.0-2.0*nu);
      fact3 = (1.0+nu)*ro*ro*ro/2.0;
    }

    ~ThickSphere() {}

    double  forceX(double x, double y, double z)
    { return 0.0; }

    double  forceY(double x, double y, double z)
    { return 0.0; }

    double  forceR(double r, double t, double p)
    { return 0.0; }

    double  dispR(double r, double t, double p)
    {
      return   ( fact1*( fact2*r + fact3/r/r) );
    }

    double  dispX(double r, double t, double p)
    {
      return   ( (fact1*( fact2*r + fact3/r/r) )*sin(p)*cos(t) );
    }

    double  dispY(double r, double t, double p)
    {
      return   ( (fact1*( fact2*r + fact3/r/r) )*sin(p)*sin(t) );
    }

    double  dispZ(double r, double t, double p)
    {
      return   ( (fact1*( fact2*r + fact3/r/r) )*cos(p) );
    }

    void  disp(double x, double y, double z, double* val)
    {
      double  r = sqrt(x*x + y*y + z*z);
      double  t = atan(y/x);
      double  p = acos(z/r);

      double ur = fact1*( fact2*r + fact3/r/r);

      val[0] = ur*sin(p)*cos(t);
      val[1] = ur*sin(p)*sin(t);
      val[2] = ur*cos(p);
    }

    double  pressure(double r, double t, double p)
    {
      return  A;
    }

    double  stressRR(double r, double t, double p)
    {
      return  (fact4*pin*(1.0-ro*ro*ro/r/r/r));
    }

    double  stressTT(double r, double t, double p)
    {
      return  (fact4*pin*(ro*ro*ro/r/r/r/2.0+1.0));
    }

    double  stressPP(double r, double t, double p)
    {
      return  (fact4*pin*(ro*ro*ro/r/r/r/2.0+1.0));
    }

    double  stressXX(double r, double t, double p)
    {
      return  fact*(1.0 - ro*ro*cos(2.0*t)/r/r);
    }

    double  stressYY(double r, double t, double p)
    {
       return  fact*(1.0 + ro*ro*cos(2.0*t)/r/r);
    }

    double  stressXY(double r, double t, double p)
    {
      return  fact*(- ro*ro*sin(2.0*t)/r/r);
    }

    void  stresses(double x, double y, double z, double* val)
    {
      /*
      Matrix3d  Srr, Sxx, Rot;
      Srr.setZero();
      Sxx.setZero();
      Rot.setZero();

      Srr(0,0) = fact4*pin*(ro*ro*ro/r/r/r-1.0);
      Srr(1,1) = fact4*pin*(ro*ro*ro/r/r/r/2.0+1.0);
      Srr(2,2) = fact4*pin*(ro*ro*ro/r/r/r/2.0+1.0);

      val[0] = fact4*pin*(ro*ro*ro/r/r/r-1.0);
      val[1] = fact4*pin*(ro*ro*ro/r/r/r/2.0+1.0);
      val[2] = val[1];
      val[3] = 0.0;
      val[4] = 0.0;
      val[5] = 0.0;
      */

      double  du[3][3];
      double  a1 = fact1*fact2;
      double  a2 = fact1*fact3;
      double  r = sqrt(x*x+y*y+z*z);
      double  rP3 = r*r*r;
      double  rP5 = r*r*r*r*r;

      du[0][0] = -3.0*a2*x*x/rP5 + a1 + a2/rP3;
      du[0][1] = -3.0*a2*x*y/rP5;
      du[0][2] = -3.0*a2*x*z/rP5;

      du[1][0] = -3.0*a2*y*x/rP5;
      du[1][1] = -3.0*a2*y*y/rP5 + a1 + a2/rP3;
      du[1][2] = -3.0*a2*y*z/rP5;

      du[2][0] = -3.0*a2*z*x/rP5;
      du[2][1] = -3.0*a2*z*y/rP5;
      du[2][2] = -3.0*a2*z*z/rP5 + a1 + a2/rP3;

      double volstrn = du[0][0]+du[1][1]+du[2][2];

      val[0] = 2.0*G*du[0][0] + L*volstrn ;
      val[1] = 2.0*G*du[1][1] + L*volstrn ;
      val[2] = 2.0*G*du[2][2] + L*volstrn ;
      val[3] = G*(du[0][1]+du[1][0]);
      val[4] = G*(du[1][2]+du[2][1]);
      val[5] = G*(du[2][0]+du[0][2]);

    }


    double  strainXX(double r, double t, double p)
    {
      return  0.0 ;
    }

    double  strainYY(double r, double t, double p)
    {
      return  0.0 ;
    }

    double  strainXY(double r, double t, double p)
    {
      return  0.0 ;
    }

    void  strains(double r, double t, double p, double* val)
    {
      val[0] =  0.0;
      val[1] =  0.0 ;
      val[2] =  0.0 ;

      return;
    }
};






#endif



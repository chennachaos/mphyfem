
#ifndef incl_MathBasic_h
#define incl_MathBasic_h


#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <limits.h>


static double PI    = acos(-1.0);
static double twoPI = 2*acos(-1.0);

#define EPSILON       1E-10
#define NEG_EPSILON  -1E-10

#define BIGD DBL_MAX
#define SMAD DBL_MIN
#define BIGI INT_MAX
#define SMAI INT_MIN

#define _NOZ +BIGD
#define _NOW -BIGD




double inline myAcos(double x)
{
  double tol = 1.e-12;

  if      (x > 1. - tol) return 0.;
  else if (x < tol - 1.) return PI;
  else                   return acos(x);
}





double inline power(double x, int n)
{
  double pwr = 1.;

  int i;

  for (i=0; i<n; i++) pwr *= x;
	
  return pwr;
}








template<typename Type> inline void simpleSwap(Type *x, int a, int b)
{
  if (x == NULL) return;

  Type tmp = x[a];
  
  x[a] = x[b];
  x[b] = tmp;

  return;
}




template<typename Type> inline Type norm(Type *x, int dim)
{
  Type y = (Type) 0;

  for (int i=0; i<dim; i++) y += x[i] * x[i];
	
  return (Type) sqrt(y);
}




template<typename Type> Type inline dot3(Type *x1, Type *x2)
{
  return x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2];
}




template<typename Type> Type inline dot(Type *v1, Type *v2, int n)
{
  double p = 0.;

  for (int i=0; i<n; i++) p += *(++v1) * *(++v2);

  return p;
}




template<typename Type> void inline cross3(Type *c, Type *x1, Type *x2)
{
  c[0] = x1[1]*x2[2]-x1[2]*x2[1];
  c[1] = x1[2]*x2[0]-x1[0]*x2[2];
  c[2] = x1[0]*x2[1]-x1[1]*x2[0];

  return;
}








#endif


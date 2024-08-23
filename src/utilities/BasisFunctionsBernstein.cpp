
#include "BasisFunctionsBernstein.h"
#include "headersEigen.h"
#include "BasisFunctionsLagrange.h"

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <vector>

using std::vector;



void Bernstein_BasisFuns1D(int npElem, double xi, double* N)
{
    switch(npElem)
    {
        case 1:
  
            N[0] = 1.0;
  
        break;
  
        case 2:
  
            N[0] = 0.5*(1.0 - xi);
            N[1] = 0.5*(1.0 + xi);
  
        break;
  
        case 3:
  
            // In Gauss point npElemarameter space [-1, 1]
  
            N[0] = 0.25*(1.0-xi)*(1.0-xi);
            N[1] = 0.25*(1.0+xi)*(1.0+xi);
            N[2] = 0.50*(1.0-xi*xi);
  
            // In Bezier element parameter space [0, 1]
  //           N[0] = (1.0-xi)*(1.0-xi);
  //           N[1] =  xi*xi;
  //           N[2] =  2.0*xi*(1.0-xi);
  
        break;
  
        default:
  
            printf("Bernstein_BasisFuns1D ...1... basis functions not defined for %5d nodes \n", npElem);
  
        break;
    }
    return;
}





void Bernstein_BasisFuns1D(int npElem, double xi, double* N, double* dN_dxi)
{
    switch(npElem)
    {
        case 1:
  
            N[0] = 1.0;
  
            dN_dxi[0] = 0.0;
  
        break;
  
        case 2:
  
            N[0] = 0.5*(1.0 - xi);
            N[1] = 0.5*(1.0 + xi);
  
            dN_dxi[0] = -0.5;
            dN_dxi[1] =  0.5;
  
        break;
  
        case 3:
  
            // In Gauss point npElemarameter space [-1, 1]
  
            N[0] = 0.25*(1.0-xi)*(1.0-xi);
            N[1] = 0.25*(1.0+xi)*(1.0+xi);
            N[2] = 0.50*(1.0-xi*xi);
  
            dN_dxi[0] =  0.5*(xi-1.0);
            dN_dxi[1] =  0.5*(xi+1.0);
            dN_dxi[2] = -xi;
  
            // In Bezier element parameter space [0, 1]
  //           N[0] = (1.0-xi)*(1.0-xi);
  //           N[1] =  xi*xi;
  //           N[2] =  2.0*xi*(1.0-xi);
  //
  //           dN_dxi[0] = -2.0*(1.0-xi);
  //           dN_dxi[1] =  2.0*xi;
  //           dN_dxi[2] =  2.0*(1.0-2.0*xi);
  
        break;
  
        default:
  
            printf("Bernstein_BasisFuns1D ...2... basis functions not defined for %5d nodes \n", npElem);

        break;
    }

   return;
}




void Bernstein_BasisFuns1D(int npElem, double xi, double* N, double* dN_dxi, double* d2N_dxi2)
{
    switch(npElem)
    {
        case 1:
  
            N[0] = 1.0;
  
            dN_dxi[0] = 0.0;
            d2N_dxi2[0] = 0.0;
  
        break;
  
        case 2:
            // In Gauss point npElemarameter space [-1, 1]
  
            N[0] = 0.5*(1.0 - xi);
            N[1] = 0.5*(1.0 + xi);
  
            dN_dxi[0] = -0.5;
            dN_dxi[1] =  0.5;
  
            d2N_dxi2[0] = 0.0;
            d2N_dxi2[1] = 0.0;
  
        break;
  
        case 3:
  
            // In Gauss point npElemarameter space [-1, 1]
  
            N[0] = 0.25*(1.0-xi)*(1.0-xi);
            N[1] = 0.25*(1.0+xi)*(1.0+xi);
            N[2] = 0.50*(1.0-xi*xi);
  
            dN_dxi[0] =  0.5*(xi-1.0);
            dN_dxi[1] =  0.5*(xi+1.0);
            dN_dxi[2] = -xi;
  
            d2N_dxi2[0] =  0.5;
            d2N_dxi2[1] =  0.5;
            d2N_dxi2[2] = -1.0;
  
  
            // In Bezier element parameter space [0, 1]
  //           N[0] = (1.0-xi)*(1.0-xi);
  //           N[1] =  xi*xi;
  //           N[2] =  2.0*xi*(1.0-xi);
  //
  //           dN_dxi[0] = -2.0*(1.0-xi);
  //           dN_dxi[1] =  2.0*xi;
  //           dN_dxi[2] =  2.0*(1.0-2.0*xi);
  //
  //           d2N_dxi2[0] =  2.0;
  //           d2N_dxi2[1] =  2.0;
  //           d2N_dxi2[2] = -4.0;
  
        break;
  
        default:
  
            printf("Bernstein_BasisFuns1D ...3... basis functions not defined for %5d nodes \n", npElem);

        break;
    }
   return;
}



void BernsteinBasisFunsTria(int npElem, double xi1, double xi2, double* N)
{
    double  xi3 = 1.0 - xi1 - xi2;
  
    switch(npElem)
    {
        case 1:
  
            N[0] = 1.0;
  
        break;
  
        case 3:
  
            N[0] = xi3;
            N[1] = xi1;
            N[2] = xi2;
  
        break;
  
        case 6:
  
            N[0] = xi3*xi3;
            N[1] = xi1*xi1;
            N[2] = xi2*xi2;
            N[3] = 2.0*xi1*xi3;
            N[4] = 2.0*xi1*xi2;
            N[5] = 2.0*xi2*xi3;
  
        break;
  
        default:
  
            printf("BernsteinBasisFunsTria ...1... basis functions not defined for %5d nodes \n", npElem);

        break;
    }
    return;
}




void BernsteinBasisFunsTria(int npElem, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
{
    double  xi3 = 1.0 - xi1 - xi2;

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

            dN_dxi1[0] = 0.0;
            dN_dxi2[0] = 0.0;

        break;

        case 3:

            N[0] = xi3;
            N[1] = xi1;
            N[2] = xi2;

            dN_dxi1[0] = -1.0;
            dN_dxi1[1] =  1.0;
            dN_dxi1[2] =  0.0;

            dN_dxi2[0] = -1.0;
            dN_dxi2[1] =  0.0;
            dN_dxi2[2] =  1.0;

        break;

        case 6:

            N[0] = xi3*xi3;
            N[1] = xi1*xi1;
            N[2] = xi2*xi2;
            N[3] = 2.0*xi1*xi3;
            N[4] = 2.0*xi1*xi2;
            N[5] = 2.0*xi2*xi3;

            dN_dxi1[0] = -2.0*xi3;
            dN_dxi1[1] =  2.0*xi1;
            dN_dxi1[2] =  0.0;
            dN_dxi1[3] =  2.0*(xi3 - xi1);
            dN_dxi1[4] =  2.0*xi2;
            dN_dxi1[5] = -2.0*xi2;

            dN_dxi2[0] = -2.0*xi3;
            dN_dxi2[1] =  0.0;
            dN_dxi2[2] =  2.0*xi2;
            dN_dxi2[3] = -2.0*xi1;
            dN_dxi2[4] =  2.0*xi1;
            dN_dxi2[5] =  2.0*(xi3 - xi2);

        break;

        default:

            printf("BernsteinBasisFunsTria ...2... basis functions not defined for %5d nodes \n", npElem);

        break;
    }

    return;
}



void BernsteinBasisFunsTria(int npElem, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2, double* d2N_dxi1dxi1, double* d2N_dxi2dxi2, double* d2N_dxi1dxi2)
{
  double  xi3 = 1.0 - xi1 - xi2;

  switch(npElem)
  {
      case 1:

          N[0] = 1.0;

          dN_dxi1[0] = 0.0;
          dN_dxi2[0] = 0.0;

          d2N_dxi1dxi1[0] = 0.0;
          d2N_dxi2dxi2[0] = 0.0;
          d2N_dxi1dxi2[0] = 0.0;

      break;

      case 3:

          N[0] = xi3;
          N[1] = xi1;
          N[2] = xi2;

          dN_dxi1[0] = -1.0;
          dN_dxi1[1] =  1.0;
          dN_dxi1[2] =  0.0;

          dN_dxi2[0] = -1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  1.0;

          d2N_dxi1dxi1[0] = 0.0;
          d2N_dxi1dxi1[1] = 0.0;
          d2N_dxi1dxi1[2] = 0.0;

          d2N_dxi2dxi2[0] = 0.0;
          d2N_dxi2dxi2[1] = 0.0;
          d2N_dxi2dxi2[2] = 0.0;

          d2N_dxi1dxi2[0] = 0.0;
          d2N_dxi1dxi2[1] = 0.0;
          d2N_dxi1dxi2[2] = 0.0;


      break;

      case 6:

          N[0] = xi3*xi3;
          N[1] = xi1*xi1;
          N[2] = xi2*xi2;
          N[3] = 2.0*xi1*xi3;
          N[4] = 2.0*xi1*xi2;
          N[5] = 2.0*xi2*xi3;

          dN_dxi1[0] = -2.0*xi3;
          dN_dxi1[1] =  2.0*xi1;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  2.0*(xi3 - xi1);
          dN_dxi1[4] =  2.0*xi2;
          dN_dxi1[5] = -2.0*xi2;

          dN_dxi2[0] = -2.0*xi3;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  2.0*xi2;
          dN_dxi2[3] = -2.0*xi1;
          dN_dxi2[4] =  2.0*xi1;
          dN_dxi2[5] =  2.0*(xi3 - xi2);

          d2N_dxi1dxi1[0] =  2.0;
          d2N_dxi1dxi1[1] =  2.0;
          d2N_dxi1dxi1[2] =  0.0;
          d2N_dxi1dxi1[3] = -4.0;
          d2N_dxi1dxi1[4] =  0.0;
          d2N_dxi1dxi1[5] =  0.0;

          d2N_dxi2dxi2[0] =  2.0;
          d2N_dxi2dxi2[1] =  0.0;
          d2N_dxi2dxi2[2] =  2.0;
          d2N_dxi2dxi2[3] =  0.0;
          d2N_dxi2dxi2[4] =  0.0;
          d2N_dxi2dxi2[5] = -4.0;

          d2N_dxi1dxi2[0] =  2.0;
          d2N_dxi1dxi2[1] =  0.0;
          d2N_dxi1dxi2[2] =  0.0;
          d2N_dxi1dxi2[3] = -2.0;
          d2N_dxi1dxi2[4] =  2.0;
          d2N_dxi1dxi2[5] = -2.0;


      break;

      default:

          printf("BernsteinBasisFunsTria ...3... basis functions not defined for %5d nodes \n", npElem);

      break;
  }
  return;
}





void BernsteinBasisFunsQuad(int npElem, double xi1, double xi2, double* N)
{
    double  u1, u2, u3, v1, v2, v3;
    double  du1, du2, du3, dv1, dv2, dv3;

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

        break;

        case 4:

            /*
             v2          u1,v2  u2,v2
             v1          u1,v1  u2,v1
               u1 u2
            */

            v1 = 0.5*(1.0 - xi2);
            v2 = 0.5*(1.0 + xi2);

            u1 = 0.5*(1.0 - xi1);
            u2 = 0.5*(1.0 + xi1);

            N[0] = u1*v1;
            N[1] = u2*v1;
            N[2] = u2*v2;
            N[3] = u1*v2;

        break;

        case 9:

            /*
             v2          u1,v2  u3,v2  u2,v2
             v3          u1,v3  u3,v3  u2,v3
             v1          u1,v1  u3,v1  u2,v1
               u1 u3 u2
            */

            v1 = 0.25*(1.0-xi2)*(1.0-xi2);
            v2 = 0.25*(1.0+xi2)*(1.0+xi2);
            v3 = 0.50*(1.0-xi2*xi2);

            u1 = 0.25*(1.0-xi1)*(1.0-xi1);
            u2 = 0.25*(1.0+xi1)*(1.0+xi1);
            u3 = 0.50*(1.0-xi1*xi1);

            N[0] = u1*v1;
            N[1] = u2*v1;
            N[2] = u2*v2;
            N[3] = u1*v2;

            N[4] = u3*v1;
            N[5] = u2*v3;
            N[6] = u3*v2;
            N[7] = u1*v3;

            N[8] = u3*v3;

        break;

        default:

          printf("BernsteinBasisFunsQuad ...1... basis functions not defined for %5d nodes \n", npElem);

        break;
    }

    return;
}




void BernsteinBasisFunsQuad(int npElem, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
{
    double  u1, u2, u3, v1, v2, v3;
    double  du1, du2, du3, dv1, dv2, dv3;

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

            dN_dxi1[0] = 0.0;
            dN_dxi2[0] = 0.0;

        break;

        case 4:

            /*
             v2          u1,v2  u2,v2
             v1          u1,v1  u2,v1
               u1 u2
            */

            v1 = 0.5*(1.0 - xi2);
            v2 = 0.5*(1.0 + xi2);

            dv1 = -0.5;
            dv2 =  0.5;

            u1 = 0.5*(1.0 - xi1);
            u2 = 0.5*(1.0 + xi1);

            du1 = -0.5;
            du2 =  0.5;

            N[0] = u1*v1;
            N[1] = u2*v1;
            N[2] = u2*v2;
            N[3] = u1*v2;

            dN_dxi1[0] = du1*v1;
            dN_dxi1[1] = du2*v1;
            dN_dxi1[2] = du2*v2;
            dN_dxi1[3] = du1*v2;

            dN_dxi2[0] = u1*dv1;
            dN_dxi2[1] = u2*dv1;
            dN_dxi2[2] = u2*dv2;
            dN_dxi2[3] = u1*dv2;

        break;

        case 9:

            /*
             v2          u1,v2  u3,v2  u2,v2
             v3          u1,v3  u3,v3  u2,v3
             v1          u1,v1  u3,v1  u2,v1
               u1 u3 u2
            */

            v1 = 0.25*(1.0-xi2)*(1.0-xi2);
            v2 = 0.25*(1.0+xi2)*(1.0+xi2);
            v3 = 0.50*(1.0-xi2*xi2);

            dv1 =  0.5*(xi2-1.0);
            dv2 =  0.5*(xi2+1.0);
            dv3 = -xi2;

            u1 = 0.25*(1.0-xi1)*(1.0-xi1);
            u2 = 0.25*(1.0+xi1)*(1.0+xi1);
            u3 = 0.50*(1.0-xi1*xi1);

            du1 =  0.5*(xi1-1.0);
            du2 =  0.5*(xi1+1.0);
            du3 = -xi1;

            N[0] = u1*v1;
            N[1] = u2*v1;
            N[2] = u2*v2;
            N[3] = u1*v2;
            N[4] = u3*v1;
            N[5] = u2*v3;
            N[6] = u3*v2;
            N[7] = u1*v3;
            N[8] = u3*v3;

            dN_dxi1[0] = du1*v1;
            dN_dxi1[1] = du2*v1;
            dN_dxi1[2] = du2*v2;
            dN_dxi1[3] = du1*v2;
            dN_dxi1[4] = du3*v1;
            dN_dxi1[5] = du2*v3;
            dN_dxi1[6] = du3*v2;
            dN_dxi1[7] = du1*v3;
            dN_dxi1[8] = du3*v3;

            dN_dxi2[0] = u1*dv1;
            dN_dxi2[1] = u2*dv1;
            dN_dxi2[2] = u2*dv2;
            dN_dxi2[3] = u1*dv2;
            dN_dxi2[4] = u3*dv1;
            dN_dxi2[5] = u2*dv3;
            dN_dxi2[6] = u3*dv2;
            dN_dxi2[7] = u1*dv3;
            dN_dxi2[8] = u3*dv3;

        break;

        default:

          printf("BernsteinBasisFunsQuad ...2... basis functions not defined for %5d nodes \n", npElem);

        break;
    }

    return;
}





void BernsteinBasisFunsTetra(int npElem, double xi1, double xi2, double xi3, double* N)
{
    double  xi4 = 1.0 - xi1 - xi2 - xi3;

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

        break;

        case 4:

            N[0] = xi4;
            N[1] = xi1;
            N[2] = xi2;
            N[3] = xi3;

        break;

        case 10:

            N[0] = xi4*xi4;
            N[1] = xi1*xi1;
            N[2] = xi2*xi2;
            N[3] = xi3*xi3;
            N[4] = 2.0*xi1*xi4;
            N[5] = 2.0*xi1*xi2;
            N[6] = 2.0*xi2*xi4;
            N[7] = 2.0*xi4*xi3;
            N[8] = 2.0*xi1*xi3;
            N[9] = 2.0*xi2*xi3;

        break;

        default:

          printf("BernsteinBasisFunsTetra ...1... basis functions not defined for %5d nodes \n", npElem);

        break;
    }
    return;
}




void BernsteinBasisFunsTetra(int npElem, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
    double  xi4 = 1.0 - xi1 - xi2 - xi3;

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

            dN_dxi1[0] = 0.0;
            dN_dxi2[0] = 0.0;
            dN_dxi3[0] = 0.0;

        break;

        case 4:

            N[0] = xi4;
            N[1] = xi1;
            N[2] = xi2;
            N[3] = xi3;

            dN_dxi1[0] = -1.0;
            dN_dxi1[1] =  1.0;
            dN_dxi1[2] =  0.0;
            dN_dxi1[3] =  0.0;

            dN_dxi2[0] = -1.0;
            dN_dxi2[1] =  0.0;
            dN_dxi2[2] =  1.0;
            dN_dxi2[3] =  0.0;

            dN_dxi3[0] = -1.0;
            dN_dxi3[1] =  0.0;
            dN_dxi3[2] =  0.0;
            dN_dxi3[3] =  1.0;

        break;

        case 10:

            N[0] = xi4*xi4;
            N[1] = xi1*xi1;
            N[2] = xi2*xi2;
            N[3] = xi3*xi3;
            N[4] = 2.0*xi1*xi4;
            N[5] = 2.0*xi1*xi2;
            N[6] = 2.0*xi2*xi4;
            N[7] = 2.0*xi4*xi3;
            N[8] = 2.0*xi1*xi3;
            N[9] = 2.0*xi2*xi3;

            dN_dxi1[0] = -2.0*xi4;
            dN_dxi1[1] =  2.0*xi1;
            dN_dxi1[2] =  0.0;
            dN_dxi1[3] =  0.0;
            dN_dxi1[4] =  2.0*(xi4-xi1);
            dN_dxi1[5] =  2.0*xi2;
            dN_dxi1[6] = -2.0*xi2;
            dN_dxi1[7] = -2.0*xi3;
            dN_dxi1[8] =  2.0*xi3;
            dN_dxi1[9] =  0.0;

            dN_dxi2[0] = -2.0*xi4;
            dN_dxi2[1] =  0.0;
            dN_dxi2[2] =  2.0*xi2;
            dN_dxi2[3] =  0.0;
            dN_dxi2[4] = -2.0*xi1;
            dN_dxi2[5] =  2.0*xi1;
            dN_dxi2[6] =  2.0*(xi4-xi2);
            dN_dxi2[7] = -2.0*xi3;
            dN_dxi2[8] =  0.0;
            dN_dxi2[9] =  2.0*xi3;

            dN_dxi3[0] = -2.0*xi4;
            dN_dxi3[1] =  0.0;
            dN_dxi3[2] =  0.0;
            dN_dxi3[3] =  2.0*xi3;
            dN_dxi3[4] = -2.0*xi1;
            dN_dxi3[5] =  0.0;
            dN_dxi3[6] = -2.0*xi2;
            dN_dxi3[7] =  2.0*(xi4-xi3);
            dN_dxi3[8] =  2.0*xi1;
            dN_dxi3[9] =  2.0*xi2;

        break;

        default:

          printf("BernsteinBasisFunsTetra ...2... basis functions not defined for %5d nodes \n", npElem);

        break;
    }
    return;
}





void BernsteinBasisFunsHexa(int npElem, double xi1, double xi2, double xi3, double* N)
{
    double  Nu[3], Nv[3], Nw[3];

    Bernstein_BasisFuns1D(npElem, xi1, Nu);
    Bernstein_BasisFuns1D(npElem, xi2, Nv);
    Bernstein_BasisFuns1D(npElem, xi3, Nw);

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

        break;

        case 8:

            N[0] = Nu[0]*Nv[0]*Nw[0];
            N[1] = Nu[1]*Nv[0]*Nw[0];
            N[2] = Nu[1]*Nv[1]*Nw[0];
            N[3] = Nu[0]*Nv[1]*Nw[0];
            N[4] = Nu[0]*Nv[0]*Nw[1];
            N[5] = Nu[1]*Nv[0]*Nw[1];
            N[6] = Nu[1]*Nv[1]*Nw[1];
            N[7] = Nu[0]*Nv[1]*Nw[1];

        break;

        case 27:

            N[0]  = Nu[0]*Nv[0]*Nw[0];
            N[1]  = Nu[1]*Nv[0]*Nw[0];
            N[2]  = Nu[1]*Nv[1]*Nw[0];
            N[3]  = Nu[0]*Nv[1]*Nw[0];
            N[8]  = Nu[2]*Nv[0]*Nw[0];
            N[11] = Nu[1]*Nv[2]*Nw[0];
            N[13] = Nu[2]*Nv[1]*Nw[0];
            N[9]  = Nu[0]*Nv[2]*Nw[0];
            N[20] = Nu[2]*Nv[2]*Nw[0];
            N[10] = Nu[0]*Nv[0]*Nw[2];
            N[12] = Nu[1]*Nv[0]*Nw[2];
            N[14] = Nu[1]*Nv[1]*Nw[2];
            N[15] = Nu[0]*Nv[1]*Nw[2];
            N[21] = Nu[2]*Nv[0]*Nw[2];
            N[23] = Nu[1]*Nv[2]*Nw[2];
            N[24] = Nu[2]*Nv[1]*Nw[2];
            N[22] = Nu[0]*Nv[2]*Nw[2];
            N[26] = Nu[2]*Nv[2]*Nw[2];
            N[4]  = Nu[0]*Nv[0]*Nw[1];
            N[5]  = Nu[1]*Nv[0]*Nw[1];
            N[6]  = Nu[1]*Nv[1]*Nw[1];
            N[7]  = Nu[0]*Nv[1]*Nw[1];
            N[16] = Nu[2]*Nv[0]*Nw[1];
            N[18] = Nu[1]*Nv[2]*Nw[1];
            N[19] = Nu[2]*Nv[1]*Nw[1];
            N[17] = Nu[0]*Nv[2]*Nw[1];
            N[25] = Nu[2]*Nv[2]*Nw[1];

        break;

        default:

          printf("BernsteinBasisFunsHexa ...1... basis functions not defined for %5d nodes \n", npElem);

        break;
    }

    return;
}




void BernsteinBasisFunsHexa(int npElem, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
/*
       v
3----------2            3----13----2           3----13----2
|\     ^   |\           |\         |\          |\         |\
| \    |   | \          | 15       | 14        |15    24  | 14
|  \   |   |  \         9  \       11 \        9  \ 20    11 \
|   7------+---6        |   7----19+---6       |   7----19+---6
|   |  +-- |-- | -> u   |   |      |   |       |22 |  26  | 23|
0---+---\--1   |        0---+-8----1   |       0---+-8----1   |
 \  |    \  \  |         \  17      \  18       \ 17    25 \  18
  \ |     \  \ |         10 |        12|        10 |  21    12|
   \|      w  \|           \|         \|          \|         \|
    4----------5            4----16----5           4----16----5
*/

    int  ind = npElem+1;
    double  Nu[ind], dNu[ind], Nv[ind], dNv[ind], Nw[ind], dNw[ind];

    Bernstein_BasisFuns1D(npElem, xi1, Nu, dNu);
    Bernstein_BasisFuns1D(npElem, xi2, Nv, dNv);
    Bernstein_BasisFuns1D(npElem, xi3, Nw, dNw);

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

            dN_dxi1[0] = 0.0;
            dN_dxi2[0] = 0.0;
            dN_dxi3[0] = 0.0;

        break;

        case 8:

            N[0] = Nu[0]*Nv[0]*Nw[0];
            N[1] = Nu[1]*Nv[0]*Nw[0];
            N[2] = Nu[1]*Nv[1]*Nw[0];
            N[3] = Nu[0]*Nv[1]*Nw[0];
            N[4] = Nu[0]*Nv[0]*Nw[1];
            N[5] = Nu[1]*Nv[0]*Nw[1];
            N[6] = Nu[1]*Nv[1]*Nw[1];
            N[7] = Nu[0]*Nv[1]*Nw[1];

            dN_dxi1[0] = dNu[0]*Nv[0]*Nw[0];
            dN_dxi1[1] = dNu[1]*Nv[0]*Nw[0];
            dN_dxi1[2] = dNu[1]*Nv[1]*Nw[0];
            dN_dxi1[3] = dNu[0]*Nv[1]*Nw[0];
            dN_dxi1[4] = dNu[0]*Nv[0]*Nw[1];
            dN_dxi1[5] = dNu[1]*Nv[0]*Nw[1];
            dN_dxi1[6] = dNu[1]*Nv[1]*Nw[1];
            dN_dxi1[7] = dNu[0]*Nv[1]*Nw[1];

            dN_dxi2[0] = Nu[0]*dNv[0]*Nw[0];
            dN_dxi2[1] = Nu[1]*dNv[0]*Nw[0];
            dN_dxi2[2] = Nu[1]*dNv[1]*Nw[0];
            dN_dxi2[3] = Nu[0]*dNv[1]*Nw[0];
            dN_dxi2[4] = Nu[0]*dNv[0]*Nw[1];
            dN_dxi2[5] = Nu[1]*dNv[0]*Nw[1];
            dN_dxi2[6] = Nu[1]*dNv[1]*Nw[1];
            dN_dxi2[7] = Nu[0]*dNv[1]*Nw[1];

            dN_dxi3[0] = Nu[0]*Nv[0]*dNw[0];
            dN_dxi3[1] = Nu[1]*Nv[0]*dNw[0];
            dN_dxi3[2] = Nu[1]*Nv[1]*dNw[0];
            dN_dxi3[3] = Nu[0]*Nv[1]*dNw[0];
            dN_dxi3[4] = Nu[0]*Nv[0]*dNw[1];
            dN_dxi3[5] = Nu[1]*Nv[0]*dNw[1];
            dN_dxi3[6] = Nu[1]*Nv[1]*dNw[1];
            dN_dxi3[7] = Nu[0]*Nv[1]*dNw[1];

        break;

        case 27:

            N[0]  = Nu[0]*Nv[0]*Nw[0];
            N[1]  = Nu[1]*Nv[0]*Nw[0];
            N[2]  = Nu[1]*Nv[1]*Nw[0];
            N[3]  = Nu[0]*Nv[1]*Nw[0];
            N[8]  = Nu[2]*Nv[0]*Nw[0];
            N[11] = Nu[1]*Nv[2]*Nw[0];
            N[13] = Nu[2]*Nv[1]*Nw[0];
            N[9]  = Nu[0]*Nv[2]*Nw[0];
            N[20] = Nu[2]*Nv[2]*Nw[0];
            N[4]  = Nu[0]*Nv[0]*Nw[1];
            N[5]  = Nu[1]*Nv[0]*Nw[1];
            N[6]  = Nu[1]*Nv[1]*Nw[1];
            N[7]  = Nu[0]*Nv[1]*Nw[1];
            N[16] = Nu[2]*Nv[0]*Nw[1];
            N[18] = Nu[1]*Nv[2]*Nw[1];
            N[19] = Nu[2]*Nv[1]*Nw[1];
            N[17] = Nu[0]*Nv[2]*Nw[1];
            N[25] = Nu[2]*Nv[2]*Nw[1];
            N[10] = Nu[0]*Nv[0]*Nw[2];
            N[12] = Nu[1]*Nv[0]*Nw[2];
            N[14] = Nu[1]*Nv[1]*Nw[2];
            N[15] = Nu[0]*Nv[1]*Nw[2];
            N[21] = Nu[2]*Nv[0]*Nw[2];
            N[23] = Nu[1]*Nv[2]*Nw[2];
            N[24] = Nu[2]*Nv[1]*Nw[2];
            N[22] = Nu[0]*Nv[2]*Nw[2];
            N[26] = Nu[2]*Nv[2]*Nw[2];

            dN_dxi1[0]  = dNu[0]*Nv[0]*Nw[0];
            dN_dxi1[1]  = dNu[1]*Nv[0]*Nw[0];
            dN_dxi1[2]  = dNu[1]*Nv[1]*Nw[0];
            dN_dxi1[3]  = dNu[0]*Nv[1]*Nw[0];
            dN_dxi1[8]  = dNu[2]*Nv[0]*Nw[0];
            dN_dxi1[11] = dNu[1]*Nv[2]*Nw[0];
            dN_dxi1[13] = dNu[2]*Nv[1]*Nw[0];
            dN_dxi1[9]  = dNu[0]*Nv[2]*Nw[0];
            dN_dxi1[20] = dNu[2]*Nv[2]*Nw[0];
            dN_dxi1[4]  = dNu[0]*Nv[0]*Nw[1];
            dN_dxi1[5]  = dNu[1]*Nv[0]*Nw[1];
            dN_dxi1[6]  = dNu[1]*Nv[1]*Nw[1];
            dN_dxi1[7]  = dNu[0]*Nv[1]*Nw[1];
            dN_dxi1[16] = dNu[2]*Nv[0]*Nw[1];
            dN_dxi1[18] = dNu[1]*Nv[2]*Nw[1];
            dN_dxi1[19] = dNu[2]*Nv[1]*Nw[1];
            dN_dxi1[17] = dNu[0]*Nv[2]*Nw[1];
            dN_dxi1[25] = dNu[2]*Nv[2]*Nw[1];
            dN_dxi1[10] = dNu[0]*Nv[0]*Nw[2];
            dN_dxi1[12] = dNu[1]*Nv[0]*Nw[2];
            dN_dxi1[14] = dNu[1]*Nv[1]*Nw[2];
            dN_dxi1[15] = dNu[0]*Nv[1]*Nw[2];
            dN_dxi1[21] = dNu[2]*Nv[0]*Nw[2];
            dN_dxi1[23] = dNu[1]*Nv[2]*Nw[2];
            dN_dxi1[24] = dNu[2]*Nv[1]*Nw[2];
            dN_dxi1[22] = dNu[0]*Nv[2]*Nw[2];
            dN_dxi1[26] = dNu[2]*Nv[2]*Nw[2];

            dN_dxi2[0]  = Nu[0]*dNv[0]*Nw[0];
            dN_dxi2[1]  = Nu[1]*dNv[0]*Nw[0];
            dN_dxi2[2]  = Nu[1]*dNv[1]*Nw[0];
            dN_dxi2[3]  = Nu[0]*dNv[1]*Nw[0];
            dN_dxi2[8]  = Nu[2]*dNv[0]*Nw[0];
            dN_dxi2[11] = Nu[1]*dNv[2]*Nw[0];
            dN_dxi2[13] = Nu[2]*dNv[1]*Nw[0];
            dN_dxi2[9]  = Nu[0]*dNv[2]*Nw[0];
            dN_dxi2[20] = Nu[2]*dNv[2]*Nw[0];
            dN_dxi2[4]  = Nu[0]*dNv[0]*Nw[1];
            dN_dxi2[5]  = Nu[1]*dNv[0]*Nw[1];
            dN_dxi2[6]  = Nu[1]*dNv[1]*Nw[1];
            dN_dxi2[7]  = Nu[0]*dNv[1]*Nw[1];
            dN_dxi2[16] = Nu[2]*dNv[0]*Nw[1];
            dN_dxi2[18] = Nu[1]*dNv[2]*Nw[1];
            dN_dxi2[19] = Nu[2]*dNv[1]*Nw[1];
            dN_dxi2[17] = Nu[0]*dNv[2]*Nw[1];
            dN_dxi2[25] = Nu[2]*dNv[2]*Nw[1];
            dN_dxi2[10] = Nu[0]*dNv[0]*Nw[2];
            dN_dxi2[12] = Nu[1]*dNv[0]*Nw[2];
            dN_dxi2[14] = Nu[1]*dNv[1]*Nw[2];
            dN_dxi2[15] = Nu[0]*dNv[1]*Nw[2];
            dN_dxi2[21] = Nu[2]*dNv[0]*Nw[2];
            dN_dxi2[23] = Nu[1]*dNv[2]*Nw[2];
            dN_dxi2[24] = Nu[2]*dNv[1]*Nw[2];
            dN_dxi2[22] = Nu[0]*dNv[2]*Nw[2];
            dN_dxi2[26] = Nu[2]*dNv[2]*Nw[2];

            dN_dxi3[0]  = Nu[0]*Nv[0]*dNw[0];
            dN_dxi3[1]  = Nu[1]*Nv[0]*dNw[0];
            dN_dxi3[2]  = Nu[1]*Nv[1]*dNw[0];
            dN_dxi3[3]  = Nu[0]*Nv[1]*dNw[0];
            dN_dxi3[8]  = Nu[2]*Nv[0]*dNw[0];
            dN_dxi3[11] = Nu[1]*Nv[2]*dNw[0];
            dN_dxi3[13] = Nu[2]*Nv[1]*dNw[0];
            dN_dxi3[9]  = Nu[0]*Nv[2]*dNw[0];
            dN_dxi3[20] = Nu[2]*Nv[2]*dNw[0];
            dN_dxi3[4]  = Nu[0]*Nv[0]*dNw[1];
            dN_dxi3[5]  = Nu[1]*Nv[0]*dNw[1];
            dN_dxi3[6]  = Nu[1]*Nv[1]*dNw[1];
            dN_dxi3[7]  = Nu[0]*Nv[1]*dNw[1];
            dN_dxi3[16] = Nu[2]*Nv[0]*dNw[1];
            dN_dxi3[18] = Nu[1]*Nv[2]*dNw[1];
            dN_dxi3[19] = Nu[2]*Nv[1]*dNw[1];
            dN_dxi3[17] = Nu[0]*Nv[2]*dNw[1];
            dN_dxi3[25] = Nu[2]*Nv[2]*dNw[1];
            dN_dxi3[10] = Nu[0]*Nv[0]*dNw[2];
            dN_dxi3[12] = Nu[1]*Nv[0]*dNw[2];
            dN_dxi3[14] = Nu[1]*Nv[1]*dNw[2];
            dN_dxi3[15] = Nu[0]*Nv[1]*dNw[2];
            dN_dxi3[21] = Nu[2]*Nv[0]*dNw[2];
            dN_dxi3[23] = Nu[1]*Nv[2]*dNw[2];
            dN_dxi3[24] = Nu[2]*Nv[1]*dNw[2];
            dN_dxi3[22] = Nu[0]*Nv[2]*dNw[2];
            dN_dxi3[26] = Nu[2]*Nv[2]*dNw[2];

        break;

        default:

          printf("BernsteinBasisFunsHexa ...2... basis functions not defined for %5d nodes \n", npElem);

        break;
    }

    return;
}





void BernsteinBasisFunsWedge(int npElem, double xi1, double xi2, double xi4, double* N)
{
    // xi1 and xi2 are the parametric coordinates for the triangle, and
    // xi4 is the parametric coordinate in the direction normal to the triangle
    double  Nt[6], Nq[3];

    switch(npElem)
    {
        case 6:

            BernsteinBasisFunsTria(npElem, xi1, xi2, Nt);
            Bernstein_BasisFuns1D(npElem, xi4, Nq);

            N[0] = Nt[0]*Nq[0];
            N[1] = Nt[1]*Nq[0];
            N[2] = Nt[2]*Nq[0];
            N[3] = Nt[0]*Nq[1];
            N[4] = Nt[1]*Nq[1];
            N[5] = Nt[2]*Nq[1];

        break;

        case 18:

            BernsteinBasisFunsTria(npElem, xi1, xi2, Nt);
            Bernstein_BasisFuns1D(npElem, xi4, Nq);

            N[0]  = Nt[0]*Nq[0];
            N[1]  = Nt[1]*Nq[0];
            N[2]  = Nt[2]*Nq[0];
            N[6]  = Nt[3]*Nq[0];
            N[9]  = Nt[4]*Nq[0];
            N[7]  = Nt[5]*Nq[0];
            N[8]  = Nt[0]*Nq[2];
            N[10] = Nt[1]*Nq[2];
            N[11] = Nt[2]*Nq[2];
            N[15] = Nt[3]*Nq[2];
            N[17] = Nt[4]*Nq[2];
            N[16] = Nt[5]*Nq[2];
            N[3]  = Nt[0]*Nq[1];
            N[4]  = Nt[1]*Nq[1];
            N[5]  = Nt[2]*Nq[1];
            N[12] = Nt[3]*Nq[1];
            N[14] = Nt[4]*Nq[1];
            N[13] = Nt[5]*Nq[1];

        break;

        default:

          printf("BernsteinBasisFunsWedge ...1... basis functions not defined for %5d nodes \n", npElem);

        break;
    }
    return;
}





void BernsteinBasisFunsWedge(int npElem, double xi1, double xi2, double xi4, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi4)
{
    // xi1 and xi2 are the parametric coordinates for the triangle, and
    // xi4 is the parametric coordinate in the direction normal to the triangle
    double  Nt[6], dNt1[6], dNt2[6], Nq[3], dNq[3];

    switch(npElem)
    {
        case 6:

            BernsteinBasisFunsTria(npElem, xi1, xi2, Nt, dNt1, dNt2);
            Bernstein_BasisFuns1D(npElem, xi4, Nq, dNq);

            N[0] = Nt[0]*Nq[0];
            N[1] = Nt[1]*Nq[0];
            N[2] = Nt[2]*Nq[0];
            N[3] = Nt[0]*Nq[1];
            N[4] = Nt[1]*Nq[1];
            N[5] = Nt[2]*Nq[1];

            dN_dxi1[0] = dNt1[0]*Nq[0];
            dN_dxi1[1] = dNt1[1]*Nq[0];
            dN_dxi1[2] = dNt1[2]*Nq[0];
            dN_dxi1[3] = dNt1[0]*Nq[1];
            dN_dxi1[4] = dNt1[1]*Nq[1];
            dN_dxi1[5] = dNt1[2]*Nq[1];

            dN_dxi2[0] = dNt2[0]*Nq[0];
            dN_dxi2[1] = dNt2[1]*Nq[0];
            dN_dxi2[2] = dNt2[2]*Nq[0];
            dN_dxi2[3] = dNt2[0]*Nq[1];
            dN_dxi2[4] = dNt2[1]*Nq[1];
            dN_dxi2[5] = dNt2[2]*Nq[1];

            dN_dxi4[0] = Nt[0]*dNq[0];
            dN_dxi4[1] = Nt[1]*dNq[0];
            dN_dxi4[2] = Nt[2]*dNq[0];
            dN_dxi4[3] = Nt[0]*dNq[1];
            dN_dxi4[4] = Nt[1]*dNq[1];
            dN_dxi4[5] = Nt[2]*dNq[1];

        break;

        case 18:

            BernsteinBasisFunsTria(npElem, xi1, xi2, Nt, dNt1, dNt2);
            Bernstein_BasisFuns1D(npElem, xi4, Nq, dNq);

            N[0]  = Nt[0]*Nq[0];
            N[1]  = Nt[1]*Nq[0];
            N[2]  = Nt[2]*Nq[0];
            N[6]  = Nt[3]*Nq[0];
            N[9]  = Nt[4]*Nq[0];
            N[7]  = Nt[5]*Nq[0];
            N[8]  = Nt[0]*Nq[2];
            N[10] = Nt[1]*Nq[2];
            N[11] = Nt[2]*Nq[2];
            N[15] = Nt[3]*Nq[2];
            N[17] = Nt[4]*Nq[2];
            N[16] = Nt[5]*Nq[2];
            N[3]  = Nt[0]*Nq[1];
            N[4]  = Nt[1]*Nq[1];
            N[5]  = Nt[2]*Nq[1];
            N[12] = Nt[3]*Nq[1];
            N[14] = Nt[4]*Nq[1];
            N[13] = Nt[5]*Nq[1];

            dN_dxi1[0]  = dNt1[0]*Nq[0];
            dN_dxi1[1]  = dNt1[1]*Nq[0];
            dN_dxi1[2]  = dNt1[2]*Nq[0];
            dN_dxi1[6]  = dNt1[3]*Nq[0];
            dN_dxi1[9]  = dNt1[4]*Nq[0];
            dN_dxi1[7]  = dNt1[5]*Nq[0];
            dN_dxi1[8]  = dNt1[0]*Nq[2];
            dN_dxi1[10] = dNt1[1]*Nq[2];
            dN_dxi1[11] = dNt1[2]*Nq[2];
            dN_dxi1[15] = dNt1[3]*Nq[2];
            dN_dxi1[17] = dNt1[4]*Nq[2];
            dN_dxi1[16] = dNt1[5]*Nq[2];
            dN_dxi1[3]  = dNt1[0]*Nq[1];
            dN_dxi1[4]  = dNt1[1]*Nq[1];
            dN_dxi1[5]  = dNt1[2]*Nq[1];
            dN_dxi1[12] = dNt1[3]*Nq[1];
            dN_dxi1[14] = dNt1[4]*Nq[1];
            dN_dxi1[13] = dNt1[5]*Nq[1];

            dN_dxi2[0]  = dNt2[0]*Nq[0];
            dN_dxi2[1]  = dNt2[1]*Nq[0];
            dN_dxi2[2]  = dNt2[2]*Nq[0];
            dN_dxi2[6]  = dNt2[3]*Nq[0];
            dN_dxi2[9]  = dNt2[4]*Nq[0];
            dN_dxi2[7]  = dNt2[5]*Nq[0];
            dN_dxi2[8]  = dNt2[0]*Nq[2];
            dN_dxi2[10] = dNt2[1]*Nq[2];
            dN_dxi2[11] = dNt2[2]*Nq[2];
            dN_dxi2[15] = dNt2[3]*Nq[2];
            dN_dxi2[17] = dNt2[4]*Nq[2];
            dN_dxi2[16] = dNt2[5]*Nq[2];
            dN_dxi2[3]  = dNt2[0]*Nq[1];
            dN_dxi2[4]  = dNt2[1]*Nq[1];
            dN_dxi2[5]  = dNt2[2]*Nq[1];
            dN_dxi2[12] = dNt2[3]*Nq[1];
            dN_dxi2[14] = dNt2[4]*Nq[1];
            dN_dxi2[13] = dNt2[5]*Nq[1];

            dN_dxi4[0]  = Nt[0]*dNq[0];
            dN_dxi4[1]  = Nt[1]*dNq[0];
            dN_dxi4[2]  = Nt[2]*dNq[0];
            dN_dxi4[6]  = Nt[3]*dNq[0];
            dN_dxi4[9]  = Nt[4]*dNq[0];
            dN_dxi4[7]  = Nt[5]*dNq[0];
            dN_dxi4[8]  = Nt[0]*dNq[2];
            dN_dxi4[10] = Nt[1]*dNq[2];
            dN_dxi4[11] = Nt[2]*dNq[2];
            dN_dxi4[15] = Nt[3]*dNq[2];
            dN_dxi4[17] = Nt[4]*dNq[2];
            dN_dxi4[16] = Nt[5]*dNq[2];
            dN_dxi4[3]  = Nt[0]*dNq[1];
            dN_dxi4[4]  = Nt[1]*dNq[1];
            dN_dxi4[5]  = Nt[2]*dNq[1];
            dN_dxi4[12] = Nt[3]*dNq[1];
            dN_dxi4[14] = Nt[4]*dNq[1];
            dN_dxi4[13] = Nt[5]*dNq[1];

        break;

        default:

          printf("BernsteinBasisFunsWedge ...2... basis functions not defined for %5d nodes \n", npElem);

        break;
    }
    return;
}


void BernsteinBasisFunsEdge2D(int npElem, double* param, double *xNode, double* yNode, double *N, double *dN_dx, double *dN_dy, double *normal, double& curvature, double& J)
{
    double  du_dx, dx, dy, d2x, d2y;
    vector<double>  dN1(npElem), dN2(npElem);

    Bernstein_BasisFuns1D(npElem, param[0], N, &dN1[0], &dN2[0]);

    if(npElem == 0)
    {
      J = 1.0;
    }
    else
    {
      dx = dy = 0.0;
      d2x = d2y = 0.0;
      for(int ii=0; ii<npElem; ii++)
      {
        dx +=  (xNode[ii] * dN1[ii]);
        dy +=  (yNode[ii] * dN1[ii]);

        d2x +=  (xNode[ii] * dN2[ii]);
        d2y +=  (yNode[ii] * dN2[ii]);
      }
      J = sqrt(dx*dx+dy*dy);
    }

    du_dx = 1.0/J ;

    // Compute the normal
    normal[0] =  dy/J;
    normal[1] = -dx/J;

    // Compute curvature
    //double  dJ  = (dx*d2x+dy*d2y)/J;
    //double  dnx =  d2y/J - dJ*dy/J/J;
    //double  dny = -d2x/J + dJ*dx/J/J;

    curvature = (dy*d2x - dx*d2y)/J/J/J;

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(int ii=0; ii<npElem; ii++)
    {
      dN_dx[ii] = dN1[ii] / dx ;
      dN_dy[ii] = dN1[ii] / dy ;
    }

    return;
}



void BernsteinBasisFunsEdge3D(int npElem, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
    double  du_dx, dx, dy, dz ;
    vector<double>  dN1(npElem);

    Bernstein_BasisFuns1D(npElem, param[0], N, &dN1[0]);

    if(npElem == 0)
    {
        Jac = 1.0;
    }
    else
    {
        dx = dy = dz = 0.0;
        for(int ii=0; ii<npElem; ii++)
        {
            dx +=  (xNode[ii] * dN1[ii]);
            dy +=  (yNode[ii] * dN1[ii]);
            dz +=  (zNode[ii] * dN1[ii]);
        }
        Jac = sqrt(dx*dx+dy*dy+dz*dz);
    }
    du_dx = 1.0/Jac ;

    // Compute derivatives of basis functions w.r.t physical coordinates
    //for(ii=0; ii<nlbf; ii++)
      //dN_dx[ii] = dN1[ii] * du_dx ;

    return;
}


void BernsteinBasisFunsFaceTria(int npElem, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double *tangent1, double *tangent2, double& curvature, double& J)
{
    vector<double>  dN1(npElem), dN2(npElem), d2N11(npElem), d2N22(npElem), d2N12(npElem);

    BernsteinBasisFunsTria(npElem, param[0], param[1], N, &dN1[0], &dN2[0], &d2N11[0], &d2N22[0], &d2N12[0]);

    double  x1=0.0, y1=0.0, z1=0.0, x2=0.0, y2=0.0, z2=0.0;
    double  x11=0.0, y11=0.0, z11=0.0, x22=0.0, y22=0.0, z22=0.0, x21=0.0, y21=0.0, z21=0.0;

    for(int ii=0; ii<npElem; ii++)
    {
      // dx_dxi1
      x1 += (xNode[ii] * dN1[ii]);
      y1 += (yNode[ii] * dN1[ii]);
      z1 += (zNode[ii] * dN1[ii]);

      // dx_dxi2
      x2 += (xNode[ii] * dN2[ii]);
      y2 += (yNode[ii] * dN2[ii]);
      z2 += (zNode[ii] * dN2[ii]);

      // d2x_dxi1dxi1
      x11 += (xNode[ii] * d2N11[ii]);
      y11 += (yNode[ii] * d2N11[ii]);
      z11 += (zNode[ii] * d2N11[ii]);

      // d2x_dxi2dxi2
      x22 += (xNode[ii] * d2N22[ii]);
      y22 += (yNode[ii] * d2N22[ii]);
      z22 += (zNode[ii] * d2N22[ii]);

      // d2x_dxi2dxi1
      x21 += (xNode[ii] * d2N12[ii]);
      y21 += (yNode[ii] * d2N12[ii]);
      z21 += (zNode[ii] * d2N12[ii]);
    }

    double  x12 = x21, y12 = y21, z12 = z21;

    // Compute the normal

    double  mx = y1*z2 - z1*y2;
    double  my = z1*x2 - x1*z2;
    double  mz = x1*y2 - y1*x2;

    J = sqrt(mx*mx + my*my + mz*mz);

    normal[0] = mx/J;
    normal[1] = my/J;
    normal[2] = mz/J;

    // Compute the tangents
    tangent1[0] = 0.0;
    tangent1[1] = 0.0;
    tangent1[2] = 0.0;

    tangent2[0] = 0.0;
    tangent2[1] = 0.0;
    tangent2[2] = 0.0;

    // Compute the curvature

    double  dmx1 = y11*z2 + y1*z12 - z11*y2 - z1*y12;
    double  dmx2 = y21*z2 + y1*z22 - z21*y2 - z1*y22;

    double  dmy1 = z11*x2 + z1*x12 - x11*z2 - x1*z12;
    double  dmy2 = z21*x2 + z1*x22 - x21*z2 - x1*z22;

    double  dmz1 = x11*y2 + x1*y12 - y11*x2 - y1*x12;
    double  dmz2 = x21*y2 + x1*y22 - y21*x2 - y1*x22;

    double  dJ1 = (mx*dmx1 + my*dmy1 + mz*dmz1)/J;
    double  dJ2 = (mx*dmx2 + my*dmy2 + mz*dmz2)/J;

    double  dnx1 = dmx1/J - mx*dJ1/J/J;
    double  dnx2 = dmx2/J - mx*dJ2/J/J;

    double  dny1 = dmy1/J - my*dJ1/J/J;
    double  dny2 = dmy2/J - my*dJ2/J/J;

    double  dnz1 = dmz1/J - mz*dJ1/J/J;
    double  dnz2 = dmz2/J - mz*dJ2/J/J;

    curvature = -0.5*(dnx1/x1 + dnx2/x2 + dny1/y1 + dny2/y2 + dnz1/z1 + dnz2/z2);

   return;
}


void BernsteinBasisFunsFaceQuad(int npElem, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double *tangent1, double *tangent2, double& curvature, double& J)
{
    vector<double>  dN1(npElem), dN2(npElem), d2N11, d2N22, d2N12;

    //BernsteinBasisFunsQuad(npElem, param[0], param[1], N, &dN1[0], &dN2[0], &d2N11[0], &d2N22[0], &d2N12[0]);
    BernsteinBasisFunsQuad(npElem, param[0], param[1], N, &dN1[0], &dN2[0]);

    double  x1=0.0, y1=0.0, z1=0.0, x2=0.0, y2=0.0, z2=0.0;
    double  x11=0.0, y11=0.0, z11=0.0, x22=0.0, y22=0.0, z22=0.0, x21=0.0, y21=0.0, z21=0.0;

    for(int ii=0; ii<npElem; ii++)
    {
      // dx_dxi1
      x1 += (xNode[ii] * dN1[ii]);
      y1 += (yNode[ii] * dN1[ii]);
      z1 += (zNode[ii] * dN1[ii]);

      // dx_dxi2
      x2 += (xNode[ii] * dN2[ii]);
      y2 += (yNode[ii] * dN2[ii]);
      z2 += (zNode[ii] * dN2[ii]);

      /*
      // d2x_dxi1dxi1
      x11 += (xNode[ii] * d2N11[ii]);
      y11 += (yNode[ii] * d2N11[ii]);
      z11 += (zNode[ii] * d2N11[ii]);

      // d2x_dxi2dxi2
      x22 += (xNode[ii] * d2N22[ii]);
      y22 += (yNode[ii] * d2N22[ii]);
      z22 += (zNode[ii] * d2N22[ii]);

      // d2x_dxi2dxi1
      x21 += (xNode[ii] * d2N12[ii]);
      y21 += (yNode[ii] * d2N12[ii]);
      z21 += (zNode[ii] * d2N12[ii]);
      */
    }

    double  x12 = x21, y12 = y21, z12 = z21;

    // Compute the normal

    double  mx = y1*z2 - z1*y2;
    double  my = z1*x2 - x1*z2;
    double  mz = x1*y2 - y1*x2;

    J = sqrt(mx*mx + my*my + mz*mz);

    normal[0] = mx/J;
    normal[1] = my/J;
    normal[2] = mz/J;

    // Compute the tangents
    tangent1[0] = 0.0;
    tangent1[1] = 0.0;
    tangent1[2] = 0.0;

    tangent2[0] = 0.0;
    tangent2[1] = 0.0;
    tangent2[2] = 0.0;

    // Compute the curvature

    double  dmx1 = y11*z2 + y1*z12 - z11*y2 - z1*y12;
    double  dmx2 = y21*z2 + y1*z22 - z21*y2 - z1*y22;

    double  dmy1 = z11*x2 + z1*x12 - x11*z2 - x1*z12;
    double  dmy2 = z21*x2 + z1*x22 - x21*z2 - x1*z22;

    double  dmz1 = x11*y2 + x1*y12 - y11*x2 - y1*x12;
    double  dmz2 = x21*y2 + x1*y22 - y21*x2 - y1*x22;

    double  dJ1 = (mx*dmx1 + my*dmy1 + mz*dmz1)/J;
    double  dJ2 = (mx*dmx2 + my*dmy2 + mz*dmz2)/J;

    double  dnx1 = dmx1/J - mx*dJ1/J/J;
    double  dnx2 = dmx2/J - mx*dJ2/J/J;

    double  dny1 = dmy1/J - my*dJ1/J/J;
    double  dny2 = dmy2/J - my*dJ2/J/J;

    double  dnz1 = dmz1/J - mz*dJ1/J/J;
    double  dnz2 = dmz2/J - mz*dJ2/J/J;

    curvature = -0.5*(dnx1/x1 + dnx2/x2 + dny1/y1 + dny2/y2 + dnz1/z1 + dnz2/z2);

    return;
}





#include "BasisFunctionsLagrange.h"
#include "headersEigen.h"
#include "ElementBase.h"

#include <math.h>
#include <cmath>
#include <stdio.h>
#include <vector>

using namespace std;




void Lagrange_BasisFuns1D(int npElem, double xi, double* N, double* dN_dxi, double* d2N_dxi2)
{
    double  fact1, fact2, val1, val2, val3, val4, xiSqr, xiCube;

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

            dN_dxi[0] = 0.0;

            d2N_dxi2[0] = 0.0;

        break;

        case 2:

            N[0] = 0.5*(1.0 - xi);
            N[1] = 0.5*(1.0 + xi);

            dN_dxi[0] = -0.5;
            dN_dxi[1] =  0.5;

            d2N_dxi2[0] = 0.0;
            d2N_dxi2[1] = 0.0;

        break;

        case 3:

            val1 = xi*xi;

            N[0] =  0.5*(val1 - xi);
            N[1] =  0.5*(val1 + xi);
            N[2] =  1.0 - val1;

            val1 = 2.0*xi;

            dN_dxi[0] =  0.5*(val1 - 1.0);
            dN_dxi[1] =  0.5*(val1 + 1.0);
            dN_dxi[2] = -val1;

            d2N_dxi2[0] =  1.0;
            d2N_dxi2[1] =  1.0;
            d2N_dxi2[2] = -2.0;

        break;

        case 4:

            fact1 = 1.0/16.0;
            xiCube = xi*xi*xi;
            xiSqr  = xi*xi;

            N[0] =  fact1 * ( -9*xiCube + 9*xiSqr +   xi - 1);
            N[1] =  fact1 * (  9*xiCube + 9*xiSqr -   xi - 1);
            N[2] =  fact1 * ( 27*xiCube - 9*xiSqr -27*xi + 9);
            N[3] =  fact1 * (-27*xiCube - 9*xiSqr +27*xi + 9);

            dN_dxi[0]   =  fact1 * (-27*xiSqr + 18*xi +  1);
            dN_dxi[1]   =  fact1 * ( 27*xiSqr + 18*xi -  1);
            dN_dxi[2]   =  fact1 * ( 81*xiSqr - 18*xi - 27);
            dN_dxi[3]   =  fact1 * (-81*xiSqr - 18*xi + 27);

            d2N_dxi2[0] =  fact1 * ( -54*xi + 18);
            d2N_dxi2[1] =  fact1 * (  54*xi + 18);
            d2N_dxi2[2] =  fact1 * ( 162*xi - 18);
            d2N_dxi2[3] =  fact1 * (-162*xi - 18);

        break;

        case 5:

            fact1 = 2.0/3.0;
            fact2 = 8.0/3.0;
            val1 = xi*xi;
            val2 = val1*xi;
            val3 = val2*xi;

            N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
            N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
            N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
            N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
            N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

            val4 = 4.0*val2;

            dN_dxi[0] =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
            dN_dxi[1] = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
            dN_dxi[2] =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
            dN_dxi[3] =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
            dN_dxi[4] = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);

            val4 = 12.0*val1;

            d2N_dxi2[0] =  fact1 * (-0.5  -  6.0*xi  +  val4);
            d2N_dxi2[1] = -fact2 * (-2.0  -  3.0*xi  +  val4);
            d2N_dxi2[2] =    4.0 * (-2.5  -  0.0     +  val4);
            d2N_dxi2[3] =  fact2 * ( 2.0  -  3.0*xi  -  val4);
            d2N_dxi2[4] = -fact1 * ( 0.5  -  6.0*xi  -  val4);

        break;

        default:

            printf("no basis functions defined for this npElem = %5d \n", npElem);

        break;
    }

    return;
}



void Lagrange_BasisFuns1D(int npElem, double xi, double* N, double* dN_dxi)
{
    double  fact1, fact2, val1, val2, val3, val4, xiSqr, xiCube;

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

            val1 = xi*xi;

            N[0] =  0.5*(val1 - xi);
            N[1] =  0.5*(val1 + xi);
            N[2] =  1.0 - val1;

            val1 = 2.0*xi;

            dN_dxi[0] =  0.5*(val1 - 1.0);
            dN_dxi[1] =  0.5*(val1 + 1.0);
            dN_dxi[2] = -val1;

        break;

        case 4:

            fact1 = 1.0/16.0;
            xiCube = xi*xi*xi;
            xiSqr  = xi*xi;

            N[0] =  fact1 * ( -9*xiCube + 9*xiSqr +   xi - 1);
            N[1] =  fact1 * (  9*xiCube + 9*xiSqr -   xi - 1);
            N[2] =  fact1 * ( 27*xiCube - 9*xiSqr -27*xi + 9);
            N[3] =  fact1 * (-27*xiCube - 9*xiSqr +27*xi + 9);

            dN_dxi[0]   =  fact1 * (-27*xiSqr + 18*xi +  1);
            dN_dxi[1]   =  fact1 * ( 27*xiSqr + 18*xi -  1);
            dN_dxi[2]   =  fact1 * ( 81*xiSqr - 18*xi - 27);
            dN_dxi[3]   =  fact1 * (-81*xiSqr - 18*xi + 27);

        break;

        case 5:

            fact1 = 2.0/3.0;
            fact2 = 8.0/3.0;
            val1 = xi*xi;
            val2 = val1*xi;
            val3 = val2*xi;

            N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
            N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
            N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
            N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
            N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

            val4 = 4.0*val2;

            dN_dxi[0] =  fact1 * (0.25 - 0.5*xi  - 3.0*val1  + val4);
            dN_dxi[1] = -fact2 * (0.5  - 2.0*xi  - 1.5*val1  + val4);
            dN_dxi[2] =    4.0 * (  0  - 2.5*xi  -   0.0     + val4);
            dN_dxi[3] =  fact2 * (0.5  + 2.0*xi  - 1.5*val1  - val4);
            dN_dxi[4] = -fact1 * (0.25 + 0.5*xi  - 3.0*val1  - val4);

        break;

        default:

            printf("no basis functions defined for this npElem = %5d \n", npElem);

        break;
    }

    return;
}



void Lagrange_BasisFuns1D(int npElem, double xi, double* N)
{
    double  fact1, fact2, val1, val2, val3, val4, xiSqr, xiCube;

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

            val1 = xi*xi;

            N[0] = 0.5 * (val1 - xi);
            N[1] = 0.5 * (val1 + xi);
            N[2] = 1.0 - val1;

        break;

        case 4:

            fact1 = 1.0/16.0;
            xiCube = xi*xi*xi;
            xiSqr  = xi*xi;

            N[0] =  fact1 * ( -9*xiCube + 9*xiSqr +   xi - 1);
            N[1] =  fact1 * (  9*xiCube + 9*xiSqr -   xi - 1);
            N[2] =  fact1 * ( 27*xiCube - 9*xiSqr -27*xi + 9);
            N[3] =  fact1 * (-27*xiCube - 9*xiSqr +27*xi + 9);

        break;

        case 5:

            fact1 = 2.0/3.0;
            fact2 = 8.0/3.0;
            val1 = xi*xi;
            val2 = val1*xi;
            val3 = val2*xi;

            N[0] =  fact1 * (0.25*xi  - 0.25*val1  -  val2     + val3);
            N[1] = -fact2 * (0.5 *xi  - val1       -  0.5*val2 + val3);
            N[2] =    4.0 * (0.25     - 1.25*val1  -  0        + val3);
            N[3] =  fact2 * (0.5*xi   + val1       -  0.5*val2 - val3);
            N[4] = -fact1 * (0.25*xi  + 0.25*val1  -  val2     - val3);

        break;

        default:

            printf("no basis functions defined for this npElem = %5d \n", npElem);

        break;
    }

    return;
}



void LagrangeBasisFunsTria(int npElem, double xi1, double xi2, double* N)
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

          N[0] = xi3*(2.0*xi3 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = 4.0*xi1*xi3;
          N[4] = 4.0*xi1*xi2;
          N[5] = 4.0*xi2*xi3;

      break;

      case 7:

          N[0] = xi3*(2.0*xi3 - 1.0) +  3.0*xi1*xi2*xi3;
          N[1] = xi1*(2.0*xi1 - 1.0) +  3.0*xi1*xi2*xi3;
          N[2] = xi2*(2.0*xi2 - 1.0) +  3.0*xi1*xi2*xi3;
          N[3] = 4.0*xi1*xi3         - 12.0*xi1*xi2*xi3;
          N[4] = 4.0*xi1*xi2         - 12.0*xi1*xi2*xi3;
          N[5] = 4.0*xi2*xi3         - 12.0*xi1*xi2*xi3;
          N[6] = 27.0*xi1*xi2*xi3;

      break;

      case 10:

          N[0] = 0.5*xi3*(2.0-9.0*xi3+9.0*xi3*xi3);
          N[1] = 0.5*xi1*(2.0-9.0*xi1+9.0*xi1*xi1);
          N[2] = 0.5*xi2*(2.0-9.0*xi2+9.0*xi2*xi2);
          N[3] = 4.5*xi1*xi3*(3*xi3-1.0);
          N[4] = 4.5*xi1*xi3*(3*xi1-1.0);
          N[5] = 4.5*xi1*xi2*(3*xi1-1.0);
          N[6] = 4.5*xi1*xi2*(3*xi2-1.0);
          N[7] = 4.5*xi2*xi3*(3*xi2-1.0);
          N[8] = 4.5*xi2*xi3*(3*xi3-1.0);
          N[9] = 27.0*xi1*xi2*xi3;

      break;

      default:

          printf("no basis functions defined for this npElem = %5d \n", npElem);

      break;
  }
  return;
}


void LagrangeBasisFunsTria(int npElem, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
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

          N[0] = xi3*(2.0*xi3 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = 4.0*xi1*xi3;
          N[4] = 4.0*xi1*xi2;
          N[5] = 4.0*xi2*xi3;

          dN_dxi1[0] = -4.0*xi3 + 1.0;
          dN_dxi1[1] =  4.0*xi1 - 1.0;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  4.0*(xi3 - xi1);
          dN_dxi1[4] =  4.0*xi2;
          dN_dxi1[5] = -4.0*xi2;

          dN_dxi2[0] = -4.0*xi3 + 1.0;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  4.0*xi2 - 1.0;
          dN_dxi2[3] = -4.0*xi1;
          dN_dxi2[4] =  4.0*xi1;
          dN_dxi2[5] =  4.0*(xi3 - xi2);

      break;

      case 7:

          N[0] = xi3*(2.0*xi3 - 1.0) +  3.0*xi1*xi2*xi3;
          N[1] = xi1*(2.0*xi1 - 1.0) +  3.0*xi1*xi2*xi3;
          N[2] = xi2*(2.0*xi2 - 1.0) +  3.0*xi1*xi2*xi3;
          N[3] = 4.0*xi1*xi3         - 12.0*xi1*xi2*xi3;
          N[4] = 4.0*xi1*xi2         - 12.0*xi1*xi2*xi3;
          N[5] = 4.0*xi2*xi3         - 12.0*xi1*xi2*xi3;
          N[6] = 27.0*xi1*xi2*xi3;

          dN_dxi1[0] = -4.0*xi3 + 1.0         +  3.0*(xi2*xi3-xi1*xi2);
          dN_dxi1[1] =  4.0*xi1 - 1.0         +  3.0*(xi2*xi3-xi1*xi2);
          dN_dxi1[2] =  0.0                   +  3.0*(xi2*xi3-xi1*xi2);
          dN_dxi1[3] =  4.0*(xi3 - xi1)       - 12.0*(xi2*xi3-xi1*xi2);
          dN_dxi1[4] =  4.0*xi2               - 12.0*(xi2*xi3-xi1*xi2);
          dN_dxi1[5] = -4.0*xi2               - 12.0*(xi2*xi3-xi1*xi2);
          dN_dxi1[6] =  27.0*(xi2*xi3-xi1*xi2);

          dN_dxi2[0] = -4.0*xi3 + 1.0         +  3.0*(xi1*xi3-xi1*xi2);
          dN_dxi2[1] =  0.0                   +  3.0*(xi1*xi3-xi1*xi2);
          dN_dxi2[2] =  4.0*xi2 - 1.0         +  3.0*(xi1*xi3-xi1*xi2);
          dN_dxi2[3] = -4.0*xi1               - 12.0*(xi1*xi3-xi1*xi2);
          dN_dxi2[4] =  4.0*xi1               - 12.0*(xi1*xi3-xi1*xi2);
          dN_dxi2[5] =  4.0*(xi3 - xi2)       - 12.0*(xi1*xi3-xi1*xi2);
          dN_dxi2[6] =  27.0*(xi1*xi3-xi1*xi2);

      break;

      case 10:

          N[0] = 0.5*xi3*(3*xi3-1)*(3*xi3-2);
          N[1] = 0.5*xi1*(3*xi1-1)*(3*xi1-2);
          N[2] = 0.5*xi2*(3*xi2-1)*(3*xi2-2);
          N[3] = 4.5*xi1*xi3*(3*xi3-1);
          N[4] = 4.5*xi1*xi3*(3*xi1-1);
          N[5] = 4.5*xi1*xi2*(3*xi1-1);
          N[6] = 4.5*xi1*xi2*(3*xi2-1);
          N[7] = 4.5*xi2*xi3*(3*xi2-1);
          N[8] = 4.5*xi2*xi3*(3*xi3-1);
          N[9] = 27.0*xi1*xi2*xi3;

          dN_dxi1[0] = -1.0+9.0*xi3-13.5*xi3*xi3;
          dN_dxi1[1] =  1.0-9.0*xi1+13.5*xi1*xi1;
          dN_dxi1[2] =  0.0;
          dN_dxi1[3] =  4.5*(xi1-xi3-6.0*xi1*xi3+3*xi3*xi3);
          dN_dxi1[4] =  4.5*(xi1-xi3+6.0*xi1*xi3-3*xi1*xi1);
          dN_dxi1[5] =  4.5*xi2*(6*xi1-1.0);
          dN_dxi1[6] =  4.5*xi2*(3*xi2-1.0);
          dN_dxi1[7] = -4.5*xi2*(3*xi2-1.0);
          dN_dxi1[8] = -4.5*xi2*(6*xi3-1.0);
          dN_dxi1[9] =  27.0*xi2*(xi3-xi1);

          dN_dxi2[0] = -1.0+9.0*xi3-13.5*xi3*xi3;
          dN_dxi2[1] =  0.0;
          dN_dxi2[2] =  1.0-9.0*xi2+13.5*xi2*xi2;
          dN_dxi2[3] = -4.5*xi1*(6*xi3-1.0);
          dN_dxi2[4] = -4.5*xi1*(3*xi1-1.0);
          dN_dxi2[5] =  4.5*xi1*(3*xi1-1.0);
          dN_dxi2[6] =  4.5*xi1*(6*xi2-1.0);
          dN_dxi2[7] =  4.5*(xi2-xi3+6.0*xi2*xi3-3*xi2*xi2);
          dN_dxi2[8] =  4.5*(xi2-xi3-6.0*xi2*xi3+3*xi3*xi3);
          dN_dxi2[9] =  27.0*xi1*(xi3-xi2);

      break;

      default:

          printf("no basis functions defined for this npElem = %5d \n", npElem);

      break;
  }
  return;
}





void LagrangeBasisFunsQuad(int npElem, double xi1, double xi2, double* N)
{
    int  ind = sqrt(npElem);
    double  Nu[ind], Nv[ind];

    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

        break;

        case 4:

            Lagrange_BasisFuns1D(2, xi1, Nu);
            Lagrange_BasisFuns1D(2, xi2, Nv);

            /*
             v2          u1,v2  u2,v2
             v1          u1,v1  u2,v1
               u1 u2
            */

            N[0] = Nv[0]*Nu[0];
            N[1] = Nv[0]*Nu[1];
            N[2] = Nv[1]*Nu[1];
            N[3] = Nv[1]*Nu[0];

        break;

        case 9:

            Lagrange_BasisFuns1D(3, xi1, Nu);
            Lagrange_BasisFuns1D(3, xi2, Nv);

            /*
             v2          u1,v2  u3,v2  u2,v2
             v3          u1,v3  u3,v3  u2,v3
             v1          u1,v1  u3,v1  u2,v1
               u1 u3 u2
            */

                  N[0] = Nv[0]*Nu[0];
                  N[4] = Nv[0]*Nu[2];
                  N[1] = Nv[0]*Nu[1];
                  //
                  N[7] = Nv[2]*Nu[0];
                  N[8] = Nv[2]*Nu[2];
                  N[5] = Nv[2]*Nu[1];
                  //
                  N[3] = Nv[1]*Nu[0];
                  N[6] = Nv[1]*Nu[2];
                  N[2] = Nv[1]*Nu[1];

        break;

        case 16:

            Lagrange_BasisFuns1D(4, xi1, Nu);
            Lagrange_BasisFuns1D(4, xi2, Nv);

            /*
             v2          u1,v2  u3,v2  u4,v2  u2,v2
             v4          u1,v4  u3,v4  u4,v4  u2,v4
             v3          u1,v3  u3,v3  u4,v3  u2,v3
             v1          u1,v1  u3,v1  u4,v1  u2,v1
               u1 u3 u4 u2
            */

                        N[ 0] = Nv[0]*Nu[0];
                        N[ 4] = Nv[0]*Nu[2];
                        N[ 5] = Nv[0]*Nu[3];
                        N[ 1] = Nv[0]*Nu[1];
                        //
                        N[11] = Nv[2]*Nu[0];
                        N[12] = Nv[2]*Nu[2];
                        N[13] = Nv[2]*Nu[3];
                        N[ 6] = Nv[2]*Nu[1];
                        //
                        N[10] = Nv[3]*Nu[0];
                        N[15] = Nv[3]*Nu[2];
                        N[14] = Nv[3]*Nu[3];
                        N[ 7] = Nv[3]*Nu[1];
                        //
                        N[ 3] = Nv[1]*Nu[0];
                        N[ 8] = Nv[1]*Nu[2];
                        N[ 9] = Nv[1]*Nu[3];
                        N[ 2] = Nv[1]*Nu[1];

        break;

        default:

            printf("no basis functions defined for this npElem = %5d \n", npElem);

        break;
    }

    return;
}



void LagrangeBasisFunsQuad(int npElem, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2)
{
/*
   3---------2     3----6----2     3---9---8--2
   |         |     |    |    |     |   |   |  |
   |         |     |    |    |    10--15--14--7
   |         |     7----8----5     |   |   |  |
   |         |     |    |    |    11--12--13--6
   |         |     |    |    |     |   |   |  |
   0---------1     0----4----1     0---4---5--1
*/

    int  ind = sqrt(npElem);
    double  Nu[ind], dNu[ind], Nv[ind], dNv[ind];


    switch(npElem)
    {
        case 1:

            N[0] = 1.0;

            dN_dxi1[0] = 0.0;
            dN_dxi2[0] = 0.0;

        break;

        case 4:

            Lagrange_BasisFuns1D(2, xi1, Nu, dNu);
            Lagrange_BasisFuns1D(2, xi2, Nv, dNv);

            /*
             v2          u1,v2  u2,v2
             v1          u1,v1  u2,v1
               u1 u2
            */

            N[0] = Nv[0]*Nu[0];
            N[1] = Nv[0]*Nu[1];
            N[3] = Nv[1]*Nu[0];
            N[2] = Nv[1]*Nu[1];

            dN_dxi1[0] = Nv[0]*dNu[0];
            dN_dxi1[1] = Nv[0]*dNu[1];
            dN_dxi1[3] = Nv[1]*dNu[0];
            dN_dxi1[2] = Nv[1]*dNu[1];

            dN_dxi2[0] = dNv[0]*Nu[0];
            dN_dxi2[1] = dNv[0]*Nu[1];
            dN_dxi2[3] = dNv[1]*Nu[0];
            dN_dxi2[2] = dNv[1]*Nu[1];

        break;

        case 9:

            Lagrange_BasisFuns1D(3, xi1, Nu, dNu);
            Lagrange_BasisFuns1D(3, xi2, Nv, dNv);

            /*
             v2          u1,v2  u3,v2  u2,v2
             v3          u1,v3  u3,v3  u2,v3
             v1          u1,v1  u3,v1  u2,v1
               u1 u3 u2
            */

                  N[0] = Nv[0]*Nu[0];
                  N[4] = Nv[0]*Nu[2];
                  N[1] = Nv[0]*Nu[1];
                  //
                  N[7] = Nv[2]*Nu[0];
                  N[8] = Nv[2]*Nu[2];
                  N[5] = Nv[2]*Nu[1];
                  //
                  N[3] = Nv[1]*Nu[0];
                  N[6] = Nv[1]*Nu[2];
                  N[2] = Nv[1]*Nu[1];

            dN_dxi1[0] = Nv[0]*dNu[0];
            dN_dxi1[4] = Nv[0]*dNu[2];
            dN_dxi1[1] = Nv[0]*dNu[1];
            //
            dN_dxi1[7] = Nv[2]*dNu[0];
            dN_dxi1[8] = Nv[2]*dNu[2];
            dN_dxi1[5] = Nv[2]*dNu[1];
            //
            dN_dxi1[3] = Nv[1]*dNu[0];
            dN_dxi1[6] = Nv[1]*dNu[2];
            dN_dxi1[2] = Nv[1]*dNu[1];

            dN_dxi2[0] = dNv[0]*Nu[0];
            dN_dxi2[4] = dNv[0]*Nu[2];
            dN_dxi2[1] = dNv[0]*Nu[1];
            //
            dN_dxi2[7] = dNv[2]*Nu[0];
            dN_dxi2[8] = dNv[2]*Nu[2];
            dN_dxi2[5] = dNv[2]*Nu[1];
            //
            dN_dxi2[3] = dNv[1]*Nu[0];
            dN_dxi2[6] = dNv[1]*Nu[2];
            dN_dxi2[2] = dNv[1]*Nu[1];

        break;

        case 16:

            Lagrange_BasisFuns1D(4, xi1, Nu, dNu);
            Lagrange_BasisFuns1D(4, xi2, Nv, dNv);

            /*
             v2          u1,v2  u3,v2  u4,v2  u2,v2
             v4          u1,v4  u3,v4  u4,v4  u2,v4
             v3          u1,v3  u3,v3  u4,v3  u2,v3
             v1          u1,v1  u3,v1  u4,v1  u2,v1
               u1 u3 u4 u2
            */

                        N[ 0] = Nv[0]*Nu[0];
                        N[ 4] = Nv[0]*Nu[2];
                        N[ 5] = Nv[0]*Nu[3];
                        N[ 1] = Nv[0]*Nu[1];
                        //
                        N[11] = Nv[2]*Nu[0];
                        N[12] = Nv[2]*Nu[2];
                        N[13] = Nv[2]*Nu[3];
                        N[ 6] = Nv[2]*Nu[1];
                        //
                        N[10] = Nv[3]*Nu[0];
                        N[15] = Nv[3]*Nu[2];
                        N[14] = Nv[3]*Nu[3];
                        N[ 7] = Nv[3]*Nu[1];
                        //
                        N[ 3] = Nv[1]*Nu[0];
                        N[ 8] = Nv[1]*Nu[2];
                        N[ 9] = Nv[1]*Nu[3];
                        N[ 2] = Nv[1]*Nu[1];

                  dN_dxi1[ 0] = Nv[0]*dNu[0];
                  dN_dxi1[ 4] = Nv[0]*dNu[2];
                  dN_dxi1[ 5] = Nv[0]*dNu[3];
                  dN_dxi1[ 1] = Nv[0]*dNu[1];
                  //
                  dN_dxi1[11] = Nv[2]*dNu[0];
                  dN_dxi1[12] = Nv[2]*dNu[2];
                  dN_dxi1[13] = Nv[2]*dNu[3];
                  dN_dxi1[ 6] = Nv[2]*dNu[1];
                  //
                  dN_dxi1[10] = Nv[3]*dNu[0];
                  dN_dxi1[15] = Nv[3]*dNu[2];
                  dN_dxi1[14] = Nv[3]*dNu[3];
                  dN_dxi1[ 7] = Nv[3]*dNu[1];
                  //
                  dN_dxi1[ 3] = Nv[1]*dNu[0];
                  dN_dxi1[ 8] = Nv[1]*dNu[2];
                  dN_dxi1[ 9] = Nv[1]*dNu[3];
                  dN_dxi1[ 2] = Nv[1]*dNu[1];

                  dN_dxi2[ 0] = dNv[0]*Nu[0];
                  dN_dxi2[ 4] = dNv[0]*Nu[2];
                  dN_dxi2[ 5] = dNv[0]*Nu[3];
                  dN_dxi2[ 1] = dNv[0]*Nu[1];
                  //
                  dN_dxi2[11] = dNv[2]*Nu[0];
                  dN_dxi2[12] = dNv[2]*Nu[2];
                  dN_dxi2[13] = dNv[2]*Nu[3];
                  dN_dxi2[ 6] = dNv[2]*Nu[1];
                  //
                  dN_dxi2[10] = dNv[3]*Nu[0];
                  dN_dxi2[15] = dNv[3]*Nu[2];
                  dN_dxi2[14] = dNv[3]*Nu[3];
                  dN_dxi2[ 7] = dNv[3]*Nu[1];
                  //
                  dN_dxi2[ 3] = dNv[1]*Nu[0];
                  dN_dxi2[ 8] = dNv[1]*Nu[2];
                  dN_dxi2[ 9] = dNv[1]*Nu[3];
                  dN_dxi2[ 2] = dNv[1]*Nu[1];

        break;


        default:

            printf("no basis functions defined for this npElem = %5d \n", npElem);

        break;
    }

    return;
}



void LagrangeBasisFunsTetra(int npElem, double xi1, double xi2, double xi3, double* N)
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

          N[0] = xi4*(2.0*xi4 - 1.0);
          N[1] = xi1*(2.0*xi1 - 1.0);
          N[2] = xi2*(2.0*xi2 - 1.0);
          N[3] = xi3*(2.0*xi3 - 1.0);
          N[4] = 4.0*xi1*xi4;
          N[5] = 4.0*xi1*xi2;
          N[6] = 4.0*xi2*xi4;
          N[7] = 4.0*xi4*xi3;
          N[8] = 4.0*xi2*xi3;
          N[9] = 4.0*xi1*xi3;

      break;


      case 20:

          N[ 0] = 0.5*xi4*(3*xi4-1)*(3*xi4-2);
          N[ 1] = 0.5*xi2*(3*xi2-1)*(3*xi2-2);
          N[ 2] = 0.5*xi3*(3*xi3-1)*(3*xi3-2);
          N[ 3] = 0.5*xi1*(3*xi1-1)*(3*xi1-2);
          N[ 4] = 4.5*xi2*xi4*(3*xi4-1);
          N[ 5] = 4.5*xi2*xi4*(3*xi2-1);
          N[ 6] = 4.5*xi2*xi3*(3*xi2-1);
          N[ 7] = 4.5*xi2*xi3*(3*xi3-1);
          N[ 8] = 4.5*xi3*xi4*(3*xi3-1);
          N[ 9] = 4.5*xi3*xi4*(3*xi4-1);
          N[10] = 4.5*xi1*xi4*(3*xi1-1);
          N[11] = 4.5*xi1*xi4*(3*xi4-1);
          N[12] = 4.5*xi1*xi3*(3*xi1-1);
          N[13] = 4.5*xi1*xi3*(3*xi3-1);
          N[14] = 4.5*xi1*xi2*(3*xi1-1);
          N[15] = 4.5*xi1*xi2*(3*xi2-1);
          N[16] = 27*xi2*xi3*xi4;
          N[17] = 27*xi1*xi2*xi4;
          N[18] = 27*xi1*xi3*xi4;
          N[19] = 27*xi1*xi2*xi3;

      break;

      default:

          printf("no basis functions defined for this npElem = %5d \n", npElem);

      break;
  }

}





void LagrangeBasisFunsTetra(int npElem, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
  double  xi4 = 1.0 - xi1 - xi2 - xi3, fact1, fact2;

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

                N[0] = xi4*(2*xi4-1);
                N[1] = xi1*(2*xi1-1);
                N[2] = xi2*(2*xi2-1);
                N[3] = xi3*(2*xi3-1);
                N[4] = 4*xi1*xi4;
                N[5] = 4*xi1*xi2;
                N[6] = 4*xi2*xi4;
                N[7] = 4*xi3*xi4;
                N[8] = 4*xi2*xi3;
                N[9] = 4*xi1*xi3;

          dN_dxi1[0] = -4*xi4 + 1;
          dN_dxi1[1] =  4*xi1 - 1;
          dN_dxi1[2] =  0;
          dN_dxi1[3] =  0;
          dN_dxi1[4] =  4*(xi4-xi1);
          dN_dxi1[5] =  4*xi2;
          dN_dxi1[6] = -4*xi2;
          dN_dxi1[7] = -4*xi3;
          dN_dxi1[8] =  0;
          dN_dxi1[9] =  4*xi3;

          dN_dxi2[0] = -4*xi4 + 1;
          dN_dxi2[1] =  0;
          dN_dxi2[2] =  4*xi2 - 1;
          dN_dxi2[3] =  0;
          dN_dxi2[4] = -4*xi1;
          dN_dxi2[5] =  4*xi1;
          dN_dxi2[6] =  4*(xi4-xi2);
          dN_dxi2[7] = -4*xi3;
          dN_dxi2[8] =  4*xi3;
          dN_dxi2[9] =  0;

          dN_dxi3[0] = -4*xi4 + 1;
          dN_dxi3[1] =  0;
          dN_dxi3[2] =  0;
          dN_dxi3[3] =  4*xi3 - 1;
          dN_dxi3[4] = -4*xi1;
          dN_dxi3[5] =  0;
          dN_dxi3[6] = -4*xi2;
          dN_dxi3[7] =  4*(xi4-xi3);
          dN_dxi3[8] =  4*xi2;
          dN_dxi3[9] =  4*xi1;

      break;

      case 11:

               fact1 = 32*xi1*xi2*xi3*xi4;
               fact2 = 64*xi1*xi2*xi3*xi4;

               N[ 0] = xi4*(2*xi4 - 1) + fact1;
               N[ 1] = xi1*(2*xi1 - 1) + fact1;
               N[ 2] = xi2*(2*xi2 - 1) + fact1;
               N[ 3] = xi3*(2*xi3 - 1) + fact1;
               N[ 4] = 4*xi1*xi4       - fact2;
               N[ 5] = 4*xi1*xi2       - fact2;
               N[ 6] = 4*xi2*xi4       - fact2;
               N[ 7] = 4*xi3*xi4       - fact2;
               N[ 8] = 4*xi2*xi3       - fact2;
               N[ 9] = 4*xi1*xi3       - fact2;
               N[10] = 4*fact2;

               fact1 = 32*xi2*xi3*(xi4-xi1);
               fact2 = 64*xi2*xi3*(xi4-xi1);

          dN_dxi1[ 0] = -4*xi4 + 1    + fact1;
          dN_dxi1[ 1] =  4*xi1 - 1    + fact1;
          dN_dxi1[ 2] =  0            + fact1;
          dN_dxi1[ 3] =  0            + fact1;
          dN_dxi1[ 4] =  4*(xi4-xi1)  - fact2;
          dN_dxi1[ 5] =  4*xi2        - fact2;
          dN_dxi1[ 6] = -4*xi2        - fact2;
          dN_dxi1[ 7] = -4*xi3        - fact2;
          dN_dxi1[ 8] =  0            - fact2;
          dN_dxi1[ 9] =  4*xi3        - fact2;
          dN_dxi1[10] =  4*fact2;

               fact1 = 32*xi1*xi3*(xi4-xi2);
               fact2 = 64*xi1*xi3*(xi4-xi2);

          dN_dxi2[ 0] = -4*xi4 + 1    + fact1;
          dN_dxi2[ 1] =  0            + fact1;
          dN_dxi2[ 2] =  4*xi2 - 1    + fact1;
          dN_dxi2[ 3] =  0            + fact1;
          dN_dxi2[ 4] = -4*xi1        - fact2;
          dN_dxi2[ 5] =  4*xi1        - fact2;
          dN_dxi2[ 6] =  4*(xi4-xi2)  - fact2;
          dN_dxi2[ 7] = -4*xi3        - fact2;
          dN_dxi2[ 8] =  4*xi3        - fact2;
          dN_dxi2[ 9] =  0            - fact2;
          dN_dxi2[10] =  4*fact2;

               fact1 = 32*xi1*xi2*(xi4-xi3);
               fact2 = 64*xi1*xi2*(xi4-xi3);

          dN_dxi3[ 0] = -4*xi4 + 1    + fact1;
          dN_dxi3[ 1] =  0            + fact1;
          dN_dxi3[ 2] =  0            + fact1;
          dN_dxi3[ 3] =  4*xi3 - 1    + fact1;
          dN_dxi3[ 4] = -4*xi1        - fact2;
          dN_dxi3[ 5] =  0            - fact2;
          dN_dxi3[ 6] = -4*xi2        - fact2;
          dN_dxi3[ 7] =  4*(xi4-xi3)  - fact2;
          dN_dxi3[ 8] =  4*xi2        - fact2;
          dN_dxi3[ 9] =  4*xi1        - fact2;
          dN_dxi3[10] =  4*fact2;

      break;

      case 20:

                N[ 0] = 0.5*xi4*(3*xi4-1)*(3*xi4-2);
                N[ 1] = 0.5*xi2*(3*xi2-1)*(3*xi2-2);
                N[ 2] = 0.5*xi3*(3*xi3-1)*(3*xi3-2);
                N[ 3] = 0.5*xi1*(3*xi1-1)*(3*xi1-2);
                N[ 4] = 4.5*xi2*xi4*(3*xi4-1);
                N[ 5] = 4.5*xi2*xi4*(3*xi2-1);
                N[ 6] = 4.5*xi2*xi3*(3*xi2-1);
                N[ 7] = 4.5*xi2*xi3*(3*xi3-1);
                N[ 8] = 4.5*xi3*xi4*(3*xi3-1);
                N[ 9] = 4.5*xi3*xi4*(3*xi4-1);
                N[10] = 4.5*xi1*xi4*(3*xi1-1);
                N[11] = 4.5*xi1*xi4*(3*xi4-1);
                N[12] = 4.5*xi1*xi3*(3*xi1-1);
                N[13] = 4.5*xi1*xi3*(3*xi3-1);
                N[14] = 4.5*xi1*xi2*(3*xi1-1);
                N[15] = 4.5*xi1*xi2*(3*xi2-1);
                N[16] = 27*xi2*xi3*xi4;
                N[17] = 27*xi1*xi2*xi4;
                N[18] = 27*xi1*xi3*xi4;
                N[19] = 27*xi1*xi2*xi3;


          dN_dxi1[ 0] = -1.0+9*xi4-13.5*xi4*xi4;
          dN_dxi1[ 1] =  0.0;
          dN_dxi1[ 2] =  0.0;
          dN_dxi1[ 3] =  1.0-9*xi1+13.5*xi1*xi1;
          dN_dxi1[ 4] = -4.5*xi2*(6*xi4-1);
          dN_dxi1[ 5] = -4.5*xi2*(3*xi2-1);
          dN_dxi1[ 6] =  0.0;
          dN_dxi1[ 7] =  0.0;
          dN_dxi1[ 8] = -4.5*xi3*(3*xi3-1);
          dN_dxi1[ 9] = -4.5*xi3*(6*xi4-1);
          dN_dxi1[10] =  4.5*xi1-4.5*xi4+27*xi1*xi4-13.5*xi1*xi1;
          dN_dxi1[11] =  4.5*xi1-4.5*xi4-27*xi1*xi4+13.5*xi4*xi4;
          dN_dxi1[12] =  4.5*xi3*(6*xi1-1);
          dN_dxi1[13] =  4.5*xi3*(3*xi3-1);
          dN_dxi1[14] =  4.5*xi2*(6*xi1-1);
          dN_dxi1[15] =  4.5*xi2*(3*xi2-1);
          dN_dxi1[16] = -27*xi2*xi3;
          dN_dxi1[17] =  27*xi2*(xi4-xi1);
          dN_dxi1[18] =  27*xi3*(xi4-xi1);
          dN_dxi1[19] =  27*xi2*xi3;


          dN_dxi2[ 0] = -1.0+9*xi4-13.5*xi4*xi4;
          dN_dxi2[ 1] =  1.0-9*xi2+13.5*xi2*xi2;
          dN_dxi2[ 2] =  0.0;
          dN_dxi2[ 3] =  0.0;
          dN_dxi2[ 4] =  4.5*xi2-4.5*xi4-27*xi2*xi4+13.5*xi4*xi4;
          dN_dxi2[ 5] =  4.5*xi2-4.5*xi4+27*xi2*xi4-13.5*xi2*xi2;
          dN_dxi2[ 6] =  4.5*xi3*(6*xi2-1);
          dN_dxi2[ 7] =  4.5*xi3*(3*xi3-1);
          dN_dxi2[ 8] = -4.5*xi3*(3*xi3-1);
          dN_dxi2[ 9] = -4.5*xi3*(6*xi4-1);
          dN_dxi2[10] = -4.5*xi1*(3*xi1-1);
          dN_dxi2[11] = -4.5*xi1*(6*xi4-1);
          dN_dxi2[12] =  0;
          dN_dxi2[13] =  0;
          dN_dxi2[14] =  4.5*xi1*(3*xi1-1);
          dN_dxi2[15] =  4.5*xi1*(6*xi2-1);
          dN_dxi2[16] =  27*xi3*(xi4-xi2);
          dN_dxi2[17] =  27*xi1*(xi4-xi2);
          dN_dxi2[18] = -27*xi1*xi3;
          dN_dxi2[19] =  27*xi1*xi3;


          dN_dxi3[ 0] = -1.0+9*xi4-13.5*xi4*xi4;
          dN_dxi3[ 1] =  0.0;
          dN_dxi3[ 2] =  1.0-9*xi3+13.5*xi3*xi3;
          dN_dxi3[ 3] =  0.0;
          dN_dxi3[ 4] = -4.5*xi2*(6*xi4-1);
          dN_dxi3[ 5] = -4.5*xi2*(3*xi2-1);
          dN_dxi3[ 6] =  4.5*xi2*(3*xi2-1);
          dN_dxi3[ 7] =  4.5*xi2*(6*xi3-1);
          dN_dxi3[ 8] =  4.5*xi3-4.5*xi4+27*xi3*xi4-13.5*xi3*xi3;
          dN_dxi3[ 9] =  4.5*xi3-4.5*xi4-27*xi3*xi4+13.5*xi4*xi4;
          dN_dxi3[10] = -4.5*xi1*(3*xi1-1);
          dN_dxi3[11] = -4.5*xi1*(6*xi4-1);
          dN_dxi3[12] =  4.5*xi1*(3*xi1-1);
          dN_dxi3[13] =  4.5*xi1*(6*xi3-1);
          dN_dxi3[14] =  0.0;
          dN_dxi3[15] =  0.0;
          dN_dxi3[16] =  27*xi2*(xi4-xi3);
          dN_dxi3[17] = -27*xi1*xi2 ;
          dN_dxi3[18] =  27*xi1*(xi4-xi3);
          dN_dxi3[19] =  27*xi1*xi2;

      break;

      default:

          printf("no basis functions defined for this npElem = %5d \n", npElem);

      break;
  }

}



void LagrangeBasisFunsHexa(int npElem, double xi1, double xi2, double xi3, double* N)
{
    int ind = pow(npElem,1/3);
    double  Nu[ind], Nv[ind], Nw[ind];

    Lagrange_BasisFuns1D(ind, xi1, Nu);
    Lagrange_BasisFuns1D(ind, xi2, Nv);
    Lagrange_BasisFuns1D(ind, xi3, Nw);

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

            printf("no basis functions defined for this npElem = %5d in 'LagrangeBasisFunsHexa' \n", npElem);

        break;
    }

    return;
}




void LagrangeBasisFunsHexa(int npElem, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
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

    int ind;
    if(npElem == 8) ind=2;
    else if(npElem == 27) ind=3;

    double  Nu[ind], dNu[ind], Nv[ind], dNv[ind], Nw[ind], dNw[ind];

    Lagrange_BasisFuns1D(ind, xi1, Nu, dNu);
    Lagrange_BasisFuns1D(ind, xi2, Nv, dNv);
    Lagrange_BasisFuns1D(ind, xi3, Nw, dNw);

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

            printf("no basis functions defined for this npElem = %5d in 'LagrangeBasisFunsHexa' \n", npElem);

        break;
    }

    return;
}






void LagrangeBasisFunsPrism(int npElem, double xi1, double xi2, double xi4, double* N)
{
  // xi1, xi2 and xi3 are the parametric coordinates for the triangle, and
  // xi4 is the parametric coordinate in the direction normal to the triangle
  double  xi3 = 1.0- xi1 - xi2;

  switch(npElem)
  {
      case 6:

          N[0] = xi3*(0.5*(1.0-xi4));
          N[1] = xi1*(0.5*(1.0-xi4));
          N[2] = xi2*(0.5*(1.0-xi4));
          N[3] = xi3*(0.5*(1.0+xi4));
          N[4] = xi1*(0.5*(1.0+xi4));
          N[5] = xi2*(0.5*(1.0+xi4));

      break;

      default:

          printf("no basis functions defined for this npElem = %5d in 'LagrangeBasisFunsPrism' \n", npElem);

      break;
  }
  return;
}




void LagrangeBasisFunsPrism(int npElem, double xi1, double xi2, double xi4, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi4)
{
  // xi1, xi2 and xi3 are the parametric coordinates for the triangle, and
  // xi4 is the parametric coordinate in the direction normal to the triangle
  double  xi3 = 1.0- xi1 - xi2;

  switch(npElem)
  {
      case 6:

          N[0] = xi3*(0.5*(1.0-xi4));
          N[1] = xi1*(0.5*(1.0-xi4));
          N[2] = xi2*(0.5*(1.0-xi4));
          N[3] = xi3*(0.5*(1.0+xi4));
          N[4] = xi1*(0.5*(1.0+xi4));
          N[5] = xi2*(0.5*(1.0+xi4));


          dN_dxi1[0]   = (-1.0)*(0.5*(1.0-xi4));
          dN_dxi1[1]   = ( 1.0)*(0.5*(1.0-xi4));
          dN_dxi1[2]   = ( 0.0)*(0.5*(1.0-xi4));
          dN_dxi1[3]   = (-1.0)*(0.5*(1.0+xi4));
          dN_dxi1[4]   = ( 1.0)*(0.5*(1.0+xi4));
          dN_dxi1[5]   = ( 0.0)*(0.5*(1.0+xi4));

          dN_dxi2[0]  = (-1.0)*(0.5*(1.0-xi4));
          dN_dxi2[1]  = ( 0.0)*(0.5*(1.0-xi4));
          dN_dxi2[2]  = ( 1.0)*(0.5*(1.0-xi4));
          dN_dxi2[3]  = (-1.0)*(0.5*(1.0+xi4));
          dN_dxi2[4]  = ( 0.0)*(0.5*(1.0+xi4));
          dN_dxi2[5]  = ( 1.0)*(0.5*(1.0+xi4));

          dN_dxi4[0]  = xi3*(-0.5);
          dN_dxi4[1]  = xi1*(-0.5);
          dN_dxi4[2]  = xi2*(-0.5);
          dN_dxi4[3]  = xi3*( 0.5);
          dN_dxi4[4]  = xi1*( 0.5);
          dN_dxi4[5]  = xi2*( 0.5);

      break;

      default:

          printf("no basis functions defined for this npElem = %5d in 'LagrangeBasisFunsPrism' \n", npElem);

      break;
  }
  return;
}




void LagrangeBasisFunsPyramid(int npElem, double xi1, double xi2, double xi3, double* N)
{
 double  v11, v12, v21, v22, v31, v32;

  switch(npElem)
  {
      case 5:

          v11 = 1.0 - xi1;
          v12 = 1.0 + xi1;
          v21 = 1.0 - xi2;
          v22 = 1.0 + xi2;
          v31 = 1.0 - xi3;
          v32 = 1.0 + xi3;

          N[0] = 0.125*v11*v21*v31;
          N[1] = 0.125*v12*v21*v31;
          N[2] = 0.125*v12*v22*v31;
          N[3] = 0.125*v11*v22*v31;
          N[4] = 0.5*v32;

      break;

      default:

          printf("no basis functions defined for this npElem = %5d in 'LagrangeBasisFunsPyramid' \n", npElem);

      break;
  }

  return;
}






void LagrangeBasisFunsPyramid(int npElem, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3)
{
 double  v11, v12, v21, v22, v31, v32;

  switch(npElem)
  {
      case 5:

          v11 = 1.0 - xi1;
          v12 = 1.0 + xi1;
          v21 = 1.0 - xi2;
          v22 = 1.0 + xi2;
          v31 = 1.0 - xi3;
          v32 = 1.0 + xi3;

          N[0] = 0.125*v11*v21*v31;
          N[1] = 0.125*v12*v21*v31;
          N[2] = 0.125*v12*v22*v31;
          N[3] = 0.125*v11*v22*v31;
          N[4] = 0.5*v32;

          dN_dxi1[0] = -0.125*v21*v31;
          dN_dxi1[1] =  0.125*v21*v31;
          dN_dxi1[2] =  0.125*v22*v31;
          dN_dxi1[3] = -0.125*v22*v31;
          dN_dxi1[4] =  0.0;

          dN_dxi2[0] = -0.125*v11*v31;
          dN_dxi2[1] = -0.125*v12*v31;
          dN_dxi2[2] =  0.125*v11*v31;
          dN_dxi2[3] =  0.125*v12*v31;
          dN_dxi2[4] =  0.0;

          dN_dxi3[0] = -0.125*v11*v21;
          dN_dxi3[1] = -0.125*v12*v21;
          dN_dxi3[2] = -0.125*v11*v22;
          dN_dxi3[3] = -0.125*v12*v22;
          dN_dxi3[4] =  0.5;

      break;

      default:

          printf("no basis functions defined for this npElem = %5d in 'LagrangeBasisFunsPyramid' \n", npElem);

      break;
  }

  return;
}




void LagrangeBasisFunsLine1D(int npElem, double uu, double *xx, double *N, double *dN_dx, double& Jac)
{
  double dx_du, dN1[npElem];

  Lagrange_BasisFuns1D(npElem, uu, N, dN1);

  if(npElem == 1)
  {
    Jac = 1.0;
  }
  else
  {
    Jac = 0.0;
    for(int ii=0; ii<npElem; ii++)
      Jac +=  (xx[ii] * dN1[ii]);
  }

  dx_du = 1.0/Jac ;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(int ii=0; ii<npElem; ii++)
    dN_dx[ii] = dN1[ii] * dx_du ;

  return;
}




void LagrangeBasisFunsLine3D(int npElem, double uu, double *xx, double* yy, double* zz, double *N, double *dN_dx, double& Jac)
{
   double  du_dx, dx, dy, dz, dN1[npElem];

   Lagrange_BasisFuns1D(npElem, uu, N, dN1);

   if(npElem == 1)
   {
     Jac = 1.0;
   }
   else
   {
     dx = dy = dz = 0.0;
     for(int ii=0; ii<npElem; ii++)
     {
       dx +=  (xx[ii] * dN1[ii]);
       dy +=  (yy[ii] * dN1[ii]);
       dz +=  (zz[ii] * dN1[ii]);
     }
     Jac = sqrt(dx*dx+dy*dy+dz*dz);
   }

   du_dx = 1.0/Jac ;

  // Compute derivatives of basis functions w.r.t physical coordinates
  for(int ii=0; ii<npElem; ii++)
    dN_dx[ii] = dN1[ii] * du_dx ;

  return;
}





void LagrangeBasisFunsEdge2D(int npElem, double* param, double *xNode, double* yNode, double *N, double *normal, double& Jac)
{
   if(npElem == 1)
   {
       cerr << "LagrangeBasisFunsEdge2D is not available for npElem = " << npElem << endl;
   }

   double  dx, dy, dN1[npElem];

   Lagrange_BasisFuns1D(npElem, param[0], N, dN1);

   dx = dy = 0.0;
   for(int ii=0; ii<npElem; ii++)
   {
     dx +=  (xNode[ii] * dN1[ii]);
     dy +=  (yNode[ii] * dN1[ii]);
   }
   Jac = sqrt(dx*dx+dy*dy);

   // Compute the normal
   normal[0] =  dy/Jac;
   normal[1] = -dx/Jac;

   return;
}




void LagrangeBasisFunsEdge3D(int p, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac)
{
  return;
}



void LagrangeBasisFunsFace(int npElem, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double* tangent1, double* tangent2, double& Jac)
{
   double  dN1[npElem], dN2[npElem];

   if( (npElem == 3) || (npElem == 6) || (npElem == 10) )
     LagrangeBasisFunsTria(npElem, param[0], param[1], N, dN1, dN2);
   else if( (npElem == 4) || (npElem == 9) )
     LagrangeBasisFunsQuad(npElem, param[0], param[1], N, dN1, dN2);

   Vector3d  du, dv;
   du.setZero();
   dv.setZero();
   for(int ii=0; ii<npElem; ii++)
   {
     du(0) += (xNode[ii] * dN1[ii]);
     du(1) += (yNode[ii] * dN1[ii]);
     du(2) += (zNode[ii] * dN1[ii]);

     dv(0) += (xNode[ii] * dN2[ii]);
     dv(1) += (yNode[ii] * dN2[ii]);
     dv(2) += (zNode[ii] * dN2[ii]);
   }

   Vector3d normal1 = dv.cross(du);

   Jac = normal1.norm();
   normal1 /= Jac;

   // Compute tangent1
   tangent1[0] = du[0];
   tangent1[1] = du[1];
   tangent1[2] = du[2];

   // Compute tangent2
   tangent2[0] = dv[0];
   tangent2[1] = dv[1];
   tangent2[2] = dv[2];

   // Compute the normal
   normal[0] = normal1[0];
   normal[1] = normal1[1];
   normal[2] = normal1[2];

   return;
}








void shape_functions_Derivatives_Shell(int type, int npElem, double* param, double *xNode, double *yNode, double *N, double *dN_dx, double *dN_dy, double& Jac)
{
    double  dN_du1[npElem], dN_du2[npElem];
    double  xx, yy, detinv, B[2][2], Binv[2][2] ;

    if(type == ELEM_SHAPE_TRIA) // triangular elements
    {
      LagrangeBasisFunsTria(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else  if(type == ELEM_SHAPE_QUAD) // Lagrange quad elements
    {
      LagrangeBasisFunsQuad(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }

    // Gradient of mapping from parameter space to physical space
    B[0][0] = B[1][0] = B[0][1] = B[1][1] = 0.0;

    for(int ii=0; ii<npElem; ii++)
    {
        xx = xNode[ii];
        yy = yNode[ii];

        B[0][0] +=  (xx * dN_du1[ii]) ;
        B[1][0] +=  (xx * dN_du2[ii]) ;
        B[0][1] +=  (yy * dN_du1[ii]) ;
        B[1][1] +=  (yy * dN_du2[ii]) ;
    }

    Jac  = B[0][0]*B[1][1] - B[0][1]*B[1][0];

    detinv = 1.0/Jac ;

    Binv[0][0] =  B[1][1] * detinv;
    Binv[0][1] = -B[0][1] * detinv;
    Binv[1][0] = -B[1][0] * detinv;
    Binv[1][1] =  B[0][0] * detinv;

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(int ii=0; ii<npElem; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv[0][0] + dN_du2[ii] * Binv[0][1];
      dN_dy[ii] = dN_du1[ii] * Binv[1][0] + dN_du2[ii] * Binv[1][1];
    }

    return;
}



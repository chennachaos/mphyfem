#ifndef incl_utilities_cfd_h
#define incl_utilities_cfd_h


#include "headersBasic.h"
#include "headersEigen.h"


int  pointInsideTria3Node(double* xNode, double* yNode, double* target);

int  pointInverseTria3node(double* xNode, double* yNode, double* target, double* param);

int  pointInverseTria6node(double* xNode, double* yNode, double* target, double* param);

int  pointInverseTetra10node(double* xNode, double* yNode, double* zNode, double* target, double* param);

int  pointInverseQuad4node(double* xNode, double* yNode, double* target, double* param);


inline void printMatrix(MatrixXd& AA)
{
    int ii, jj;
    printf("\n\n");
    for(ii=0;ii<AA.rows();ii++)
    {
        for(jj=0;jj<AA.cols();jj++)
           printf("\t%14.8f", AA(ii,jj));
        printf("\n");
    }
    printf("\n\n");

    return;
}


inline void printVector(VectorXd& AA)
{
    printf("\n\n");
    for(int ii=0;ii<AA.rows();ii++)
        printf("\t%6d\t%12.8f\n", ii, AA(ii));
    printf("\n\n");

   return;
}



inline void printVector(vector<int>&  vec)
{
    printf("\n\n");
    for(int ii=0;ii<vec.size();ii++)
        printf("\t%6d ", vec[ii]);
    printf("\n\n");

   return;
}


inline void printVector(double* data, int nn)
{
    printf("\n\n");
    for(int ii=0;ii<nn;ii++)
      printf("\t%6d\t%12.8f\n", ii, data[ii]);
    printf("\n\n");

   return;
}



inline  void aitken_accelerator(const VectorXd&  x1, const VectorXd&  x2, VectorXd&  x)
{
  // Applies Aitken acceleration to the sequence in x,
  // returning a sequence of n-2 new approximations,
  // where n is the length of the vector x.

  for(int k=0; k<x.rows(); k++)
  {
    double val = (x(k) - 2*x1(k) + x2(k));
    if(abs(val) > 1.0e-10)
      x(k) = x(k) - (x(k) - x1(k))*(x(k) - x1(k))/val;;
  }

  return;
}




#endif

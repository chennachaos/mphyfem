
#ifndef incl_BasisFunctionsBernstein_h
#define incl_BasisFunctionsBernstein_h


/**
  Computes function values of
  univariate Bernstein polynomials of a given degree
*/
void Bernstein_BasisFuns1D(int degree, double xi, double* N);

/**
  Computes function values and first derivative of
  univariate Bernstein polynomials of a given degree
*/
void Bernstein_BasisFuns1D(int degree, double xi, double* N, double* dN_dxi);

/**
  Computes function values and first & second derivative of
  univariate Bernstein polynomials of a given degree
*/
void Bernstein_BasisFuns1D(int degree, double xi, double* N, double* dN_dxi, double* d2N_dxi2);

/**
  Computes function values of
  Bernstein polynomials of a given degree for Triangular elements
*/
void BernsteinBasisFunsTria(int degree, double xi1, double xi2, double* N);

/**
  Computes function values and first derivative of
  Bernstein polynomials of a given degree for Triangular elements
*/
void BernsteinBasisFunsTria(int degree, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2);

/**
  Computes function values of
  Bernstein polynomials of a triangular face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsTria(int degree, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2, double* d2N_dxi1dxi1, double* d2N_dxi2dxi2, double* d2N_dxi1dxi2);

/**
  Computes function values of
  Bernstein polynomials of a given degree for Quadrilateral elements
*/
void BernsteinBasisFunsQuad(int degree, double xi1, double xi2, double* N);

/**
  Computes function values and first derivative of
  Bernstein polynomials of a given degree for Quadrilateral elements
*/
void BernsteinBasisFunsQuad(int degree, double xi1, double xi2, double* N, double* dN_dxi1, double* dN_dxi2);

/**
  Computes function values of
  Bernstein polynomials of a given degree for Tetrahedral elements
*/
void BernsteinBasisFunsTetra(int degree, double xi1, double xi2, double xi3, double* N);

/**
  Computes function values and first derivative of
  Bernstein polynomials of a given degree for Tetrahedral elements
*/
void BernsteinBasisFunsTetra(int degree, double xi1, double xi2, double xi3, double* N, double* dN_dxi, double* dN_dxi2, double* dN_dxi3);

/**
  Computes function values of
  Bernstein polynomials of a given degree for Hexahedral elements
*/
void BernsteinBasisFunsHexa(int degree, double xi1, double xi2, double xi3, double* N);

/**
  Computes function values and first derivative of
  Bernstein polynomials of a given degree for Hexahedral elements
*/
void BernsteinBasisFunsHexa(int degree, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3);

/**
  Computes function values of
  Bernstein polynomials of a given degree for Wedge elements
*/
void BernsteinBasisFunsWedge(int degree, double xi1, double xi2, double xi3, double* N);

/**
  Computes function values and first derivative of
  Bernstein polynomials of a given degree for Wedge elements
*/
void BernsteinBasisFunsWedge(int degree, double xi1, double xi2, double xi3, double* N, double* dN_dxi1, double* dN_dxi2, double* dN_dxi3);

/**
  Computes function values of
  Bernstein polynomials of an edge in 2D.
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsEdge2D(int degree, double* param, double *xNode, double* yNode, double *N, double *dNdx, double *dNdy, double *normal, double& curvature, double& Jac);

/**
  Computes function values of
  Bernstein polynomials of an edge in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsEdge3D(int degree, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);

/**
  Computes function values of
  Bernstein polynomials of a triangular face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsFaceTria(int degree, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double *tangent1, double *tangent2, double& curvature, double& Jac);

/**
  Computes function values of
  Bernstein polynomials of a quadrilateral face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void BernsteinBasisFunsFaceQuad(int degree, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double *tangent1, double *tangent2, double& curvature, double& Jac);



#endif


#ifndef incl_BasisFunctionsLagrange_h
#define incl_BasisFunctionsLagrange_h

/**
  Computes function values and the first derivative (wrt to spatial coordinates)
  for basis functions for the 2D problem
*/
int computeBasisFunctions2D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double*  N, double*  dN_dx, double* dN_dy, double&  Jac);


/**
  Computes function values and the first derivative (wrt to spatial coordinates)
  for basis functions for the 3D problem
*/
int computeBasisFunctions3D(bool flag, int ETYPE, int degree, double* param, double* xNode, double* yNode, double* zNode, double*  N, double*  dN_dx, double* dN_dy, double* dN_dz, double&  Jac);


/**
  Computes function values, first derivative and second derivative of
  univariate Lagrange polynomials based on given number of nodes per element
*/
void Lagrange_BasisFuns1D(int npElem, double xi, double* N, double* dN_dxi, double* d2N_dxi2);

/**
  Computes function values and first derivative of
  univariate Lagrange polynomials based on given number of nodes per element
*/
void Lagrange_BasisFuns1D(int npElem, double xi, double* N, double* dN_dxi);

/**
  Computes function values of
  univariate Lagrange polynomials based on given number of nodes per element
*/
void Lagrange_BasisFuns1D(int npElem, double xi, double* N);


/**
  Computes function values of
  Lagrange Triangular elements based on given number of nodes per element
*/
void LagrangeBasisFunsTria(int npElem, double xi, double zeta, double* N);

/**
  Computes function values and first derivative of
  Lagrange Triangular elements based on given number of nodes per element
*/
void LagrangeBasisFunsTria(int npElem, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);


/**
  Computes function values of
  Lagrange Quadrilateral elements based on given number of nodes per element
*/
void LagrangeBasisFunsQuad(int npElem, double xi, double zeta, double* N);

/**
  Computes function values and first derivative of
  Lagrange Quadrilateral elements based on given number of nodes per element
*/
void LagrangeBasisFunsQuad(int npElem, double xi, double zeta, double* N, double* dN_dxi, double* dN_dzeta);


/**
  Computes function values of
  Lagrange Tetrahedral elements based on given number of nodes per element
*/
void LagrangeBasisFunsTetra(int npElem, double xi1, double xi2, double xi3, double* N);

/**
  Computes function values and first derivative of
  Lagrange Tetrahedral elements based on given number of nodes per element
*/
void LagrangeBasisFunsTetra(int npElem, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values of
  Lagrange Hexahedral elements based on given number of nodes per element
*/
void LagrangeBasisFunsHexa(int npElem, double xi, double zeta, double eta, double* N);

/**
  Computes function values and first derivative of
  Lagrange Hexahedral elements based on given number of nodes per element
*/
void LagrangeBasisFunsHexa(int npElem, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values of
  Lagrange Penta/wedge elements based on given number of nodes per element
*/
void LagrangeBasisFunsPrism(int npElem, double xi, double zeta, double eta, double* N);

/**
  Computes function values and first derivative of
  Lagrange Penta/wedge elements based on given number of nodes per element
*/
void LagrangeBasisFunsPrism(int npElem, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values and first derivative of
  Lagrange Pyramid elements based on given number of nodes per element
*/
void LagrangeBasisFunsPyramid(int npElem, double xi, double zeta, double eta, double* N, double* dN_dxi, double* dN_dzeta, double* dN_deta);

/**
  Computes function values of
  Lagrange Pyramid elements based on given number of nodes per element
*/
void LagrangeBasisFunsPyramid(int npElem, double xi, double zeta, double eta, double* N);


void LagrangeBasisFunsLine1D(int npElem, double uu, double *xNode, double *N, double *dN_dx, double& Jac);


/**
  Computes function values of
  Lagrange polynomials of an edge in 2D.
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void LagrangeBasisFunsEdge2D(int npElem, double* param, double *xNode, double* yNode, double *N, double *normal, double& Jac);

/**
  Computes function values of
  Lagrange polynomials of an edge in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void LagrangeBasisFunsEdge3D(int npElem, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double& Jac);

/**
  Computes function values of
  Lagrange element face in 3D
  This function is for computing the required basis functions for calculating
  the contribution from traction forces and/or applying Dirichlet BCs using Penalty/Nitsche method
*/
void LagrangeBasisFunsFace(int npElem, double* param, double *xNode, double* yNode, double *zNode, double *N, double *normal, double* tangent1, double* tangent2, double& Jac);


/**
 Computes the values and derivatives of basis functions
 for shell elements in local coordinate system
*/
void shape_functions_Derivatives_Shell(int type, int npElem, double* param, double *xNode, double *yNode, double *N, double *dN_dx, double *dN_dy, double& Jac);



#endif


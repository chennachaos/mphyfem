#include "BasisFunctionsLagrange.h"
#include "BasisFunctionsBernstein.h"
#include "headersBasic.h"
#include "headersEigen.h"
#include "ElementBase.h"




int computeBasisFunctions2D(bool flag, int ELEM_TYPE, int npElem, double* param, double* xNode, double* yNode, double*  N, double*  dN_dx, double* dN_dy, double&  Jac)
{
    double  dN_du1[npElem], dN_du2[npElem];
    double  xx, yy;
    MatrixXd B(2,2), Binv(2,2) ;

    //cout << "ELEM_TYPE = " << ELEM_TYPE << '\t' << npElem << '\t' << param[0] << '\t' << param[1] << endl;

    if( ELEM_TYPE == ELEM_SHAPE_TRIA )
    {
      LagrangeBasisFunsTria(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    if( ELEM_TYPE == ELEM_SHAPE_TRIA_BERNSTEIN )
    {
      BernsteinBasisFunsTria(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else if( ELEM_TYPE == ELEM_SHAPE_QUAD )
    {
      LagrangeBasisFunsQuad(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else if( ELEM_TYPE == ELEM_SHAPE_QUAD_BERNSTEIN )
    {
      BernsteinBasisFunsQuad(npElem, param[0], param[1], N, dN_du1, dN_du2);
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions2D... " << endl;
      exit(-2);
    }

    // Gradient of mapping from parameter space to physical space

    B.setZero();
    for(int ii=0; ii<npElem; ii++)
    {
      xx = xNode[ii];
      yy = yNode[ii];

      //cout << dN_du1[ii] << '\t' << dN_du2[ii] << endl;

      B(0,0) +=  (xx * dN_du1[ii]) ;
      B(1,0) +=  (xx * dN_du2[ii]) ;
      B(0,1) +=  (yy * dN_du1[ii]) ;
      B(1,1) +=  (yy * dN_du2[ii]) ;
    }

    Jac  = B.determinant();
    Binv = B.inverse();

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(int ii=0; ii<npElem; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv(0,0) + dN_du2[ii] * Binv(0,1);
      dN_dy[ii] = dN_du1[ii] * Binv(1,0) + dN_du2[ii] * Binv(1,1);
    }

    return 0;
}



int computeBasisFunctions3D(bool flag, int ELEM_TYPE, int npElem, double* param, double* xNode, double* yNode, double* zNode, double*  N, double*  dN_dx, double* dN_dy, double* dN_dz, double&  Jac)
{
    double  dN_du1[npElem], dN_du2[npElem], dN_du3[npElem];
    double  xx, yy, zz;
    MatrixXd B(3,3), Binv(3,3) ;

    if( ELEM_TYPE == ELEM_SHAPE_TETRA)
    {
      LagrangeBasisFunsTetra(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    if( ELEM_TYPE == ELEM_SHAPE_TETRA_BERNSTEIN)
    {
      BernsteinBasisFunsTetra(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(ELEM_TYPE == ELEM_SHAPE_WEDGE)
    {
      LagrangeBasisFunsPrism(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(ELEM_TYPE == ELEM_SHAPE_HEXA)
    {
      LagrangeBasisFunsHexa(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else if(ELEM_TYPE == ELEM_SHAPE_HEXA_BERNSTEIN)
    {
      BernsteinBasisFunsHexa(npElem, param[0], param[1], param[2], N, dN_du1, dN_du2, dN_du3);
    }
    else
    {
      cerr << " Invalid Element Type in computeBasisFunctions3D... " << endl;
      exit(-2);
    }

    // Gradient of mapping from parameter space to physical space

    B.setZero();
    for(int ii=0; ii<npElem; ii++)
    {
      xx = xNode[ii];
      yy = yNode[ii];
      zz = zNode[ii];

      B(0,0) +=  (xx * dN_du1[ii]) ;
      B(1,0) +=  (xx * dN_du2[ii]) ;
      B(2,0) +=  (xx * dN_du3[ii]) ;

      B(0,1) +=  (yy * dN_du1[ii]) ;
      B(1,1) +=  (yy * dN_du2[ii]) ;
      B(2,1) +=  (yy * dN_du3[ii]) ;

      B(0,2) +=  (zz * dN_du1[ii]) ;
      B(1,2) +=  (zz * dN_du2[ii]) ;
      B(2,2) +=  (zz * dN_du3[ii]) ;
    }

    Jac  = B.determinant();
    Binv = B.inverse();

    // Compute derivatives of basis functions w.r.t physical coordinates
    for(int ii=0; ii<npElem; ii++)
    {
      dN_dx[ii] = dN_du1[ii] * Binv(0,0) + dN_du2[ii] * Binv(0,1) + dN_du3[ii] * Binv(0,2);
      dN_dy[ii] = dN_du1[ii] * Binv(1,0) + dN_du2[ii] * Binv(1,1) + dN_du3[ii] * Binv(1,2);
      dN_dz[ii] = dN_du1[ii] * Binv(2,0) + dN_du2[ii] * Binv(2,1) + dN_du3[ii] * Binv(2,2);
    }

    return 0;
}










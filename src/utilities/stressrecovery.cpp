
#include "stressrecovery.h"
#include <vector>
#include <math.h>
#include <iostream>
#include "BasisFunctionsBernstein.h"


using namespace std;


void gausspoints_extrapolate_Triangle(int ngp, vector<double>& gpts1, vector<double>& gpts2)
{
    switch(ngp)
    {
        case 3:

            gpts1[0] = 0.0;    gpts2[0] = 0.0;
            gpts1[1] = 1.0;    gpts2[1] = 0.0;
            gpts1[2] = 0.0;    gpts2[2] = 1.0;

        break;

        case 6:

            gpts1[0] = -1.0/6.0;    gpts2[0] = -1.0/6.0;
            gpts1[1] =  11.0/6.0;   gpts2[1] = -1.0/6.0;
            gpts1[2] = -1.0/6.0;    gpts2[2] =  11.0/6.0;

            gpts1[3] =  5.0/6.0;    gpts2[3] = -1.0/6.0;
            gpts1[4] =  5.0/6.0;    gpts2[4] =  5.0/6.0;
            gpts1[5] = -1.0/6.0;    gpts2[5] =  5.0/6.0;

        break;

        default:
            cerr << " 'gausspoints_extrapolate_Triangle' not implemented for ngp = " << ngp << endl;
            exit(1);
        break;
    }
    return;
}



void stressrecovery_extrapolate_Triangle(int degree, double* outval, double* vals2project)
{
    int ii, jj, ngp;
    vector<double>  gpts1, gpts2;
    vector<vector<double> >  N;

    switch(degree)
    {
        case 1:
            // value at the centroid of the triangle is projected to 
            // the 3 nodes of the triangle
            for(ii=0;ii<3;ii++)
              vals2project[ii] = outval[ii];
        break;

      case 2:
            // values at the 3 Gauss points are projected to 
            // the 6 nodes of the triangle

            //get the coordinates for nodes in the Gauss quadrature element space
            ngp = 6;
            gpts1.resize(ngp);
            gpts2.resize(ngp);
            gausspoints_extrapolate_Triangle(ngp, gpts1, gpts2);

            //compute the basis function values at the coordinates obtained above

            N.resize(6);
            for(ii=0; ii<6; ii++)
            {
              N[ii].resize(3);

              // Here, 3 point quadrature is used. So, the degree of basis functions
              // is one less than that of the actual element
              BernsteinBasisFunsTria(degree-1, gpts1[ii], gpts2[ii], &(N[ii][0]));
            }

            //project values from Gauss points to nodes
            for(ii=0; ii<6; ii++)
            {
              vals2project[ii] = 0.0;
              for(jj=0; jj<3; jj++)
                vals2project[ii] += N[ii][jj]*outval[jj];
            }

        break;

        default:
            cerr << " 'stressrecovery_extrapolate_Triangle' not implemented for degree = " << degree << endl;
            exit(1);
        break;
    }
    return;
}




void gausspoints_extrapolate_Quadrilateral(int ngp, vector<double>& gpts1, vector<double>& gpts2)
{
    double  pp, mm, zz;
    switch(ngp)
    {
        case 4:

            pp = sqrt(3.0); mm = -pp;

            gpts1[0]  = mm;  gpts2[0] = mm;
            gpts1[1]  = pp;  gpts2[1] = mm;
            gpts1[2]  = pp;  gpts2[2] = pp;
            gpts1[3]  = mm;  gpts2[3] = pp;

        break;

        case 9:

            pp = sqrt(5.0/3.0); zz = 0.0; mm = -pp;

            gpts1[0]  = mm;  gpts2[0] = mm;
            gpts1[1]  = pp;  gpts2[1] = mm;
            gpts1[2]  = pp;  gpts2[2] = pp;
            gpts1[3]  = mm;  gpts2[3] = pp;

            gpts1[4]  = zz;  gpts2[4] = mm;
            gpts1[5]  = pp;  gpts2[5] = zz;
            gpts1[6]  = zz;  gpts2[6] = pp;
            gpts1[7]  = mm;  gpts2[7] = zz;

            gpts1[8]  = zz;  gpts2[8] = zz;

        break;

        default:
            cerr << " 'gausspoints_extrapolate_Quadrilateral' not implemented for ngp = " << ngp << endl;
            exit(1);
        break;
    }
    return;
}




void stressrecovery_extrapolate_Quadrilateral(int degree, double* outval, double* vals2project)
{
    int ii, jj, ngp;
    if(degree ==  1) ngp = 4;
    else if(degree ==  2) ngp = 9;
    else
    {
        cerr << " 'stressrecovery_extrapolate_Quadrilateral' not implemented for degree = " << degree << endl;
        exit(1);
    }

    vector<double>  gpts1(ngp), gpts2(ngp), N(ngp);

    //get the coordinates for nodes in the Gauss quadrature element space
    gausspoints_extrapolate_Quadrilateral(ngp, gpts1, gpts2);

    for(ii=0; ii<ngp; ii++)
    {
        //compute the basis function values at the coordinates obtained above
        BernsteinBasisFunsQuad(degree, gpts1[ii], gpts2[ii], &(N[0]));

        //project values from Gauss points to nodes
        vals2project[ii] = 0.0;
        for(jj=0; jj<ngp; jj++)
        {
            vals2project[ii] += N[jj]*outval[jj];
        }
    }

    return;
}



void gausspoints_extrapolate_Hexahedron(int ngp, vector<double>& gpts1, vector<double>& gpts2, vector<double>& gpts3)
{
    double  pp, mm, zz;
    switch(ngp)
    {
        case 8:

            pp = sqrt(3.0); mm = -pp;

            gpts1[0]  = mm;    gpts2[0] = mm;    gpts3[0] = mm;
            gpts1[1]  = pp;    gpts2[1] = mm;    gpts3[1] = mm;
            gpts1[3]  = mm;    gpts2[3] = pp;    gpts3[3] = mm;
            gpts1[2]  = pp;    gpts2[2] = pp;    gpts3[2] = mm;

            gpts1[4]  = mm;    gpts2[4] = mm;    gpts3[4] = pp;
            gpts1[5]  = pp;    gpts2[5] = mm;    gpts3[5] = pp;
            gpts1[7]  = mm;    gpts2[7] = pp;    gpts3[7] = pp;
            gpts1[6]  = pp;    gpts2[6] = pp;    gpts3[6] = pp;

        break;

        case 27:

            pp = sqrt(5.0/3.0); zz = 0.0; mm = -pp;

            gpts1[0]  = mm;    gpts2[0]  = mm;    gpts3[0]  = mm;
            gpts1[8]  = zz;    gpts2[8]  = mm;    gpts3[8]  = mm;
            gpts1[1]  = pp;    gpts2[1]  = mm;    gpts3[1]  = mm;
            gpts1[9]  = mm;    gpts2[9]  = zz;    gpts3[9]  = mm;
            gpts1[20] = zz;    gpts2[20] = zz;    gpts3[20] = mm;
            gpts1[11] = pp;    gpts2[11] = zz;    gpts3[11] = mm;
            gpts1[3]  = mm;    gpts2[3]  = pp;    gpts3[3]  = mm;
            gpts1[13] = zz;    gpts2[13] = pp;    gpts3[13] = mm;
            gpts1[2]  = pp;    gpts2[2]  = pp;    gpts3[2]  = mm;

            gpts1[10] = mm;    gpts2[10] = mm;    gpts3[10] = zz;
            gpts1[21] = zz;    gpts2[21] = mm;    gpts3[21] = zz;
            gpts1[12] = pp;    gpts2[12] = mm;    gpts3[12] = zz;
            gpts1[22] = mm;    gpts2[22] = zz;    gpts3[22] = zz;
            gpts1[26] = zz;    gpts2[26] = zz;    gpts3[26] = zz;
            gpts1[23] = pp;    gpts2[23] = zz;    gpts3[23] = zz;
            gpts1[15] = mm;    gpts2[15] = pp;    gpts3[15] = zz;
            gpts1[24] = zz;    gpts2[24] = pp;    gpts3[24] = zz;
            gpts1[14] = pp;    gpts2[14] = pp;    gpts3[14] = zz;

            gpts1[4]  = mm;    gpts2[4]  = mm;    gpts3[4]  = pp;
            gpts1[16] = zz;    gpts2[16] = mm;    gpts3[16] = pp;
            gpts1[5]  = pp;    gpts2[5]  = mm;    gpts3[5]  = pp;
            gpts1[17] = mm;    gpts2[17] = zz;    gpts3[17] = pp;
            gpts1[25] = zz;    gpts2[25] = zz;    gpts3[25] = pp;
            gpts1[18] = pp;    gpts2[18] = zz;    gpts3[18] = pp;
            gpts1[7]  = mm;    gpts2[7]  = pp;    gpts3[7]  = pp;
            gpts1[19] = zz;    gpts2[19] = pp;    gpts3[19] = pp;
            gpts1[6]  = pp;    gpts2[6]  = pp;    gpts3[6]  = pp;

        break;

        default:
            cerr << " 'gausspoints_extrapolate_Hexahedron' not implemented for ngp = " << ngp << endl;
            exit(1);

        break;
    }
    return;
}





void stressrecovery_extrapolate_Hexahedron(int degree, double* outval, double* vals2project)
{
    int ii, jj, ngp;
    if(degree ==  1) ngp = 8;
    else if(degree ==  2) ngp = 27;
    else
    {
        cerr << " 'stressrecovery_extrapolate_Hexahedron' not implemented for degree = " << degree << endl;
        exit(1);
    }

    vector<double>  gpts1(ngp), gpts2(ngp), gpts3(ngp), N(ngp);

    //get the coordinates for nodes in the Gauss quadrature element space
    gausspoints_extrapolate_Hexahedron(ngp, gpts1, gpts2, gpts3);

    for(ii=0; ii<ngp; ii++)
    {
        //compute the basis function values at the coordinates obtained above
        BernsteinBasisFunsHexa(degree, gpts1[ii], gpts2[ii], gpts3[ii], &(N[0]));

        //project values from Gauss points to nodes
        vals2project[ii] = 0.0;
        for(jj=0; jj<ngp; jj++)
        {
            //cout << ii << '\t' << jj << '\t' << N[jj] << '\t' << outval[jj] << endl;
            vals2project[ii] += N[jj]*outval[jj];
        }
    }

    return;
}




void gausspoints_extrapolate_Tetrahedron(int ngp, vector<double>& gpts1, vector<double>& gpts2, vector<double>& gpts3)
{
    switch(ngp)
    {
        case 4:

            gpts1[0]  = 0.0; gpts2[0] = 0.0; gpts3[0] = 0.0;
            gpts1[1]  = 1.0; gpts2[1] = 0.0; gpts3[1] = 0.0;
            gpts1[2]  = 0.0; gpts2[2] = 1.0; gpts3[2] = 0.0;
            gpts1[3]  = 0.0; gpts2[3] = 0.0; gpts3[3] = 1.0;

        break;

        case 10:

            gpts1[0]  = -1.0/6.0;   gpts2[0] = -1.0/6.0;   gpts3[0] = -1.0/6.0;
            gpts1[1]  = 11.0/6.0;   gpts2[1] = -1.0/6.0;   gpts3[1] = -1.0/6.0;
            gpts1[2]  = -1.0/6.0;   gpts2[2] = 11.0/6.0;   gpts3[2] = -1.0/6.0;
            gpts1[3]  = -1.0/6.0;   gpts2[3] = -1.0/6.0;   gpts3[3] = 11.0/6.0;
            gpts1[4]  =  5.0/6.0;   gpts2[4] = -1.0/6.0;   gpts3[4] = -1.0/6.0;
            gpts1[5]  =  5.0/6.0;   gpts2[5] =  5.0/6.0;   gpts3[5] = -1.0/6.0;
            gpts1[6]  = -1.0/6.0;   gpts2[6] =  5.0/6.0;   gpts3[6] = -1.0/6.0;
            gpts1[7]  = -1.0/6.0;   gpts2[7] = -1.0/6.0;   gpts3[7] =  5.0/6.0;
            gpts1[8]  =  5.0/6.0;   gpts2[8] = -1.0/6.0;   gpts3[8] =  5.0/6.0;
            gpts1[9]  = -1.0/6.0;   gpts2[9] =  5.0/6.0;   gpts3[9] =  5.0/6.0;

        break;

        default:
            cerr << " 'gausspoints_extrapolate_Quadrilateral' not implemented for ngp = " << ngp << endl;
            exit(1);
        break;
    }
    return;
}



void stressrecovery_extrapolate_Tetrahedron(int degree, double* outval, double* vals2project)
{
    int ii, jj, ngp;
    vector<double>  gpts1, gpts2, gpts3;
    vector<double>  N;

    switch(degree)
    {
        case 1:
            // value at the centroid of the triangle is projected to 
            // the 4 nodes of the triangle
            for(ii=0;ii<4;ii++)
              vals2project[ii] = outval[ii];
        break;

        case 2:
            // values at the 4 Gauss points are projected to 
            // the 10 nodes of the triangle

            //get the coordinates for nodes in the Gauss quadrature element space
            ngp = 10;
            gpts1.resize(ngp);
            gpts2.resize(ngp);
            gpts3.resize(ngp);
            gausspoints_extrapolate_Tetrahedron(ngp, gpts1, gpts2, gpts3);

            //compute the basis function values at the coordinates obtained above

            N.resize(4);
            for(ii=0; ii<10; ii++)
            {
                // Here, 4 point quadrature is used. So, the degree of basis functions
                // is one less than that of the actual element
                BernsteinBasisFunsTetra(degree-1, gpts1[ii], gpts2[ii], gpts3[ii], &(N[0]));

                //project values from Gauss points to nodes
                vals2project[ii] = 0.0;
                for(jj=0; jj<4; jj++)
                  vals2project[ii] += N[jj]*outval[jj];
            }

        break;

        default:
            cerr << " 'stressrecovery_extrapolate_Tetrahedron' not implemented for degree = " << degree << endl;
            exit(1);
        break;
    }
    return;
}




void gausspoints_extrapolate_Wedge(int ngp, vector<double>& gpts1, vector<double>& gpts2, vector<double>& gpts3)
{
    double  pp = sqrt(5.0/3.0), zz = 0.0, mm = -pp;

    switch(ngp)
    {
        case 2:

            gpts1[0]  = 0.0;    gpts2[0] = 0.0;    gpts3[0] = -sqrt(3);
            gpts1[1]  = 0.0;    gpts2[1] = 0.0;    gpts3[1] =  sqrt(3);

        break;

        case 18:

            pp = sqrt(5.0/3.0); zz = 0.0; mm = -pp;

            gpts1[0]  = -1.0/6.0;    gpts2[0] = -1.0/6.0;    gpts3[0] = mm;
            gpts1[1]  =  11.0/6.0;   gpts2[1] = -1.0/6.0;    gpts3[1] = mm;
            gpts1[2]  = -1.0/6.0;    gpts2[2] =  11.0/6.0;   gpts3[2] = mm;
            gpts1[6]  =  5.0/6.0;    gpts2[6] = -1.0/6.0;    gpts3[6] = mm;
            gpts1[9]  =  5.0/6.0;    gpts2[9] =  5.0/6.0;    gpts3[9] = mm;
            gpts1[7]  = -1.0/6.0;    gpts2[7] =  5.0/6.0;    gpts3[7] = mm;

            gpts1[8]  = -1.0/6.0;    gpts2[8]  = -1.0/6.0;    gpts3[8]  = zz;
            gpts1[10] =  11.0/6.0;   gpts2[10] = -1.0/6.0;    gpts3[10] = zz;
            gpts1[11] = -1.0/6.0;    gpts2[11] =  11.0/6.0;   gpts3[11] = zz;
            gpts1[15] =  5.0/6.0;    gpts2[15] = -1.0/6.0;    gpts3[15] = zz;
            gpts1[17] =  5.0/6.0;    gpts2[17] =  5.0/6.0;    gpts3[17] = zz;
            gpts1[16] = -1.0/6.0;    gpts2[16] =  5.0/6.0;    gpts3[16] = zz;

            gpts1[3]  = -1.0/6.0;    gpts2[3]  = -1.0/6.0;    gpts3[3]  = pp;
            gpts1[4]  =  11.0/6.0;   gpts2[4]  = -1.0/6.0;    gpts3[4]  = pp;
            gpts1[5]  = -1.0/6.0;    gpts2[5]  =  11.0/6.0;   gpts3[5]  = pp;
            gpts1[12] =  5.0/6.0;    gpts2[12] = -1.0/6.0;    gpts3[12] = pp;
            gpts1[14] =  5.0/6.0;    gpts2[14] =  5.0/6.0;    gpts3[14] = pp;
            gpts1[13] = -1.0/6.0;    gpts2[13] =  5.0/6.0;    gpts3[13] = pp;

        break;

        default:
            cerr << " 'gausspoints_extrapolate_Wedge' not implemented for ngp = " << ngp << endl;
            exit(1);
        break;
    }
    return;
}



void stressrecovery_extrapolate_Wedge(int degree, double* outval, double* vals2project)
{
    int ii, jj, ngp;
    vector<double>  gpts1, gpts2, gpts3;
    vector<double>  N;

    switch(degree)
    {
        case 1:
            // value at the centroid of the triangle is projected to 
            // the 3 nodes of the triangle
            vals2project[0] = outval[0];
            vals2project[1] = outval[0];
            vals2project[2] = outval[0];
            vals2project[3] = outval[1];
            vals2project[4] = outval[1];
            vals2project[5] = outval[1];

        break;

      case 2:
            // values at the 3 Gauss points are projected to 
            // the 6 nodes of the triangle

            //get the coordinates for nodes in the Gauss quadrature element space
            ngp = 18;
            gpts1.resize(ngp);
            gpts2.resize(ngp);
            gpts3.resize(ngp);
            gausspoints_extrapolate_Wedge(ngp, gpts1, gpts2, gpts3);

            //compute the basis function values at the coordinates obtained above

            N.resize(6);
            for(ii=0; ii<18; ii++)
            {
                // Here, 3 point quadrature is used. So, the degree of basis functions
                // is one less than that of the actual element
                BernsteinBasisFunsWedge(degree-1, gpts1[ii], gpts2[ii], gpts3[ii], &(N[0]));

                //project values from Gauss points to nodes
                vals2project[ii] = 0.0;
                for(jj=0; jj<6; jj++)
                  vals2project[ii] += N[jj]*outval[jj];
            }

        break;

        default:
            cerr << " 'stressrecovery_extrapolate_Wedge' not implemented for degree = " << degree << endl;
            exit(1);
        break;
    }
    return;
}




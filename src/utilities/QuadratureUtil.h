/*=============================================================================
        File: QuadratureUtil.h
  Created by: Chennakesava Kadapa          (12 Sept 2015)
 Purpose    : Header file for the definitions of
              Quadrature utility functions

 ============================================================================*/
#ifndef incl_QuadratureUtil_H
#define incl_QuadratureUtil_H

#include <vector>

using std::vector;


void getGaussPoints1D(int ngp, vector<double>& gausspoints, vector<double>& gaussweights);


int  getGaussPoints2D(int npElem, vector<double>& gps1, vector<double>& gps2, vector<double>& gws);


void getGaussPointsTriangle(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gws);


void getGaussPointsQuad(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gws);


int  getGaussPoints3D(int npElem, vector<double>& gps1, vector<double>& gps2, vector<double>& gps3, vector<double>& gws);


void getGaussPointsTetra(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gps3, vector<double>& gws);


void getGaussPointsPyramid(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gps3, vector<double>& gws);


void getGaussPointsWedge(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gps3, vector<double>& gws);


void getGaussPointsHexa(int ngp, vector<double>& gps1, vector<double>& gps2, vector<double>& gps3, vector<double>& gws);


#endif  //incl_QuadratureUtil_H






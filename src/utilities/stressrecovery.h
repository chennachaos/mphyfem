
#ifndef incl_stressrecovery_h
#define incl_stressrecovery_h


#include <vector>
#include <math.h>

using namespace std;


void gausspoints_extrapolate_Triangle(int degree, vector<double>& gpts1, vector<double>& gpts2);

void stressrecovery_extrapolate_Triangle(int degree, double* outval, double* vals2project);

void gausspoints_extrapolate_Quadrilateral(int degree, vector<double>& gpts1, vector<double>& gpts2);

void stressrecovery_extrapolate_Quadrilateral(int degree, double* outval, double* vals2project);

void gausspoints_extrapolate_Tetrahedron(int degree, vector<double>& gpts1, vector<double>& gpts2, vector<double>& gpts3);

void stressrecovery_extrapolate_Tetrahedron(int degree, double* outval, double* vals2project);

void gausspoints_extrapolate_Wedge(int degree, vector<double>& gpts1, vector<double>& gpts2, vector<double>& gpts3);

void stressrecovery_extrapolate_Wedge(int degree, double* outval, double* vals2project);

void gausspoints_extrapolate_Hexahedron(int degree, vector<double>& gpts1, vector<double>& gpts2, vector<double>& gpts3);

void stressrecovery_extrapolate_Hexahedron(int degree, double* outval, double* vals2project);



#endif

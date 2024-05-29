#ifndef incl_growthfunctions_h
#define incl_growthfunctions_h


#include "headersEigen.h"


int  growthfunction_isotropic(int  ftype, int  ndim,  double  growthFactor,  MatrixXd&  Fg);

// growth function for the plane-strain beam with a single layer
// for comparison with numerical examples presented in
// Wang et al. IJNLM, 106:280–287, 2018
//
int  growthfunction_beam_singlelayer(double  timeFact, double* geom,  MatrixXd&  Fg);


int  growthfunction_beam_singlelayer3D(double  timeFact, double* geom,  MatrixXd&  Fg);


int  growthfunction_circularplate_singlelayer(double  timeFact, double* geom,  MatrixXd&  Fg);


// growth function for the plane-strain beam with multiple layers
// for comparison with numerical examples presented in
// Du et al. IJNLM, 119:103370, 2020
//
int  growthfunction_beam_multilayer(double  timeFact, double* geom,  MatrixXd&  Fg);


int  growthfunction_anisotropic(int  ftype, int  ndim,  double  growthFactor,  MatrixXd&  Fg);


int  growthfunction_rod(int ftype, int ndim, double timeFact, double* geom,  Eigen::MatrixXd& Fg);



#endif

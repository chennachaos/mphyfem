
#ifndef incl_utilities_material_hyperelastic_h
#define incl_utilities_material_hyperelastic_h

#include <vector>
#include <cmath>
#include "headersEigen.h"
#include "InternalVariables.h"


using std::vector;


// Compute the first derivative of volumetric energy function
double  volumetricfunctions_getFirstDerivative(int type, double K, double J);


// Compute the second derivative of volumetric energy function
double  volumetricfunctions_getSecondDerivative(int type, double K, double J);


// Compute the constants in the mixed formulation
void  compute_constants_volumetric_functions(bool finite, int Utype, double Jn, double BULK, double& Jhat, double& thetahat);


// Compute stress and tangent for the Saint Venant-Kirchhoff model
int computeStressAndTangent_SaintVenantKirchhoff(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat);


// Compute stress and tangent for the Neo-Hookean model
int computeStressAndTangent_NeoHooke(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat);


// Compute stress and tangent for the Mooney-Rivlin model
int computeStressAndTangent_MooneyRivlin(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat);


// Compute stress and tangent for the Gent model
int computeStressAndTangent_Gent(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat);


// Compute stress and tangent for the Arruda-Boyce model
int computeStressAndTangent_ArrudaBoyce(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat);


// Compute stress and tangent for the Eight-chain Longevin model
int computeStressAndTangent_Longevin8chain(int sss, bool MIXED_ELEMENT, int Utype, double Kinv, vector<double>& matData, MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat);


// Compute the stress and tangent tensors for the viscoelastic part
int computeStressAndTangent_Viscoelasticity_Model1(vector<double>& data_viscoelastic, int gp, MatrixXd& F, VectorXd& td, InternalVariables& ivar, double dt, VectorXd&  stre, MatrixXd&  Cmat);




#endif

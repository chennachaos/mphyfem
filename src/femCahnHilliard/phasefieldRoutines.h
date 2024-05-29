#ifndef MY_PHASEFIELD_ROUTINES_H
#define MY_PHASEFIELD_ROUTINES_H



double  chemicalPotential(int type, double WellHeight, double phi);


double  chemicalPotentialDerivative(int type, double WellHeight, double phi);


double  mobilityFunction(int  type, double constD, double phi);


double  mobilityFunctionDerivative(int  type, double constD, double phi);






#endif


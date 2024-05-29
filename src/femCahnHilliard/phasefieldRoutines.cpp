#include "phasefieldRoutines.h"
#include <cmath>
#include <iostream>




double  chemicalPotential(int type, double WellHeight, double phi)
{
    return (4.0*WellHeight*phi*(phi*phi-1.0) );
}




double  chemicalPotentialDerivative(int type, double WellHeight, double phi)
{
    return (4.0*WellHeight*(3.0*phi*phi-1.0) );
}




double  mobilityFunction(int type, double constD, double phi)
{
    double  value = 0.0;

    switch(type)
    {
        case 1:

            value = constD;

        break;

        case 2:
            
            if(abs(phi) <= 1.0)
              value = constD*(1.0-phi*phi);
            else
              value = 0.0;

        break;

        case 3:

            if( phi >= -1.0)
              value = constD*(-2.0*phi*phi*phi-3.0*phi*phi+1.0);
            else if( phi <= 1.0)
              value = constD*( 2.0*phi*phi*phi-3.0*phi*phi+1.0);
            else
              value = 0.0;

        break;

        default:
            std::cerr << " mobilityFunction ... invalid value of type!" << std::endl;

        break;
    }
    
    return  value;
}







double mobilityFunctionDerivative(int type, double constD, double phi)
{
    double  value = 0.0;

    switch(type)
    {
        case 1:

            value = 0.0;

        break;

        case 2:
            
            if(abs(phi) <= 1.0)
              value = constD*(-2.0*phi);
            else
              value = 0.0;

        break;

        case 3:

            if( phi >= -1.0)
              value = constD*(-6.0*phi*phi-6.0*phi);
            else if( phi <= 1.0)
              value = constD*( 6.0*phi*phi-6.0*phi);
            else
              value = 0.0;

        break;

        default:
            std::cerr << " mobilityFunctionDerivative ... invalid value of type!" << std::endl;

        break;
    }

    return  value;
}




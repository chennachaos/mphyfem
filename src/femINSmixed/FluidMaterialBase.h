#ifndef incl_FluidMaterialBase_h
#define incl_FluidMaterialBase_h

#include <vector>
#include "headersEigen.h"
#include <iostream>

using namespace std;


using Eigen::VectorXd;
using Eigen::MatrixXd;





class FluidMaterialBase
{
  public:

    //member variables
    int  id;
    double  rho, mu, nu;

    vector<double>  matData;

    //SolutionData  *SolnData;

    //member functions

    FluidMaterialBase();

    virtual ~FluidMaterialBase();

    virtual int getMaterialTypeNameNumber()
    {  return -1; }

    void setID(int val)
    {
        id = val;
        return;
    }

    double  getDensity()
    {
        return rho;
    }

    double  getDynamicViscosity()
    {
        return mu;
    }

    double  getKinematicViscosity()
    {
        return (mu/rho);
    }


    virtual int readInput(ifstream& infile, string& line);

    virtual void printData();

};

#endif


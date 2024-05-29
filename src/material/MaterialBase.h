#ifndef incl_MaterialBase_h
#define incl_MaterialBase_h

#include <vector>
#include "headersBasic.h"
#include "headersEigen.h"
#include "InternalVariables.h"
#include "SolutionDataSolid.h"
#include <iostream>

using namespace std;


using Eigen::VectorXd;
using Eigen::MatrixXd;





class MaterialBase
{
  public:

    //member variables
    bool MIXED_ELEMENT;
    int  id, Utype;
    string  tis;
    double  rho0, Kinv, spectralRadius;
    vector<double>  matData;
    vector<double>  data_deviatoric, data_viscoelastic, data_magnfield, data_elecfield;
    VectorXd  ResiMagnfield, ApplMagnfield, td;

    SolutionDataSolid  *SolnData;

    //member functions

    MaterialBase();

    virtual ~MaterialBase();

    virtual int getMaterialTypeNameNumber()
    {  return -1; }

    void setID(int val)
    {
        id = val;
        return;
    }

    int  getUtype()
    {
        return Utype;
    }

    virtual  double  getKinv()
    {
        return Kinv;
    }

    double  getDensity()
    {
        return rho0;
    }

    VectorXd  getResiMagnfield()
    {
        return ResiMagnfield;
    }

    VectorXd  getApplMagnfield()
    {
        return ApplMagnfield;
    }

    virtual bool isFiniteStrain()
    {
        return true;
    }

    virtual int readInput(ifstream& infile, string& line);

    virtual void printData();

    virtual double  getPermittivity()
    { cout << "   'getPermittivity' is not defined for this material!\n\n"; return -1; }

    virtual double computeValue(int sss,  MatrixXd&  F)
    { cout << "   'computeValue' is not defined for this material!\n\n"; return -1; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat)
    { cout << "   'computeStressAndTangent1' is not defined for this material!\n\n"; return -1; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, double& pres, VectorXd&  stre, MatrixXd&  Cmat, InternalVariables& ivar, int gp, double dt)
    { cout << "   'computeStressAndTangent2' is not defined for this material!\n\n"; return -1; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt)
    { cout << "   'computeStressAndTangent3' is not defined for this material!\n\n"; return -1; }

    virtual int computeStressAndTangent(bool tangFlag, int sss,  MatrixXd&  Fn, MatrixXd&  F, VectorXd&  elecField, double& temp, double& pres, VectorXd&  stre, VectorXd&  elecDisp, MatrixXd&  Cmat, MatrixXd&  Bmat, MatrixXd&  Amat, InternalVariables& ivar, int gp, double dt)
    { cout << "   'computeStressAndTangent4' is not defined for this material!\n\n"; return -1; }

    virtual int computeElectricComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat)
    { cout << "   'computeElectricComponents' is not defined for this material!\n\n"; return -1; }

    virtual int computeElectricDisplacement(VectorXd&  elecField, VectorXd&  elecDisp)
    { cout << "   'computeElectricDisplacement' is not defined for this material!\n\n"; return -1; }

    virtual int computeMagneticComponents(VectorXd&  elecField, VectorXd&  elecDisp, MatrixXd&  Amat)
    { cout << "   'computeMagneticComponents' is not defined for this material!\n\n"; return -1; }

    virtual int computeMagneticDisplacement(VectorXd&  elecField, VectorXd&  elecDisp)
    { cout << "   'computeMagneticDisplacement' is not defined for this material!\n\n"; return -1; }

    virtual int computeMechanicalStress(MatrixXd&  F, double&  pres, VectorXd&  stre)
    { cout << "   'computeMechanicalStress' is not defined for this material!\n\n"; return -1; }

    virtual int computeMaxwellStress(VectorXd&  elecField, VectorXd&  stre)
    { cout << "   'computeMaxwellStress' is not defined for this material!\n\n"; return -1; }

    virtual int computeTotalStress(MatrixXd&  F, VectorXd&  elecField, double& pres, VectorXd&  stre)
    { cout << "   'computeTotalStress' is not defined for this material!\n\n"; return -1; }

    virtual int getStressAndJacobianForInverseGrowth(int sss,  MatrixXd&  F, MatrixXd&  Fg, VectorXd&  stre, MatrixXd&  Jgp, MatrixXd&  Hgp)
    { cout << "   'computeStressAndTangent1' is not defined for this material!\n\n"; return -1; }
};

#endif


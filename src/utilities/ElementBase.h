#ifndef incl_LagrangeElement_h
#define incl_LagrangeElement_h

#include <vector>
#include "SolverEigen.h"
#include "SolverPetsc.h"
#include "MaterialBase.h"
#include "InternalVariables.h"
#include "GeomDataLagrange.h"
#include "SolutionDataSolid.h"
#include "ElementTypeData.h"
#include "myIDmaps.h"
#include "BasisFunctionsLagrange.h"
#include "BasisFunctionsBernstein.h"

using std::vector;
using std::cout;
using std::endl;

using Eigen::VectorXd;
using Eigen::MatrixXd;





class ElementBase
{
  private:
    static int  elemcount;

  public:

    //member variables

    bool finite, axsy;
    ElementShape  ELEM_SHAPE;

    // sss=stressStrainState
    // 1.) plane stress 2.) plane strain 3.) axisymmetric

    int  elmType, matType, secType, processorId;
    int  matId, finiteInt, sss, degree, npElem, ndim;
    int  ndof, nsize, nivGP, nGP, elenum, nlbfU, nlbfF, nlbfP, nlbfT, nlbfC, AlgoType;

    double  elemError, elemVol, elemVolCur, elemVolOrig, pres;

    vector<double>  vals2project;
    InternalVariables  ivar;

    vector<int>  nodeNums, nodeNumsPres, forAssyVec, forAssyVec2, forAssyVecPres;
    vector<int>  forAssyVecCpot, forAssyVecMpot, forAssyVecEpot, forAssyVecTemp, globalDOFnums;

    SolutionDataSolid  *SolnData;
    GeomDataLagrange  *GeomData;
    MaterialBase  *MatlData;
    ElementTypeData  *ElemTypeData;

    // entries for Kpu and Kup matrices for mixed formulation 
    // to compute pressure variable
    VectorXd  Nc, dNc_dx, dNc_dy, dNc_dz, FlocalPres;
    VectorXd  dNc_dx2, dNc_dy2, dNc_dz2;
    VectorXd  presDOF, presDOFprev, presDOFprev2;

    MatrixXd  Kup, Kpu, Kpp;

    //member functions

    ElementBase();

    virtual ~ElementBase();

    int getDimension()
    { return ndim; }

    int getPolynomialDegree()
    { return degree; }

    void  setSubdomainId(int sid)
    {  processorId = sid; return;  }

    int getSubdomainId()
    {  return  processorId;  }

    int getNodesPerElement()
    {  return  npElem; }

    int getNdofPerNode()
    { return ndof; }

    int  getNdofPerElement()
    {  return  nsize;  }

    int  getElementTypeNumber()
    {  return  elmType;  }

    int  getMaterialTypeNumber()
    {  return  matType;  }

    int  getSectionTypeNumber()
    {  return  secType;  }

    std::vector<int>&  getNodeNumbers()
    {  return  nodeNums; }

    std::vector<int>&  getVectorForAssembly()
    {  return  forAssyVec; }

    virtual int getElmTypeNameNum()
    {  return -1;     }

    virtual void prepareElemData();

    double getError()
    { return  elemError;   }

    double getVolume()
    { return  elemVol;  }

    virtual void initialiseDOFvalues()
    { cout << "   'initialiseDOFvalues' is not defined for this element!\n\n"; return; }

    virtual int calcOutput(double u1, double v1)
    { cout << "   'calcOutput' is not defined for this element!\n\n"; return 0; }

    virtual void setnivGP();

    virtual void initialiseIntVar();
    //{ cout << "   'initialiseIntVar' is not defined for this element!\n\n"; }

    virtual void createTractionDataVariable()
    { cout << "  'createTractionDataVariable' is not available for this element!\n\n"; return; }

    virtual void diffStiffTest(double,int,int,bool);
    //{ cout << "   'diffStiffTest' is not defined for this element!\n\n"; return; }

    virtual int calcInternalForces()
    { cout << "   'calcAndAssyIntForceVec' is not defined for this element!\n\n"; return 0; }

    virtual int calcLoadVector(VectorXd& Flocal)
    { cout << "   'calcLoadVector' is not defined for this element!\n\n"; return 0; }

    virtual int calcForceVectorGrowthModel(VectorXd& Flocal1, VectorXd& Flocal2)
    { cout << "   'calcForceVectorGrowthModel' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidual(MatrixXd& Klocal, VectorXd& Flocal)
    { cout << "   'calcStiffnessAndResidual' is not defined for this element!\n\n"; return 0; }

    virtual int  calcMassMatrix(MatrixXd& Mlocal, bool MassLumping)
    { cout << "   'calcMassMatrix' is not defined for this element!\n\n"; return 0; }

    virtual int  calcMassMatrixPressure(MatrixXd& Mlocal, bool MassLumping)
    { cout << "   'calcMassMatrixPressure' is not defined for this element!\n\n"; return 0; }

    virtual int  calcResidual(VectorXd& Flocal)
    { cout << "   'calcResidual' is not defined for this element!\n\n"; return 0; }

    virtual int  calcResidualPressure(VectorXd& Flocal)
    { cout << "   'calcResidualPressure' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidual(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalF)
    { cout << "   'calcStiffnessAndResidualMixed' is not defined for this element!\n\n"; return 0; }

    int  calcStiffnessAndResidualMixed2D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalP);

    int  calcStiffnessAndResidualMixed3D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalP);

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalP)
    { cout << "   'calcStiffnessAndResidualMixed' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidualMixed2(MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2)
    { cout << "   'calcStiffnessAndResidualMixed2' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidualMixed(MatrixXd& Kuu, MatrixXd& Kuf, MatrixXd& Kfu, MatrixXd& Kff, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& FlocalU, VectorXd& FlocalF, VectorXd& FlocalP)
    { cout << "   'calcStiffnessAndResidualMixed' is not defined for this element!\n\n"; return 0; }

    virtual int  calcStiffnessAndResidualElecField(MatrixXd& Kff, VectorXd& FlocalF)
    { cout << "   'calcStiffnessAndResidualElecField' is not defined for this element!\n\n"; return 0; }

    virtual double  calcCriticalTimeStep(bool flag)
    { cout << "   'calcCriticalTimeStep' is not defined for this element!\n\n"; return 0; }

    virtual int  solveForPressure()
    { cout << "   'solveForPressure' is not defined for this element!\n\n"; return 0; }

    virtual int  solveForPressureMixed(VectorXd& matK1, double beta, double dt, VectorXd& Flocal)
    { cout << "   'solveForPressureMixed' is not defined for this element!\n\n"; return 0; }

    virtual int  toComputeInfSupCondition(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);
    //{ cout << "   'toComputeInfSupCondition' is not defined for this element!\n\n"; return 0; }

    int  toComputeInfSupCondition2D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);
    int  toComputeInfSupCondition3D(MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpp);

    virtual int  applyDirichletBCs(MatrixXd& Klocal, VectorXd& Flocal);

    virtual int  applyDirichletBCsMixed(int offset, MatrixXd& Kuu, MatrixXd& Kup, MatrixXd& Kpu, MatrixXd& Kpp, VectorXd& Flocal1, VectorXd& Flocal2);

    virtual int  applyDirichletBCs2field(int var3_offset, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, VectorXd& Flocal1, VectorXd& Flocal2, vector<int>& forAssyVar1, vector<int>& forAssyVar2, VectorXd& var1applied, VectorXd& var2applied);

    virtual int  applyDirichletBCs3field(int schemeType, MatrixXd& K11, MatrixXd& K12, MatrixXd& K21, MatrixXd& K22, MatrixXd& K13, MatrixXd& K31, MatrixXd& K33, VectorXd& Flocal1, VectorXd& Flocal2, VectorXd& Flocal3, vector<int>& forAssyVar1, vector<int>& forAssyVar2, vector<int>& forAssyVar3, VectorXd& var1applied, VectorXd& var2applied, VectorXd& var3applied);

    virtual int  applyDirichletBCsElecField(MatrixXd& Kff, VectorXd& FlocalF);


    virtual void contourplot(int, int, double, double)
      { cout << "   'contourPlot' is not defined for this element!\n\n"; return; }

    virtual void elementContourplot(int, int, int)
      { cout << "   'elementContourplot' is not defined for this element!\n\n"; return; }

    virtual double computeVolume(bool init = false)
      { cout << "   'computeVolume' is not defined for this element!\n\n"; return 0.; }

    virtual void plotGaussPoints(int, bool defFlg = false)
      { cout << "  'plotGaussPoints' is not available for this element!\n\n"; return; }

    virtual void projectToNodes(bool, int, int, int)
      { cout << "  'projectToNodes' is not available for this element!\n\n"; return; }

    virtual void projectInternalVariable(bool, int, int, int, double*)
      { cout << "  'projectInternalVariable' is not available for this element!\n\n"; return; }

    virtual void projectStrain(bool, int, int, int, double*)
      { cout << "  'projectStrain' is not available for this element!\n\n"; return; }

    virtual void projectStress(bool, int, int, int, double*)
      { cout << "  'projectStress' is not available for this element!\n\n"; return; }

    virtual  void computeEnergy(int, int, VectorXd& )
      { cout << "  'computeEnergy' is not available for this element!\n\n"; return; }

    void computeMomentum(int, int, VectorXd&);

    virtual  void  MatrixToPostprocess(MatrixXd& Klocal)
      { cout << "  'MatrixToPostprocess' is not available for this element!\n\n"; return; }

    virtual  void  RhsToPostprocess(int, int, int, VectorXd& Flocal)
      { cout << "  'RhsToPostprocess' is not available for this element!\n\n"; return; }

    virtual  int calcError(int index, double* val);

    virtual  int calcError2D(int index, double* val);

    virtual  int calcError3D(int index, double* val);

    double  computeGeomOrig(int dir, VectorXd& NN);
    double  computeGeomNew(int dir, VectorXd& NN);
    double  computeGeomCur(int dir, VectorXd& NN);

    void  computeGeomNew(double* param, double* geom);

    double  computeDisplacement(int dir, VectorXd& NN);
    double  computeVelocity(int dir, VectorXd& NN);
    double  computeAcceleration(int dir, VectorXd& NN);

    double  computeDisplacementCur(int dir, VectorXd& NN);
    double  computeVelocityCur(int dir, VectorXd& NN);
    double  computeAccelerationCur(int dir, VectorXd& NN);

    double  computeValue(int dir, VectorXd& NN);
    double  computeValueIncr(int dir, VectorXd& NN);
    double  computeValuePrev(int dir, VectorXd& NN);

    double  computeValueCur(int dir, VectorXd& NN);
    double  computeValueDot(int dir, VectorXd& NN);
    double  computeValueDotCur(int dir, VectorXd& NN);

    double  computeValue2(int dir, VectorXd& NN);
    double  computeValue2Prev(int dir, VectorXd& NN);
    double  computeValue2Cur(int dir, VectorXd& NN);

    double  computeValue3(int dir, VectorXd& NN);
    double  computeValue3Prev(int dir, VectorXd& NN);
    double  computeValue3Cur(int dir, VectorXd& NN);

    double  computeValue4(int dir, VectorXd& NN);
    double  computeValue4Prev(int dir, VectorXd& NN);
    double  computeValue4Cur(int dir, VectorXd& NN);
    double  computeValue4Dot(int dir, VectorXd& NN);
    double  computeValue4DotCur(int dir, VectorXd& NN);

    void  computeDefGrad2D(VectorXd& dN_dx, VectorXd& dN_dy, MatrixXd& F);
    void  computeDefGrad2DPrev(VectorXd& dN_dx, VectorXd& dN_dy, MatrixXd& F);
    void  computeDefGrad2DCur(VectorXd& dN_dx, VectorXd& dN_dy, MatrixXd& F);

    void  computeDefGrad(VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz, MatrixXd& F);
    void  computeDefGradPrev(VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz, MatrixXd& F);
    void  computeDefGradCur(VectorXd& dN_dx, VectorXd& dN_dy, VectorXd& dN_dz, MatrixXd& F);

    double  computeForce(int dir, VectorXd& NN);
    double  computeForcePrev(int dir, VectorXd& NN);
    double  computeForceCur(int dir, VectorXd& NN);
};

#endif //incl_Lagrange_Element_h


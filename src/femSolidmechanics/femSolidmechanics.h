#ifndef incl_femSolidmechanics_CLASS_h
#define incl_femSolidmechanics_CLASS_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "headersVTK.h"
#include "femBase.h"

#include "ElementBase.h"


class femSolidmechanics : public femBase
{
    //private:

    public:

        femSolidmechanics();

        ~femSolidmechanics();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        int  setSolverDataForFullyImplicit();

        int  solveStep(int algotype, int iter_max, int subiter_max, double tol1);

        int  solveFullyImplicit();

        int  solveWithNewtonRaphson();

        int  solveWithArclength();

        int  solveExplicitStep(int*);

        int  solveStepNewtonRaphson(int iter_max, double tol1);

        int  solveStepArcLength(int iter_max, double tol1);

        int  solveExplicitStepNIC(int*);
        int  solveExplicitStepFIC(int*);

        int  saveSolution();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  printComputerTime(bool reset = true, int detailFlg = 1);

        int  setSpecifiedDOFs_Displacement(vector<vector<bool> >& NodeDofType);

        int  setBoundaryConditions();

        int  prepareInputData();

        int  updateGeometry();

        int  printInfo();

        int  writeNodalData();

        int  writeNodalDataExplicitScheme();

        int  writeReadResult(int, string &);

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////


        virtual int setSolver(int);

        int  prepareMatrixPattern();

        int  copyElemInternalVariables();

        int  calcMassMatrixForExplicitDynamics();

        virtual int calcStiffnessAndResidual();

        virtual int factoriseSolveAndUpdate();

        virtual int elementDiffStiffTest(double ddd, int elnum, int dig, int dig2, bool gfrmt);

        virtual bool converged();

        virtual bool diverging(double);

        virtual int timeUpdate();

        virtual int updateIterStep();

        virtual int reset();

        int  addExternalForces();

        int  applyBoundaryConditions();

        int  computeElementErrors(int);

        int  setInitialConditions();

        int  InfSupNumber();

        int  ModalAnalysis(int nn=10, bool flag=true, double fact=1.0);

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  projectFromElemsToNodes(bool, int, int, int);

        void  projectStrains(bool, int, int, int);

        void  projectStresses(bool, int, int, int);

        void  projectInternalVariables(bool, int, int, int);

        virtual  void  plotGeom();

        virtual  void  postProcess();

        void  postProcessPressureForDiscontinuousElements();

        int  processForBernsteinElements();

        int  setSolverDataForTIC();

        int  prepareDataForMixedElements();
};





#endif







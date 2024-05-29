#ifndef incl_femMagnetomech_CLASS_h
#define incl_femMagnetomech_CLASS_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "headersVTK.h"
#include "ElementBase.h"
#include "MaterialBase.h"
#include "femBase.h"


class femMagnetomech : public femBase
{
    //private:

    public:

        femMagnetomech();

        ~femMagnetomech();

        ///////////////////////////////////////////////////////////
        //
        // DATA related member functions
        //
        ///////////////////////////////////////////////////////////

        void  printLogo();

        int  setBoundaryConditions();

        ///////////////////////////////////////////////////////////
        //
        // PRE-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        int  printComputerTime(bool reset = true, int detailFlg = 1);

        int  setSpecifiedDOFs_Displacement(vector<vector<bool> >& NodeDofType);

        int  prepareInputData();

        int  updateGeometry();

        int  processForBernsteinElements();

        int  setSolverDataForSemiImplicit();

        int  prepareDataForMagneticPotential();

        int  prepareDataForTemperature();

        int  prepareDataForPressure();

        int  prepareDataForMixedElements();

        ///////////////////////////////////////////////////////////
        //
        // SOLUTION PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        virtual int setSolver(int);

        int  prepareMatrixPattern();

        int  solveFullyImplicit();

        int  solveWithNewtonRaphson();

        int  solveWithArclength();

        int  copyElemInternalVariables();

        virtual int calcStiffnessAndResidual();

        virtual int factoriseSolveAndUpdate();

        virtual bool converged();

        virtual bool diverging(double);

        virtual int timeUpdate();

        virtual int updateIterStep();

        virtual int saveSolution();

        virtual int reset();

        virtual int addExternalForces();

        int  setInitialConditions();

        int  assignBoundaryConditions();

        ///////////////////////////////////////////////////////////
        //
        // POST-PROCESSOR PHASE member functions
        //
        ///////////////////////////////////////////////////////////

        void  writeNodalData();

        void  writeReadResult(int, string &);

        void  plotGeom();

        void  postProcess();

        void  projectFromElemsToNodes(bool, int, int, int);

        void  projectStrains(bool, int, int, int, VectorXd&  output);

        void  projectStresses(bool, int, int, int, VectorXd&  output);

        void  projectInternalVariables(bool, int, int, int, VectorXd&  output);

};





#endif







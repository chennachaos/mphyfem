#ifndef incl_femSolidmechanics_CLASS_h
#define incl_femSolidmechanics_CLASS_h


#include "util.h"
#include "GeomDataLagrange.h"
#include "SolutionData.h"
#include "headersVTK.h"
#include "femBase.h"

#include "ElementBase.h"


class PropertyItem;


#define LAGRANGE_ELEMENT_TYPE_NAMES {"LagrangeElem1Dbar2Node",                  /*  0 */ \
                                    "EulerBernoulliBeamElement2D",              /*  1 */ \
                                    "FrameElement2D",                           /*  2 */ \
                                    "ElementGeomExactTruss2D",                  /*  3 */ \
                                    "ElementGeomExactBeam2D",                   /*  4 */ \
                                    "LagrangeElem2DPoissonTria3Node",           /*  5 */ \
                                    "LagrangeElem2DPoissonQuad4Node",           /*  6 */ \
                                    "LagrangeElem2DSolidTria3Node",             /*  7 */ \
                                    "LagrangeElem2DSolidQuad4Node",             /*  8 */ \
                                    "LagrangeElem2DSolidQuad4NodeMixed",        /*  9 */ \
                                    "LagrangeElem2DSolidBbarFbarQuad4Node",     /* 10 */ \
                                    "LagrangeElem2DStokesTria3Node",            /* 11 */ \
                                    "LagrangeElem2DStokesQuad4Node",            /* 12 */ \
                                    "LagrangeElem2DNavierStokesTria3Node",      /* 13 */ \
                                    "LagrangeElem2DNavierStokesQuad4Node",      /* 14 */ \
                                    "EulerBernoulliBeamElement3D",              /* 15 */ \
                                    "LagrangeElem3DPoissonTetra4Node",          /* 16 */ \
                                    "LagrangeElem3DPoissonHexa8Node",           /* 17 */ \
                                    "LagrangeElem3DSolidTetra4Node",            /* 18 */ \
                                    "LagrangeElem3DSolidHexa8Node",             /* 19 */ \
                                    "LagrangeElem3DSolidHexa8NodeMixed",        /* 20 */ \
                                    "LagrangeElem3DSolidBbarFbarHexa8Node",     /* 21 */ \
                                    "LagrangeElem3DStokesTetra4Node",           /* 22 */ \
                                    "LagrangeElem3DStokesHexa8Node",            /* 23 */ \
                                    "LagrangeElem3DNavierStokesTetra4Node",     /* 24 */ \
                                    "LagrangeElem3DNavierStokesHexa8Node",      /* 25 */ \
                                    "LagrangeElem2DSolidTria3NodeStab",         /* 26 */ \
                                    "LagrangeElem2DSolidQuad4NodeStab",         /* 27 */ \
                                    "LagrangeElem3DSolidTetra4NodeStab",        /* 28 */ \
                                    "LagrangeElem3DSolidHexa8NodeStab",         /* 29 */ \
                                    "MindlinPlateElement",                      /* 30 */ \
                                    "KirchhoffPlateElement",                    /* 31 */ \
                                    "LagrangeElem3DShellQuad4Node",             /* 32 */ \
                                    "ContactElement2D1nodedContactAlongXaxis",  /* 33 */ \
                                    "ContactElement2D1nodedContactAlongYaxis",  /* 34 */ \
                                    "ContactElement3D1nodedContactAlongXaxis",  /* 35 */ \
                                    "ContactElement3D1nodedContactAlongYaxis",  /* 36 */ \
                                    "ContactElement3D1nodedContactAlongZaxis",  /* 37 */ \
                                    "LagrangeElem3DSolidPyramid5Node",          /* 38 */ \
                                    "LagrangeElem3DSolidPyramid5NodeStab",      /* 39 */ \
                                    "LagrangeElem3DSolidPrism6Node",            /* 40 */ \
                                    "LagrangeElem3DSolidPrism6NodeStab",        /* 41 */ \
                                    "LagrangeElem2DSolidTria6Node",             /* 42 */ \
                                    "LagrangeElem2DSolidQuad9Node",             /* 43 */ \
                                    "LagrangeElem2DSolidBbarFbarTria6Node",     /* 44 */ \
                                    "LagrangeElem2DSolidBbarFbarQuad9Node",     /* 45 */ \
                                    "LagrangeElem3DSolidTetra10Node",           /* 46 */ \
                                    "LagrangeElem3DSolidPyramid14Node",         /* 47 */ \
                                    "LagrangeElem3DSolidPrism18Node",           /* 48 */ \
                                    "LagrangeElem3DSolidHexa27Node",            /* 49 */ \
                                    "LagrangeElem3DSolidBbarFbarTetra10Node",   /* 50 */ \
                                    "LagrangeElem3DSolidBbarFbarPyramid14Node", /* 51 */ \
                                    "LagrangeElem3DSolidBbarFbarPrism18Node",   /* 52 */ \
                                    "LagrangeElem3DSolidBbarFbarHexa27Node",    /* 53 */ \
                                    "LagrangeElem2DPoissonTria6Node",           /* 54 */ \
                                    "LagrangeElem3DPoissonTetra10Node",         /* 55 */ \
                                    "LagrangeElem2DSolidTria3NodeStabSI",       /* 56 */ \
                                    "LagrangeElem2DSolidQuad4NodeStabSI",       /* 57 */ \
                                    "LagrangeElem3DSolidTetra4NodeStabSI",      /* 58 */ \
                                    "LagrangeElem3DSolidHexa8NodeStabSI", NULL};



#define BERNSTEIN_ELEMENT_TYPE_NAMES {"BernsteinElem2DPoissonTria6Node",         /*  201 */ \
                                      "BernsteinElem2DSolidTria6Node",           /*  202 */ \
                                      "BernsteinElem2DSolidBbarFbarTria6Node",   /*  203 */ \
                                      "BernsteinElem2DSolidMixedTria6Node",      /*  204 */ \
                                      "BernsteinElem3DPoissonTetra10Node",       /*  205 */ \
                                      "BernsteinElem3DSolidTetra10Node",         /*  206 */ \
                                      "BernsteinElem3DSolidBbarFbarTetra10Node", /*  207 */ \
                                      "BernsteinElem3DSolidMixedTetra10Node",    /*  208 */ \
                                      "dummy",                                   /*  209 */ \
                                      "dummy",                                   /*  210 */ \
                                      "BernsteinElem2DSolidMixedTria6NodeSI1",   /*  211 */ \
                                      "BernsteinElem2DSolidMixedTria6NodeSI3",   /*  212 */ \
                                      "BernsteinElem3DSolidMixedTetra10NodeSI1", /*  213 */ \
                                      "BernsteinElem3DSolidMixedTetra10NodeSI4", /*  214 */ \
                                      "BernsteinElem2DSolidMixedTria6Node2",     /*  215 */ \
                                      "BernsteinElem3DSolidMixedTetra10Node2",   /*  216 */ \
                                      "dummy",                                   /*  217 */ \
                                      "dummy",                                   /*  218 */ \
                                      "dummy",                                   /*  219 */ \
                                      "dummy",                                   /*  220 */ \
                                      "dummy",                                   /*  221 */ \
                                      "BernsteinElem2DSolidQuad9Node",           /*  222 */ \
                                      "BernsteinElem2DSolidMixedQuad9NodeSI4",   /*  223 */ \
                                      "BernsteinElem3DSolidWedge18Node",         /*  224 */ \
                                      "BernsteinElem3DSolidMixedWedge18NodeSI6", /*  225 */ \
                                      "BernsteinElem3DSolidHexa27Node",          /*  226 */ \
                                      "BernsteinElem3DSolidMixedHexa27NodeSI8",  NULL};


#define LAGRANGE_ELEMENT_TYPE_NAMES_FACE {\
                                     "LagrangeElem2DEdge2Node", \
                                     "LagrangeElem2DEdge3Node", \
                                     "LagrangeElem3DFaceTria3Node", \
                                     "LagrangeElem3DFaceTria6Node", \
                                     "LagrangeElem3DFaceQuad4Node", \
                                     "LagrangeElem3DFaceQuad9Node", NULL};


#define BERNSTEIN_ELEMENT_TYPE_NAMES_FACE {\
                                      "BernsteinElem2DEdge3Node", \
                                      "BernsteinElem3DFaceTria6Node", \
                                      "BernsteinElem3DFaceQuad9Node", NULL};



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

        virtual int calcStiffnessAndResidual(int printRes=2, bool zeroMtx=true, bool zeroRes=true);

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

        int  processForBernsteinElements();

        int  setSolverDataForTIC();

        int  prepareDataForMixedElements();
};





#endif







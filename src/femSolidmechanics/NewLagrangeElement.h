
#ifndef  incl_newLagrangeElement_h
#define  incl_newLagrangeElement_h

#include "ElementBase.h"
#include "myIDmaps.h"

#include "Elem_Solid_2D_Disp.h"
#include "Elem_Solid_2D_Mixed0.h"
#include "Elem_Solid_2D_Mixed1.h"
#include "Elem_Solid_2D_Tria7_Mixed1.h"
#include "Elem_Solid_2D_Tria10_Mixed1.h"


#include "Elem_Solid_3D_Disp.h"
#include "Elem_Solid_3D_Mixed0.h"
#include "Elem_Solid_3D_Mixed1.h"
#include "Elem_Solid_3D_Tetra11_Mixed1.h"
#include "Elem_Solid_3D_Tetra20_Mixed1.h"


inline  ElementBase*  NewLagrangeElement(int type)
{

  switch (type)
  {
    case  ELEM_SOLID_2D_DISP:  return (ElementBase*) new Elem_Solid_2D_Disp; break;

    case  ELEM_SOLID_2D_MIXED0:  return (ElementBase*) new Elem_Solid_2D_Mixed0; break;

    case  ELEM_SOLID_2D_MIXED1:  return (ElementBase*) new Elem_Solid_2D_Mixed1; break;

    case  ELEM_SOLID_2D_TRIA7_MIXED1:  return (ElementBase*) new Elem_Solid_2D_Tria7_Mixed1; break;

    case  ELEM_SOLID_2D_TRIA10_MIXED1:  return (ElementBase*) new Elem_Solid_2D_Tria10_Mixed1; break;


    case  ELEM_SOLID_3D_DISP:  return (ElementBase*) new Elem_Solid_3D_Disp; break;

    case  ELEM_SOLID_3D_MIXED0:  return (ElementBase*) new Elem_Solid_3D_Mixed0; break;

    case  ELEM_SOLID_3D_MIXED1:  return (ElementBase*) new Elem_Solid_3D_Mixed1; break;

    case  ELEM_SOLID_3D_TETRA11_MIXED1:  return (ElementBase*) new Elem_Solid_3D_Tetra11_Mixed1; break;

    case  ELEM_SOLID_3D_TETRA20_MIXED1:  return (ElementBase*) new Elem_Solid_3D_Tetra20_Mixed1; break;

    /*
    case  0: return (ElementBase*) new LagrangeElem1Dbar2Node; break;

    case  1: return (ElementBase*) new EulerBernoulliBeamElement2D; break;

    case  2: return (ElementBase*) new FrameElement2D; break;

    case  3: return (ElementBase*) new ElementGeomExactTruss2D; break;

    case  4: return (ElementBase*) new ElementGeomExactBeam2D; break;

    case  7: return (ElementBase*) new LagrangeElem2DSolidTria3Node; break;

    case  8: return (ElementBase*) new LagrangeElem2DSolidQuad4Node; break;

    case  9: return (ElementBase*) new LagrangeElem2DSolidQuad4NodeMixed; break;

    case 10: return (ElementBase*) new LagrangeElem2DSolidBbarFbarQuad4Node; break;

    case 15: return (ElementBase*) new EulerBernoulliBeamElement3D; break;

    case 18: return (ElementBase*) new LagrangeElem3DSolidTetra4Node; break;

    case 19: return (ElementBase*) new LagrangeElem3DSolidHexa8Node; break;

    case 20: return (ElementBase*) new LagrangeElem3DSolidHexa8NodeMixed; break;

    case 21: return (ElementBase*) new LagrangeElem3DSolidBbarFbarHexa8Node; break;

    case 26: return (ElementBase*) new LagrangeElem2DSolidTria3NodeStab; break;

    case 27: return (ElementBase*) new LagrangeElem2DSolidQuad4NodeStab; break;

    case 28: return (ElementBase*) new LagrangeElem3DSolidTetra4NodeStab; break;

    case 29: return (ElementBase*) new LagrangeElem3DSolidHexa8NodeStab; break;

    case 30: return (ElementBase*) new MindlinPlateElement; break;

    case 31: return (ElementBase*) new KirchhoffPlateElement; break;

    case 32: return (ElementBase*) new LagrangeElem3DShellQuad4Node; break;

    case 33: return (ElementBase*) new ContactElement2D1nodedContactAlongXaxis; break;

    case 34: return (ElementBase*) new ContactElement2D1nodedContactAlongYaxis; break;

    case 35: return (ElementBase*) new ContactElement3D1nodedContactAlongXaxis; break;

    case 36: return (ElementBase*) new ContactElement3D1nodedContactAlongYaxis; break;

    case 37: return (ElementBase*) new ContactElement3D1nodedContactAlongZaxis; break;

    //case 38: return (ElementBase*) new LagrangeElem3DSolidPyramid5Node; break;

    //case 39: return (ElementBase*) new LagrangeElem3DSolidPyramid5NodeStab; break;

    //case 40: return (ElementBase*) new LagrangeElem3DSolidPrism6Node; break;

    //case 41: return (ElementBase*) new LagrangeElem3DSolidPrism6NodeStab; break;

    case 42: return (ElementBase*) new LagrangeElem2DSolidTria6Node; break;

    case 43: return (ElementBase*) new LagrangeElem2DSolidQuad9Node; break;

    case 44: return (ElementBase*) new LagrangeElem2DSolidBbarFbarTria6Node; break;

    case 45: return (ElementBase*) new LagrangeElem2DSolidBbarFbarQuad9Node; break;

    case 46: return (ElementBase*) new LagrangeElem3DSolidTetra10Node; break;

    //case 47: return (ElementBase*) new LagrangeElem3DSolidPyramid14Node; break;

    //case 48: return (ElementBase*) new LagrangeElem3DSolidPrism18Node; break;

    case 49: return (ElementBase*) new LagrangeElem3DSolidHexa27Node; break;

    case 50: return (ElementBase*) new LagrangeElem3DSolidBbarFbarTetra10Node; break;

    //case 51: return (ElementBase*) new LagrangeElem3DSolidBbarFbarPyramid14Node; break;

    //case 52: return (ElementBase*) new LagrangeElem3DSolidBbarFbarPrism18Node; break;

    case 53: return (ElementBase*) new LagrangeElem3DSolidBbarFbarHexa27Node; break;

    case 54: return (ElementBase*) new LagrangeElem2DPoissonTria6Node; break;

    case 55: return (ElementBase*) new LagrangeElem3DPoissonTetra10Node; break;

    case 56: return (ElementBase*) new LagrangeElem2DSolidTria3NodeStabSI; break;

    case 57: return (ElementBase*) new LagrangeElem2DSolidQuad4NodeStabSI; break;

    case 58: return (ElementBase*) new LagrangeElem3DSolidTetra4NodeStabSI; break;

    case 59: return (ElementBase*) new LagrangeElem3DSolidHexa8NodeStabSI; break;

    // Elements based on Bernstein polynomials. Only quadratic elements are considered.
    //case 201: return (ElementBase*) new BernsteinElem2DPoissonTria6Node; break;

    case 202: return (ElementBase*) new BernsteinElem2DSolidTria6Node; break;

    case 203: return (ElementBase*) new BernsteinElem2DSolidBbarFbarTria6Node; break;

    case 204: return (ElementBase*) new BernsteinElem2DSolidMixedTria6Node; break;

    //case 205: return (ElementBase*) new BernsteinElem3DPoissonTetra10Node; break;

    case 206: return (ElementBase*) new BernsteinElem3DSolidTetra10Node; break;

    case 207: return (ElementBase*) new BernsteinElem3DSolidBbarFbarTetra10Node; break;

    case 208: return (ElementBase*) new BernsteinElem3DSolidMixedTetra10Node; break;

    case 211: return (ElementBase*) new BernsteinElem2DSolidMixedTria6NodeSI1; break;

    case 212: return (ElementBase*) new BernsteinElem2DSolidMixedTria6NodeSI3; break;

    case 213: return (ElementBase*) new BernsteinElem3DSolidMixedTetra10NodeSI1; break;

    case 214: return (ElementBase*) new BernsteinElem3DSolidMixedTetra10NodeSI4; break;


    //case 221: return (ElementBase*) new BernsteinElem2DPoissonQuad9Node; break;
    //case 224: return (ElementBase*) new BernsteinElem3DPoissonHex27Node; break;

    case 222: return (ElementBase*) new BernsteinElem2DSolidQuad9Node; break;

    case 223: return (ElementBase*) new BernsteinElem2DSolidMixedQuad9NodeSI4; break;

    case 224: return (ElementBase*) new BernsteinElem3DSolidWedge18Node; break;

    case 225: return (ElementBase*) new BernsteinElem3DSolidMixedWedge18NodeSI6; break;

    case 226: return (ElementBase*) new BernsteinElem3DSolidHexa27Node; break;

    case 227: return (ElementBase*) new BernsteinElem3DSolidMixedHexa27NodeSI8; break;


    // Edge/Face elements for applying Neumann Boundary conditions
    case 401: return (ElementBase*) new LagrangeElem2DEdge2Node; break;

    case 402: return (ElementBase*) new LagrangeElem2DEdge3Node; break;

    case 403: return (ElementBase*) new LagrangeElem3DFaceTria3Node; break;

    case 404: return (ElementBase*) new LagrangeElem3DFaceTria6Node; break;

    case 405: return (ElementBase*) new LagrangeElem3DFaceQuad4Node; break;

    case 406: return (ElementBase*) new LagrangeElem3DFaceQuad9Node; break;



    case 601: return (ElementBase*) new BernsteinElem2DEdge3Node; break;

    case 602: return (ElementBase*) new BernsteinElem3DFaceTria6Node; break;

    case 603: return (ElementBase*) new BernsteinElem3DFaceQuad9Node; break;
    */

    default: cerr << "NewLagrangeElement ... unknown element type" << endl; return NULL;
  }

}


#endif



#ifndef  incl_newElement_h
#define  incl_newElement_h

#include "ElementBase.h"
#include "FunctionsProgram.h"

#include "BernsteinElem2DElecMechTri11.h"                   // 1001
#include "BernsteinElem2DElecMechTri22.h"                   // 1002
#include "BernsteinElem2DElecMechTri22F.h"                  // 1003
#include "BernsteinElem2DElecMechTri210.h"                  // 1004
#include "BernsteinElem2DElecMechTri220.h"                  // 1005
#include "BernsteinElem2DElecMechTri211.h"                  // 1006
#include "BernsteinElem2DElecMechTri221.h"                  // 1007
#include "BernsteinElem2DElecMechQua11.h"                   // 1008
#include "BernsteinElem2DElecMechQua11F.h"                  // 1009
#include "BernsteinElem2DElecMechQua110.h"                  // 1010
#include "BernsteinElem2DElecMechQua1100.h"                 // 1011
#include "BernsteinElem2DElecMechQua221.h"                  // 1012
#include "BernsteinElem2DElecMechQua221P.h"                 // 1013


# include "BernsteinElem3DElecMechTet11.h"                  // 1041
#include "BernsteinElem3DElecMechTet22.h"                   // 1042
#include "BernsteinElem3DElecMechTet22F.h"                  // 1043
//#include "BernsteinElem3DElecMechTet210.h"                  // 1044
#include "BernsteinElem3DElecMechTet220.h"                  // 1045
//#include "BernsteinElem3DElecMechTet211.h"                  // 1046
#include "BernsteinElem3DElecMechTet221.h"                  // 1047
#include "BernsteinElem3DElecMechHex11.h"                   // 1048
#include "BernsteinElem3DElecMechHex11F.h"                  // 1049
#include "BernsteinElem3DElecMechHex110.h"                  // 1050
#include "BernsteinElem3DElecMechHex1100.h"                 // 1051
#include "BernsteinElem3DElecMechHex221.h"                  // 1052
#include "BernsteinElem3DElecMechHex221P.h"                 // 1053


#include "BernsteinElem2DEdge2Node.h"
#include "BernsteinElem2DEdge3Node.h"
#include "BernsteinElem3DFace.h"

/*
#include "TherMechElem2DQua11.h"                   // 2001
#include "TherMechElem2DQua110.h"                  // 2002
#include "TherMechElem2DTri22.h"                   // 2003
#include "TherMechElem2DTri220.h"                  // 2004
#include "TherMechElem2DTri221.h"                  // 2005

#include "TherMechElem3DHex11.h"                   // 2041
#include "TherMechElem3DHex110.h"                  // 2042
#include "TherMechElem3DTet22.h"                   // 2043
#include "TherMechElem3DTet220.h"                  // 2044
#include "TherMechElem3DTet221.h"                  // 2045


#include "TherMechElem2DEdge2Node.h"
#include "TherMechElem2DEdge3Node.h"
#include "TherMechElem3DFaceQua4Node.h"
#include "TherMechElem3DFaceTri6Node.h"
*/


inline  ElementBase*  NewElement(int type)
{
  switch (type)
  {
    // Elements for coupled electromechanical problems

    case 1001: return (ElementBase*) new BernsteinElem2DElecMechTri11;   break;
    case 1002: return (ElementBase*) new BernsteinElem2DElecMechTri22;   break;
    case 1003: return (ElementBase*) new BernsteinElem2DElecMechTri22F;  break;
    case 1004: return (ElementBase*) new BernsteinElem2DElecMechTri210;  break;
    case 1005: return (ElementBase*) new BernsteinElem2DElecMechTri220;  break;
    case 1006: return (ElementBase*) new BernsteinElem2DElecMechTri211;  break;
    case 1007: return (ElementBase*) new BernsteinElem2DElecMechTri221;  break;
    case 1008: return (ElementBase*) new BernsteinElem2DElecMechQua11;   break;
    case 1009: return (ElementBase*) new BernsteinElem2DElecMechQua11F;  break;
    case 1010: return (ElementBase*) new BernsteinElem2DElecMechQua110;  break;
    case 1011: return (ElementBase*) new BernsteinElem2DElecMechQua1100; break;
    case 1012: return (ElementBase*) new BernsteinElem2DElecMechQua221;  break;
    case 1013: return (ElementBase*) new BernsteinElem2DElecMechQua221P; break;


    case 1041: return (ElementBase*) new BernsteinElem3DElecMechTet11;  break;
    case 1042: return (ElementBase*) new BernsteinElem3DElecMechTet22;  break;
    case 1043: return (ElementBase*) new BernsteinElem3DElecMechTet22F; break;
//    case 1044: return (ElementBase*) new BernsteinElem3DElecMechTet210; break;
    case 1045: return (ElementBase*) new BernsteinElem3DElecMechTet220; break;
//    case 1046: return (ElementBase*) new BernsteinElem3DElecMechTet211; break;
    case 1047: return (ElementBase*) new BernsteinElem3DElecMechTet221; break;
    case 1048: return (ElementBase*) new BernsteinElem3DElecMechHex11;  break;
    case 1049: return (ElementBase*) new BernsteinElem3DElecMechHex11F; break;
    case 1050: return (ElementBase*) new BernsteinElem3DElecMechHex110; break;
    case 1051: return (ElementBase*) new BernsteinElem3DElecMechHex1100;break;
    case 1052: return (ElementBase*) new BernsteinElem3DElecMechHex221; break;
    case 1053: return (ElementBase*) new BernsteinElem3DElecMechHex221P; break;

    case 1081: return (ElementBase*) new BernsteinElem2DEdge2Node; break;
    case 1082: return (ElementBase*) new BernsteinElem2DEdge3Node; break;
    case 1083: return (ElementBase*) new BernsteinElem3DFace; break;


/*
    case 2001: return (ElementBase*) new TherMechElem2DQua11;   break;
    case 2002: return (ElementBase*) new TherMechElem2DQua110;  break;
    case 2003: return (ElementBase*) new TherMechElem2DTri22;   break;
    case 2004: return (ElementBase*) new TherMechElem2DTri220;  break;
    case 2005: return (ElementBase*) new TherMechElem2DTri221;  break;

    case 2041: return (ElementBase*) new TherMechElem3DHex11;   break;
    case 2042: return (ElementBase*) new TherMechElem3DHex110;  break;
    case 2043: return (ElementBase*) new TherMechElem3DTet22;   break;
    case 2044: return (ElementBase*) new TherMechElem3DTet220;  break;
    case 2045: return (ElementBase*) new TherMechElem3DTet221;  break;

    case 2081: return (ElementBase*) new TherMechElem2DEdge2Node;  break;
    case 2082: return (ElementBase*) new TherMechElem2DEdge3Node;  break;
    case 2083: return (ElementBase*) new TherMechElem3DFaceTri6Node;  break;
*/

    default: prgError(1,"NewElement","unknown element type name!"); return NULL;
  }
}


#endif


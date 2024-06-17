
#ifndef  incl_newElementGrowthFEM_h
#define  incl_newElementGrowthFEM_h

#include "ElementBase.h"
#include "FunctionsProgram.h"

#include "GrowthElem2DQua1.h"                    // 6001
#include "GrowthElem2DQua10.h"                   // 6002
#include "GrowthElem2DQua21.h"                   // 6003
#include "GrowthElem2DTri21.h"                   // 6004
#include "GrowthElem2DQua21J.h"                  // 6005
#include "GrowthElem2DQua21P.h"                  // 6006


//#include "GrowthElem3DHex1.h"                    // 6051
//#include "GrowthElem3DHex10.h"                   // 6052
#include "Elem_Growth_3D_Mixed21.h"                 // 6053
//#include "GrowthElem3DHex21.h"                   // 6054
//#include "GrowthElem3DTet21.h"                   // 6055



inline  ElementBase*  NewElementGrowthFEM(int type)
{
  switch (type)
  {
    case 6001: return (ElementBase*) new GrowthElem2DQua1;     break;
    case 6002: return (ElementBase*) new GrowthElem2DQua10;    break;
    case 6003: return (ElementBase*) new GrowthElem2DQua21;    break;
    case 6004: return (ElementBase*) new GrowthElem2DTri21;    break;
    case 6005: return (ElementBase*) new GrowthElem2DQua21J;   break;
    case 6006: return (ElementBase*) new GrowthElem2DQua21P;   break;

    //case 6004: return (ElementBase*) new GrowthElem2DTri22;    break;
    //case 6005: return (ElementBase*) new GrowthElem2DTri22F;   break;
    //case 6006: return (ElementBase*) new GrowthElem2DQua22;    break;
    //case 6007: return (ElementBase*) new GrowthElem2DTri221;   break;
    //case 6008: return (ElementBase*) new GrowthElem2DQua221;   break;

    //case 6051: return (ElementBase*) new GrowthElem3DHex1;       break;
    //case 6052: return (ElementBase*) new GrowthElem3DHex10;      break;
    case ELEM_GROWTH_3D_MIXED21: return (ElementBase*) new Elem_Growth_3D_Mixed21.h;    break;
    //case 6054: return (ElementBase*) new GrowthElem3DHex21;      break;
    //case 6055: return (ElementBase*) new GrowthElem3DTet21;      break;

    default: prgError(1,"NewElementGrowthFEM","unknown element type name!"); return NULL;
  }
}


#endif


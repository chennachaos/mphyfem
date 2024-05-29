
#ifndef  incl_newElementMagnetoFEM_h
#define  incl_newElementMagnetoFEM_h

//#include "ElementBase.h"
#include "FunctionsProgram.h"

#include "MagnetoMech2DTri21.h"                  // 2001
#include "MagnetoMech2DQua21.h"                  // 2002
#include "MagnetoMech2DTri221.h"                 // 2003
#include "MagnetoMech2DQua221.h"                 // 2004

#include "MagnetoMech2DQua10.h"                  // 2010

#include "MagnetoMech3DTet21.h"                  // 2051
#include "MagnetoMech3DWedge21.h"                // 2052
#include "MagnetoMech3DHex21.h"                  // 2053
#include "MagnetoMech3DTet221.h"                 // 2051
#include "MagnetoMech3DWedge221.h"               // 2052
#include "MagnetoMech3DHex221.h"                 // 2053

#include "MagnetoMech3DHex10.h"                  // 2060


# include "BernsteinElem3DFace.h"                 // 2082



inline  ElementBase*  NewElementMagnetomechFEM(int type)
{
  switch (type)
  {
    case 2001: return (ElementBase*) new MagnetoMech2DTri21;     break;
    case 2002: return (ElementBase*) new MagnetoMech2DQua21;     break;
    case 2003: return (ElementBase*) new MagnetoMech2DTri221;    break;
    case 2004: return (ElementBase*) new MagnetoMech2DQua221;    break;

    case 2010: return (ElementBase*) new MagnetoMech2DQua10;    break;

    case 2051: return (ElementBase*) new MagnetoMech3DTet21;     break;
    case 2052: return (ElementBase*) new MagnetoMech3DWedge21;   break;
    case 2053: return (ElementBase*) new MagnetoMech3DHex21;     break;
    case 2054: return (ElementBase*) new MagnetoMech3DTet221;    break;
    case 2055: return (ElementBase*) new MagnetoMech3DWedge221;  break;
    case 2056: return (ElementBase*) new MagnetoMech3DHex221;    break;

    case 2060: return (ElementBase*) new MagnetoMech3DHex10;     break;

    case 2083: return (ElementBase*) new BernsteinElem3DFace; break;

    default: prgError(1,"NewElementMagnetomechFEM","unknown element type name!"); return NULL;
  }
}


#endif


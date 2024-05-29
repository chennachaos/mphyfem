
#ifndef  incl_newElementMagnetoFEM_h
#define  incl_newElementMagnetoFEM_h

//#include "ElementBase.h"
#include "FunctionsProgram.h"

//#include "MagnetoMech2DTri21.h"                  // 2001
//#include "MagnetoMech2DQua21.h"                  // 2002
//#include "MagnetoMech2DTri221.h"                 // 2003
//#include "MagnetoMech2DQua221.h"                 // 2004
//#include "MagnetoMech2DQua10.h"                  // 2010

#include "Elem_Magnmech_3D_HM_10.h"
#include "Elem_Magnmech_3D_HM_21.h"
#include "Elem_Magnmech_3D_SM_221.h"


inline  ElementBase*  NewElementMagnetomechFEM(int type)
{
  switch (type)
  {
    //case 2001: return (ElementBase*) new MagnetoMech2DTri21;     break;
    //case 2002: return (ElementBase*) new MagnetoMech2DQua21;     break;
    //case 2003: return (ElementBase*) new MagnetoMech2DTri221;    break;
    //case 2004: return (ElementBase*) new MagnetoMech2DQua221;    break;

    //case 2010: return (ElementBase*) new MagnetoMech2DQua10;    break;

    case ELEM_MAGNMECH_3D_HM_10:  return (ElementBase*) new Elem_Magnmech_3D_HM_10;   break;
    case ELEM_MAGNMECH_3D_HM_21:  return (ElementBase*) new Elem_Magnmech_3D_HM_21;   break;
    case ELEM_MAGNMECH_3D_SM_221: return (ElementBase*) new Elem_Magnmech_3D_SM_221;  break;

    default: prgError(1,"NewElementMagnetomechFEM","unknown element type name!"); return NULL;
  }
}


#endif


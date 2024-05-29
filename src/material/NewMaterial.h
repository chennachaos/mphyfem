
#ifndef  incl_newMaterial_h
#define  incl_newMaterial_h

#include "MaterialBase.h"

#include "Matl_LinearElastic.h"
#include "Matl_SaintVenantKirchhoff.h"
#include "Matl_NeoHookean.h"
#include "Matl_MooneyRivlin.h"
#include "Matl_Gent.h"
#include "Matl_ArrudaBoyce.h"
#include "Matl_LinearElastic_Viscoelastic.h"
#include "Matl_NeoHookean_Viscoelastic.h"
#include "Matl_MooneyRivlin_Viscoelastic.h"
#include "Matl_Gent_Viscoelastic.h"
#include "Matl_Longevin8chain_Viscoelastic.h"


#include "Matl_ElecMech_Piezoelectric.h"
#include "Matl_ElecMech_NeoHookean.h"
#include "Matl_ElecMech_Gent.h"
#include "Matl_ElecMech_ArrudaBoyce.h"
#include "Matl_ElecMech_MooneyRivlin.h"

#include "Matl_MagnMech_NeoHookean.h"
#include "Matl_MagnMech_Gent.h"
#include "Matl_MagnMech_NeoHookean_Viscoelastic.h"
#include "Matl_MagnMech_Gent_Viscoelastic.h"
#include "Matl_MagnMech_Haldar.h"


#include "Matl_TherMech_Linear.h"

#include "Matl_Gel_Basic.h"



inline  MaterialBase*  NewMaterial(int type)
{
  cout << "type = " << type << endl;
  switch (type)
  {
    case    1: return (MaterialBase*) new Matl_LinearElastic; break;

    case    2: return (MaterialBase*) new Matl_SaintVenantKirchhoff; break;

    case    3: return (MaterialBase*) new Matl_NeoHookean; break;

    case    4: return (MaterialBase*) new Matl_MooneyRivlin; break;

    case    5: return (MaterialBase*) new Matl_Gent; break;

    case    6: return (MaterialBase*) new Matl_ArrudaBoyce; break;


    case    101: return (MaterialBase*) new Matl_LinearElastic_Viscoelastic; break;

    case    103: return (MaterialBase*) new Matl_NeoHookean_Viscoelastic; break;

    case    104: return (MaterialBase*) new Matl_MooneyRivlin_Viscoelastic; break;

    case    105: return (MaterialBase*) new Matl_Gent_Viscoelastic; break;

    case    106: return (MaterialBase*) new Matl_Longevin8chain_Viscoelastic; break;

    // electromechanical problem

    case  1001: return (MaterialBase*) new Matl_ElecMech_Piezoelectric; break;

    case  1002: return (MaterialBase*) new Matl_ElecMech_NeoHookean; break;

    case  1003: return (MaterialBase*) new Matl_ElecMech_Gent; break;

    case  1004: return (MaterialBase*) new Matl_ElecMech_ArrudaBoyce; break;

    case  1005: return (MaterialBase*) new Matl_ElecMech_MooneyRivlin; break;

    // magnetoomechanical problem

    case  2001: return (MaterialBase*) new Matl_MagnMech_NeoHookean; break;

    case  2002: return (MaterialBase*) new Matl_MagnMech_Gent; break;

    case  2005: return (MaterialBase*) new Matl_MagnMech_Haldar; break;

    case  2006: return (MaterialBase*) new Matl_MagnMech_NeoHookean_Viscoelastic; break;

    case  2007: return (MaterialBase*) new Matl_MagnMech_Gent_Viscoelastic; break;

    // thermomechanical problem

    case   50: return (MaterialBase*) new Matl_TherMech_Linear; break;

    // hydrogel problem

    case  5001: return (MaterialBase*) new Matl_Gel_Basic; break;


    default: cout << "NewMaterial ... unknown material type name!" << endl; return NULL;
  }
}


#endif


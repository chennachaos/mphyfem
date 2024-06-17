#ifndef  incl_myIDmaps_h
#define  incl_myIDmaps_h



inline  int  getDOFfromString(string& matkey)
{
    unordered_map<string,int>   map_keys_dofs = {
                        {"Xvelocity",       0},
                        {"Yvelocity",       1},
                        {"Zvelocity",       2},
                        //
                        {"phi",             0},
                        {"eta",             1},
                        //
                        {"VX",              0},
                        {"VY",              1},
                        {"VZ",              2},
                        //
                        {"Xdisplacement",   0},
                        {"Ydisplacement",   1},
                        {"Zdisplacement",   2},
                        //
                        {"UX",              0},
                        {"UY",              1},
                        {"UZ",              2},
                        };

    unordered_map<string,int>::const_iterator got = map_keys_dofs.find(matkey);

    if( got == map_keys_dofs.end() )
    {
        throw runtime_error("Key not found in getDOFfromString ...");
    }

    return got->second;
}


enum  {
     ELEM_SOLID_2D_DISP            = 1,
     ELEM_SOLID_2D_BBAR            = 2,
     ELEM_SOLID_2D_MIXED0          = 3,
     ELEM_SOLID_2D_MIXED1          = 4,
     ELEM_SOLID_2D_TRIA7_MIXED1    = 5,
     ELEM_SOLID_2D_TRIA10_MIXED1   = 6,
     //
     ELEM_SOLID_3D_DISP            = 51,
     ELEM_SOLID_3D_BBAR            = 52,
     ELEM_SOLID_3D_MIXED0          = 53,
     ELEM_SOLID_3D_MIXED1          = 54,
     ELEM_SOLID_3D_TETRA11_MIXED1  = 55,
     ELEM_SOLID_3D_TETRA20_MIXED1  = 56,
     //
     ELEM_TRUSS_2D_NODES2    = 101,
     ELEM_BEAM_2D_NODES2     = 102,
     ELEM_BEAM_2D_NODES3     = 103,
     ELEM_TRUSS_3D_NODES2    = 151,
     ELEM_BEAM_3D_NODES2     = 152,
     ELEM_BEAM_3D_NODES3     = 153,
     //
     ELEM_SHELL_FLAT_QUAD4   = 201,
     //
     //
     //
     ELEM_MAGNMECH_2D_HM_10  = 2001,
     ELEM_MAGNMECH_2D_HM_21  = 2002,
     ELEM_MAGNMECH_2D_SM_221 = 2003,
     //
     ELEM_MAGNMECH_3D_HM_10  = 2051,
     ELEM_MAGNMECH_3D_HM_21  = 2052,
     ELEM_MAGNMECH_3D_SM_221 = 2053,
     //
     //
     ELEM_GROWTH_3D_MIXED21  = 6053

};



inline  int  getElementID_Standard(std::string& matkey)
{
    unordered_map<string,int>   map_elements = {
     {"ELEM_SOLID_2D_DISP",              ELEM_SOLID_2D_DISP},
     {"ELEM_SOLID_2D_BBAR",              ELEM_SOLID_2D_BBAR},
     {"ELEM_SOLID_2D_MIXED0",            ELEM_SOLID_2D_MIXED0},
     {"ELEM_SOLID_2D_MIXED1",            ELEM_SOLID_2D_MIXED1},
     {"ELEM_SOLID_2D_TRIA7_MIXED1",      ELEM_SOLID_2D_TRIA7_MIXED1},
     {"ELEM_SOLID_2D_TRIA10_MIXED1",     ELEM_SOLID_2D_TRIA10_MIXED1},
     //
     {"ELEM_SOLID_3D_DISP",              ELEM_SOLID_3D_DISP},
     {"ELEM_SOLID_3D_BBAR",              ELEM_SOLID_3D_BBAR},
     {"ELEM_SOLID_3D_MIXED0",            ELEM_SOLID_3D_MIXED0},
     {"ELEM_SOLID_3D_MIXED1",            ELEM_SOLID_3D_MIXED1},
     {"ELEM_SOLID_3D_TETRA11_MIXED1",    ELEM_SOLID_3D_TETRA11_MIXED1},
     {"ELEM_SOLID_3D_TETRA20_MIXED1",    ELEM_SOLID_3D_TETRA20_MIXED1},
     //
     {"ELEM_TRUS_2D",                    0101},
     {"ELEM_BEAM_2D",                    0102},
     {"ELEM_TRUS_3D",                    0151},
     {"ELEM_BEAM_3D",                    0152},
     {"dummy",0},
     //
     {"ELEM_MAGNMECH_2D_HM_10",          ELEM_MAGNMECH_2D_HM_10},
     {"ELEM_MAGNMECH_2D_HM_21",          ELEM_MAGNMECH_2D_HM_21},
     {"ELEM_MAGNMECH_2D_SM_221",         ELEM_MAGNMECH_2D_SM_221},
     {"dumm1",                           2005},
     {"ELEM_MAGNMECH_3D_HM_10",          ELEM_MAGNMECH_3D_HM_10},
     {"ELEM_MAGNMECH_3D_HM_21",          ELEM_MAGNMECH_3D_HM_21},
     {"ELEM_MAGNMECH_3D_SM_221",         ELEM_MAGNMECH_3D_SM_221},
     {"dumm2",                           2055},
     //
     {"ELEM_GROWTH_3D_MIXED21",          6053}
    };

    unordered_map<string,int>::const_iterator got = map_elements.find(matkey);

    if( got == map_elements.end() )
    {
      throw runtime_error("Element key not found...");
      //return -1;
    }

    return got->second;
}


inline  int  getElementID_MagnetoMech(string& matkey)
{
    unordered_map<string,int>   map_elements = {
     {"ELEM_MAGNMECH_2D_HM_10",          ELEM_MAGNMECH_2D_HM_10},
     {"ELEM_MAGNMECH_2D_HM_21",          ELEM_MAGNMECH_2D_HM_21},
     {"ELEM_MAGNMECH_2D_SM_221",         ELEM_MAGNMECH_2D_SM_221},
     {"dumm1",                           2005},
     {"ELEM_MAGNMECH_3D_HM_10",          ELEM_MAGNMECH_3D_HM_10},
     {"ELEM_MAGNMECH_3D_HM_21",          ELEM_MAGNMECH_3D_HM_21},
     {"ELEM_MAGNMECH_3D_SM_221",         ELEM_MAGNMECH_3D_SM_221},
     {"dumm2",                           2055},
     {"dummy",0}
    };

    unordered_map<string,int>::const_iterator got = map_elements.find(matkey);

    if( got == map_elements.end() )
    {
      throw runtime_error("Element key not found...");
      //return -1;
    }

    return got->second;
}



inline  int getElementID_ElectroMech(std::string& matkey)
{
    unordered_map<string,int>   map_elements = {
     {"MagnetoMech2DTri21",              2001},
     {"MagnetoMech2DQua21",              2002},
     {"MagnetoMech2DTri221",             2003},
     {"MagnetoMech2DQua221",             2004},
     {"dumm1",                           2005},
     {"MagnetoMech3DTet21",              2051},
     {"MagnetoMech3DWedge21",            2052},
     {"MagnetoMech3DHex21",              2053},
     {"MagnetoMech3DTet221",             2054},
     {"MagnetoMech3DWedge221",           2055},
     {"MagnetoMech3DHex221",             2056},
     {"dumm2",                           2005},
     {"BernsteinElem3DFace",             2083},
     {"dummy",0}
    };

    unordered_map<string,int>::const_iterator got = map_elements.find(matkey);

    if( got == map_elements.end() )
    {
      throw runtime_error("Element key not found...");
      //return -1;
    }

    return got->second;
}


inline  int  getElementID_Growth(string& matkey)
{
    unordered_map<string,int>   map_elements = {
     {"GrowthElem2DQua1",          6001},
     {"GrowthElem2DQua10",         6002},
     {"GrowthElem2DQua21",         6003},
     {"GrowthElem2DTri21",         6004},
     {"GrowthElem2DQua21P",        6005},
     {"dummy",                     6006},
     {"dummy",                     6007},
     {"dummy",                     6008},
     {"GrowthElem3DHex1",          6051},
     {"GrowthElem3DHex10",         6052},
     {"GrowthElem3DWedge21",       6053},
     {"GrowthElem3DHex21",         6054},
     {"GrowthElem3DTet21",         6055},
     {"dummy",                     6056},
     {"dummy",                     6057},
     {"dummy",                     6058},
     {"dummy",                     6059},
     {"InverseGrowthElem3DHex1",   6551},
     {"InverseGrowthElem3DHex2",   6554}
    };

    unordered_map<string,int>::const_iterator got = map_elements.find(matkey);

    if( got == map_elements.end() )
    {
      throw runtime_error("Element key not found...");
      //return -1;
    }

    return got->second;
}





inline  int  getMaterialID(string& matkey)
{
    unordered_map<string,int>   map_materials = {
     {"Matl_LinearElastic",              1},
     {"Matl_SaintVenantKirchhoff",       2},
     {"Matl_NeoHookean",                 3},
     {"Matl_MooneyRivlin",               4},
     {"Matl_Gent",                       5},
     {"Matl_ArrudaBoyce",                6},
     {"dummy1",                          0},
     {"Matl_LinearElastic_Viscoelastic", 101},
     {"Matl_NeoHookean_Viscoelastic",    103},
     {"Matl_MooneyRivlin_Viscoelastic",  104},
     {"Matl_Gent_Viscoelastic",          105},
     {"Matl_Longevin8chain_Viscoelastic",106},
     {"dummy2",0},
     {"Matl_ElecMech_Piezoelectric",     1001},
     {"Matl_ElecMech_NeoHookean",        1002},
     {"Matl_ElecMech_Gent",              1003},
     {"Matl_ElecMech_ArrudaBoyce",       1004},
     {"Matl_ElecMech_MooneyRivlin",      1005},
     {"dummy3",0},
     {"Matl_MagnMech_NeoHookean",        2001},
     {"Matl_MagnMech_Gent",              2002},
     {"Matl_MagnMech_ArrudaBoyce",       2003},
     {"Matl_MagnMech_MooneyRivlin",      2004},
     {"Matl_MagnMech_Haldar",            2005},
     {"Matl_MagnMech_NeoHookean_Viscoelastic",        2006},
     {"Matl_MagnMech_Gent_Viscoelastic",              2007},
     {"dummy4",                            0},
     {"Matl_TherMech_Linear",            5001},
     {"Matl_TherMech_NeoHooke",          5002},
     {"Matl_TherMech_MooneyRivlin",      5003},
     {"dummy5",0},
     {"Matl_Gel_Basic",                  3001},
     {"dummy6",0},
     {"dummy",0}
    };

    unordered_map<string,int>::const_iterator got = map_materials.find(matkey);

    if( got == map_materials.end() )
    {
      throw runtime_error("Material key not found...");
      //return -1;
    }

    return got->second;
}


#endif

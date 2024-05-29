
#ifndef incl_headersBasic_h
#define incl_headersBasic_h


#include <iostream>
#include <limits.h>
#include <float.h>
#include <vector>
#include <assert.h>
#include <math.h>
#include <cmath>
#include <stdio.h>
#include <set>
#include <fstream>
#include <string.h>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include <memory>
#include <stdexcept>

#include <numeric>
#include <iterator>
#include <functional>

using std::vector;

typedef  vector<int>  vecIntSTL;
typedef  vector<float>  vecFltSTL;
typedef  vector<double>  vecDblSTL;


enum  PhysicsType  {PHYSICS_TYPE_FLUID=0, PHYSICS_TYPE_SOLID=1};

enum  SolverLibrary  {SOLVER_LIB_EIGEN=0, SOLVER_LIB_EIGEN_PARDISO=1, SOLVER_LIB_PETSC=2, SOLVER_LIB_PETSC_PARDISO=3};

enum  SolverType  {SOLVER_TYPE_NEWTONRAPHSON=0, SOLVER_TYPE_ARCLENGTH=1, SOLVER_TYPE_DEFLATION=2};

enum  FemModelBehaviour  {PSTRESS=1, PSTRAIN=2, AXISYM=3};

enum MeshConfiguration
{
  CONFIG_ORIGINAL=0, CONFIG_DEFORMED=1
};


enum ElementShape
{
  ELEM_SHAPE_TRIA    = 1,
  ELEM_SHAPE_QUAD    = 2,
  ELEM_SHAPE_TETRA   = 4,
  ELEM_SHAPE_PYRAMID = 5,
  ELEM_SHAPE_WEDGE   = 6,
  ELEM_SHAPE_HEXA    = 8,
  ELEM_TRIA_P2bP1dc  = 21,
  ELEM_SHAPE_TRIA_BERNSTEIN    = 101,
  ELEM_SHAPE_QUAD_BERNSTEIN    = 102,
  ELEM_SHAPE_TETRA_BERNSTEIN   = 103,
  ELEM_SHAPE_PYRAMID_BERNSTEIN = 104,
  ELEM_SHAPE_WEDGE_BERNSTEIN   = 105,
  ELEM_SHAPE_HEXA_BERNSTEIN    = 106
};


enum  class  TISFLUID {STEADY=0, BE=1, BDF1=1, BDF2=2, BDF3=3, BDF4=4, Midpoint=5, Galpha=6};


enum  class  TISSOLID {STEADY=0, STATIC=0, BE=1, BDF1=1, BDF2=2, Midpoint=5, Newmark=11, HHTalpha=12, CHalpha=13, KDPalpha=14, SEMIIMPL=400, CHalphaExplicit=500};




#endif



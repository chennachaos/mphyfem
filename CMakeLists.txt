cmake_minimum_required(VERSION 3.10)

#project(femMagnetomech)
project(femSolidmechanics)
#project(femGrowth)

message("<< ${CMAKE_PROJECT_NAME} / ${PROJECT_NAME}")

enable_language(CXX)

##################################
# project directories
##################################


#PROJECT_BINARY_DIR
#PROJECT_SOURCE_DIR

##################################
# compiler options
##################################

# C++ compiler

set(CMAKE_CXX_COMPILER "/usr/bin/g++")

#if(CXX_COMPILER_NAME MATCHES "gcc*")
message ("Using GNU C++ compiler")
#set(COPT "-std=c++11 -g -w -funroll-loops -ftree-vectorize -ftree-vectorizer-verbose=1  -Wno-write-strings -Wno-uninitialized -Wno-deprecated  -Wno-format -fpermissive -frounding-math -fstack-protector-all")

set(COPT "-std=c++17 -O3 -g -Wno-uninitialized -Wno-sign-compare -Wno-write-strings  -Wno-deprecated  -Wno-unused  -Wno-format -Wno-narrowing  -Wno-reorder  -Wreturn-type -Wno-ignored-attributes -fopenmp  -fpermissive -frounding-math -fstack-protector-all")
#endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${COPT}")


##################################
# include directories
##################################
include_directories(
src/lib
src/mySolvers
src/femSolidmechanics
#src/electromechFEM
src/femMagnetomech
#src/thermomechFEM
#src/gelFEM
src/femGrowth
#src/inversegrowthFEM/
src/utilities
src/material
/usr/include/
/usr/include/c++/13
/usr/include/eigen3
/usr/include/eigen3/unsupported
/usr/include/vtk-9.1
/usr/include/x86_64-linux-gnu/openmpi
/usr/include/petsc
/home/chenna/Documents/myCode/pardiso/panua-pardiso-20230908-linux/include
#/home/chenna/.spack/Spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/openmpi-4.1.4-sqt7euln7hmjk6lhssizrvf7qy7wfiq5/include
#/home/chenna/.spack/Spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/petsc-3.17.2-oxjyalypn3qx3gdx4koxmqsjye3bpcng/include
#/home/chenna/.spack/Spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/mumps-5.4.1-ajk5s6rsvgrbd3v5glzvibnanmcwpmz2/include
/usr/include/superlu
/usr/include/boost
)



##################################
# lib directories
##################################
link_directories(
/usr/lib/gcc/x86_64-linux-gnu/13
/usr/lib/x86_64-linux-gnu 
/usr/lib/blas
/usr/lib/lapack 
/usr/lib
#/home/chenna/Documents/myCode/mphyfem/lib
/home/chenna/Documents/myCode/pardiso/panua-pardiso-20230908-linux/lib
#/usr/lib/x86_64-linux-gnu/openmpi/lib 
#/usr/lib/petsc/lib
#/home/chenna/.spack/Spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/openmpi-4.1.4-sqt7euln7hmjk6lhssizrvf7qy7wfiq5/lib
#/home/chenna/.spack/Spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/petsc-3.17.2-oxjyalypn3qx3gdx4koxmqsjye3bpcng/lib
#/home/chenna/.spack/Spack/opt/spack/linux-ubuntu20.04-icelake/gcc-9.4.0/mumps-5.4.1-ajk5s6rsvgrbd3v5glzvibnanmcwpmz2/lib
)


##################################
# source files
##################################



file(GLOB sources1
src/program/*h
src/program/*cpp
src/electromechFEM/*h
src/electromechFEM/*cpp
src/material/*h
src/material/*cpp
src/mySolvers/*h
src/mySolvers/*cpp
src/utilities/*h
src/utilities/*cpp
src/lib/*h
src/lib/*cpp
src/macro/*h
src/macro/*cpp
)

file(GLOB sources1 src/*/*h)
file(GLOB sources2 src/*/*cpp)


file(GLOB sourcesGrowth1
src/femGrowth/*h
#src/inversegrowthFEM/*h
src/material/*h
src/mySolvers/*h
src/utilities/*h
)

file(GLOB sourcesGrowth2
src/femGrowth/*.cpp
#src/inversegrowthFEM/*.cpp
src/material/*.cpp
src/mySolvers/*.cpp
src/utilities/*.cpp
)


file(GLOB sourcesSolidmechanics1
src/femSolidmechanics/*h
src/material/*h
src/mySolvers/*h
src/utilities/*h
)


file(GLOB sourcesSolidmechanics2
src/femSolidmechanics/*.cpp
src/material/*.cpp
src/mySolvers/*.cpp
src/utilities/*.cpp
)



file(GLOB sourcesMagneto1
#src/growthFEM/*h
src/femMagnetomech/*h
src/material/*h
src/mySolvers/*h
src/utilities/*h
)

file(GLOB sourcesMagneto2
src/femMagnetomech/*.cpp
src/material/*.cpp
src/mySolvers/*.cpp
src/utilities/*.cpp
)


##################################
##################################

#if(${CMAKE_PROJECT_NAME} EQUAL "femSolidmechanics")
add_executable(${CMAKE_PROJECT_NAME}  src/femSolidmechanics.cpp ${sourcesSolidmechanics1} ${sourcesSolidmechanics2})
#elseif(${CMAKE_PROJECT_NAME} EQUAL "femMagnetomech")
#add_executable(${CMAKE_PROJECT_NAME}  src/femMagnetomech.cpp ${sourcesMagneto1} ${sourcesMagneto2})
#add_executable(${CMAKE_PROJECT_NAME}  src/femGrowth.cpp ${sourcesGrowth1} ${sourcesGrowth2})
#endif()

#pardiso600-GNU800-X86-64

target_link_libraries(${CMAKE_PROJECT_NAME} stdc++  gfortran gomp pthread m  blas  lapack  metis parmetis  pardiso  petsc  mpi  mpi_cxx  superlu  vtkCommonCore-9.1  vtkIOCore-9.1  vtkIOGeometry-9.1  vtkIOXML-9.1  vtkIOImage-9.1  vtkIOExport-9.1   vtkFiltersCore-9.1  vtkFiltersGeneral-9.1   vtkFiltersModeling-9.1  vtkFiltersSources-9.1  vtkCommonExecutionModel-9.1  vtkCommonDataModel-9.1)


install(TARGETS ${CMAKE_PROJECT_NAME} RUNTIME DESTINATION /home/chenna/Documents/myCode/mphyfem/bin)







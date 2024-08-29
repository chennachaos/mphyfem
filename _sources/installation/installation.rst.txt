
Installation and Execution
==========================


Dependencies
------------

* Programming languages: **C++**
* C++ standard: **C++14** (or above)
* Required third-party libraries are:

  1. CMake
  2. Blas
  3. Lapack
  4. Boost
  5. MPI (OpenMPI, MPICH or Intel MPI. Your choice!)
  6. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page)
  7. [METIS](http://glaros.dtc.umn.edu/gkhome/metis/metis/overview)
  8. [PETSc](https://www.mcs.anl.gov/petsc/)
  9. [VTK](https://vtk.org/)
  10. Metis and ParMetis


Steps for compiling and building
--------------------------------

1. Clone the repository or download the zip file and extract its contents.
2. Go to the directory of the repository in a terminal.
3. Create **build** and **bin** directories.
    * `mkdir build`
    * `mkdir bin`
4. Modify the CMake file accordingly.
    * Copy `CMakeLists-chennalaptop.txt` to a new file for your machine, say `CMakeLists-local.txt`.
    * Change the paths to the compilers, Eigen, CGAL, PETSc and VTK libraries.
    * Change the path in the `install` function.
    * Create a symbolic link to the local CMake file.
      * `ln -sf CMakeLists-local.txt CMakeLists.txt`
5. Enter the `build` directory.
    * `cd build`
6. Configure using the CMake file
    * `cmake ..`
7. Compile, build and install the executable `mpap`. This step will also copy the exe to the `bin` folder.
    * `make install`



Execution
---------
* Simulations are usually run from any folder with the appropriate input files.
* To do so, you need to add the `./mphyfem/bin` folder to the `PATH` environment variable.
* To avoid doing this every time you open a terminal, you can add it to `.bashrc` file. For example,

    * `export PATH=$PATH:/home/chenna/Documents/myCode/mphyfem/bin`

* Copy `petsc_options.dat` file from the `project/sampleinputs` folder to the `bin` folder.
* For executing the basic Solid Mechanics solver: `femSolidmechanics`
* For executing Growth solver: `femGrowth`
* For executing Magneto-Mechanics solver: `femMagnetomech`


Understanding and using output(s)
---------------------------------

For each simulation the output files are geneted in the VTU format which can be visualised using ParaView.


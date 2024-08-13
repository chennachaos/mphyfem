# MPHYFEM
Finite Element Method based simulation framework for multiphysics problems in solid and fluid mechanics.

A variety of finite element formulations are implemented.


* Programming languages: **C++**.
* C++ standard: **C++14** (or above)
* Required third-party libraries:
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
  11. SuperLU


## Compilation and building
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
8. Add the path to the project `bin` to the global `PATH` variable to detect the executable from any other location.
    * `export PATH=$PATH:/home/chenna/Documents/myCode/mphyfem/bin`
    * You can add the above command to `.bashrc` to avoid doing it everytime.

## Execution
* Simulations are usually run from any folder with a sub-folder named `inputs`.
* The project directory and file structure is as shown below.
    ```
    project
    |
    ----inputs
    |   |    config
    |   |    meshfile.msh
    |   |    petsc_options.dat
    ```
* To run the simulation, run the executable from the terminal.
* Example 1: To run using the default `config` file in the `inputs` sub-directory.
  * `./femSolidmechanics` 
* Example 2: To run using a different configuration file in the `inputs` sub-directory.
    * `./femSolidmechanics  configfile.config`


## Citing MPHYFEM
If you use MPHYFEM in your research, please cite it as:

```
@misc{mphyfem,
  author = {Chennakesava Kadapa},
  title = {{MPHYFEM}: An FEM-based simulation framework for multiphysics problems},
  year = {2024},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/chennachaos/mphyfem}},
}
```


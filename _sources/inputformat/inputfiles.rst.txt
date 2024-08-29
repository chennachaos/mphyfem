Simulation case
---------------
Each case requires three files stored in `inputs` directory.

1. Mesh file.
2. Configure file named `config`.
3. Petsc options file named `petsc_options.dat`.


Mesh file
---------

* The meshes are generated in the `GMSH` software and the input meshes are stored in the `.msh` format.
* All the domains (volumes in 3D, surfaces in 2D and curves/lines in 1D) should be tagged using `Physical groups` feature in GMSH.
* No need to tag all the boundaries though. Only those boundaries on which boundary conditions are applied need to be tagged.


Config file
-----------
The `config` file contains the details of the simulation such as element type, element properties, material properties, loading, solver, time step, etc.. A sample `config` file is shown below.

::

    Files
    {
      mesh : block3d-P2-nelem2
    }


    Model Type
    {
      dimension     : 3
    }


    Domains
    {
      block
      {
        Material   : 1
        Element    : 1
      }
    }


    Material
    {
      id               : 1
      name             : Matl_LinearElastic
      density          : 1.0
      data             : 240.565  0.3
    }


    Material
    {
      id               : 2
      name             : Matl_NeoHookean
      density          : 1.0
      data_deviatoric :  92.525
      data_volumetric :  3  0.0049883
    }

    Element
    {
      id           : 1
      type         : ELEM_SOLID_3D_MIXED1

    }


    Body Force
    {
        value        :  0  0  0
        TimeFunction :  1
    }


    Boundary Conditions
    {
        symX
        {
            type          : specified
            dof           : UX
            value         : 0
        }

        symY
        {
            type          : specified
            dof           : UY
            value         : 0
        }

        symZ
        {
            type          : specified
            dof           : UZ
            value         : 0
        }

        top
        {
            type          : specified
            dof           : UX
            value         : 0
        }

        top
        {
            type          : specified
            dof           : UY
            value         : 0
        }
    }


    Tractions
    {
        pressure
        {
            type          : traction
            value         : -320  0.0  0.0
            timefunction  : 1
        }
    }


    Time Functions
    {

    ! lam(t) = p1 + p2*t + p3*sin(p4*t+p5) + p6*cos(p7*t+p8)
    !
    ! id   t0      t1     p1   p2     p3    p4    p5    p6    p7    p8
       1   0.0   1000.0  0.0  1.0    0.0   0.0   0.0   0.0   0.0   0.0

    !   1   0.0      1.0  0.0  1.0    0.0   0.0   0.0   0.0   0.0   0.0
    !   1   1.0   1000.0  1.0  0.0    0.0   0.0   0.0   0.0   0.0   0.0

    }


    Solver
    {
        solvertype         newton

        !timescheme         BDF1
        timescheme         STEADY

        spectralRadius     0.0

        finalTime          1.0

        timeStep           0.2

        maximumSteps       10

        maximumIterations  20

        tolerance          1.0e-7

        debug              1

    }



As shown in the sample `config` file, the input is specified in the form of blocks.

Comments
---------
Any line that starts with a `!` in the config file is not effective; it is treated as a comment.
It helps the user to write useful comments or to try multiple options without having to remember the options.


Data blocks
------------
Different blocks are discribed below. Each block starts with `{` and ends with `}` on separate lines.


Files
******
Specifies the input mesh file and previous solution file for in case of `restart`. Different fields are

* **msh**: Prefix of the mesh file. *Compulsory*.
* **restart**: the name of the restart file with the solution. *Optional*.


Model Type
***********
Specifies the dimension and the model type for the 2D problem, e.g., plane stress, plane strain, axisymmetric.

* dimension: dimension of the model. Accepted values are 1, 2 or 3 depending on the problem.
* behaviour: TODO

Domains
********
Defines the element and material types for each of the (solid) domains. A problem can have more than one domain. This is especially the case for problems with multiple materials.

Each domain (volume in 3D, surface in 2D, line in 1D) should have a tag in the `.msh` file generated with GMSH. Material and element are assigned to each domain in the similar block format. For example, the following block assigns material with id 1 and element with id 1 to the domain named `block`. The material and elements are befined further below.


::

    block
    {
      Material     :  1
      Element      :  1
    }



Material
*********
Specifies the material type along with the relevant material property values such as Young's modulus, Poisson's ratio, density, etc..

* **id**: Identification number to keep track of list of materials.
* **name**: Material type name. Available material types are discribed here. TODO.
* **density**: Density of the material.
* **data**: Material property data
* **data_deviatoric**:
* **data_deviatoric**:


Element
********
Specifies the element type along with the relevant element data, if any.

* **id**: Identification number to keep track of list of elements.
* **type**: Element type name. Available element types are discribed here. TODO.


Body Force
**************
Specifies the body force.

* **value**: Value of the body force as a vector given as three values.
* **TimeFunction**: Time function id in case the body force is to be applied as a function of time. *Optional field*. Applied as the constant field if the `TimeFunction` is not specified.


Boundary Conditions
********************
Specifies the Dirichlet boundary conditions on the boundary patches.

The syntax is ::

    <patch-tag>
    {
        type          : <bc type>
        dof           : <DOF>
        value         : <value>
        TimeFunction  : <value>
    }


Available fields are:

* **type**: Boundary condition type.
    * `fixed`: to constrain the patch in all DOFs.
    * `specified`: to prescribe a value for a particular DOF.
* **dof**: Degree of freedom.
    * UX: X-displacement
    * UY: Y-displacement
    * UZ: Z-displacement
* **value**: The value to be applied.
* **timefunction**: To apply the boundary condition as a function of time. If not specified, boundary condition is applied as a constant value over time.


The following block specifies a fixed constraint on the patch tagged as `leftface`. ::

    leftface
    {
        type          : fixed
    }

The following block specifies a value of 0 (zero) for X-displacement (Dirichlet boundary condition) on the patch tagged as `symX`. ::

    symX
    {
        type          : specified
        dof           : UX
        value         : 0
    }


The following block specifies a value of 2.0 for Z-displacement on the patch tagged as `topface`, and it will be applied as a time-varying value as defined by the `time function with id 2`. ::

    topface
    {
        type          : specified
        dof           : UZ
        value         : 2.0
        timefunction  : 2
    }


Tractions
**********
Specifies traction boundary conditions on boundary patches. Required fields are
* **type**: 
* **value**: A vector of specified traction values in the normal and trangential directions. The first value is in the normal direction, and the other values are in the trangential directions.
* **timefunction**: The function id if the traction value is to be applied as a time-varying function.


::

    pressure
    {
        value         : 700.0  0.0  0.0
        timefunction  : 1
    }



Time functions
**************
This block contains the time functions for time-dependent boundary/loading or conditions or other time-dependent parameters. The time functions are specified as shown below by providing the coefficients p1, p2, p3, p4, p5, p6, p7 and p8.

::

    ! f(t) = p1 + p2*t + p3*sin(p4*t+p5) + p6*cos(p7*t+p8)
    !
    ! id   t0      t1     p1   p2     p3    p4    p5    p6    p7    p8
       1   0.0   1000.0  0.0  1.0    0.0   0.0   0.0   0.0   0.0   0.0

       2   0.0      1.0  0.0  1.0    0.0   0.0   0.0   0.0   0.0   0.0
       2   1.0   1000.0  1.0  0.0    0.0   0.0   0.0   0.0   0.0   0.0

       3   0.0   1000.0  0.0  0.0    2.0   6.2832   3.1416   0.0   0.0   0.0


The above block specifies three different time functions as given by the ids 1, 2 and 3.
* Time function 1: Constant value until 1000 time units.
* Time function 2: Linear variation with slope 1 until 1 time unit and then constant after 1 second until 1000 time units.
* Time function 3: Sinusoidal variation with an amplitude of two frequency f=1.0 Hz and phase angle of 3.1416 radians (180 degrees).


Solver
********
This block specifies what type of solver, nonlinear scheme, time step size, output frequency, max iterations, end time etc. Available fields are

* **solvertype**: newton (Newton-Raphson) or arclength (Arc-Length) method.
* **timescheme**: time integration scheme to be used.
    * STEADY: Quasi-static problem. No time integration.
    * BDF1: Backward-Euler scheme. This is a first-order scheme.
    * BDF2: Backward difference formula, which is of second-order accuracy.
    * CHalpha: CH-alpha scheme for second-order equations for dynamics. Second-order accurate and unconditionally stable.
    * KDPalpha: KDP-alpha scheme. Second-order accurate and unconditionally stable.
    * Newmark: Newmark-beta scheme. Second-order accurate and unconditionally stable.
* **spectralRadius**: Spectral radius value for CHalpha and KDPalpha schemes. Range is between 0 and 1, inclusve. A value of zero annihialites all high-frequency modes. A value of one add no numerical damping, similar to the Newmark-beta method.
* **finalTime**: The final time until which the simulation should be run.
* **timeStep**: The time step size.
* **maximumSteps**: Maximum number of time/load steps.
* **maximumIterations**: Maximum number of iterations per time/load step.
* **tolerance**: Convergence tolerance for iterations.
* **debug**: A flag used to print extra output for debugging purposes. 0 for deactivating debugging info and 1 for activating debugging info.


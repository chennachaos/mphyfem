
Example 1 - Cook's membrane in 2D with linear material
======================================================

In this example, we simulate the benchmark example of Cook's membrane in plane-strain condition with compressible linear elastic material using pure displacement formulation.

The configuration file is shown below.

::

    Files
    {
      mesh : cookmembrane2d-P2-nelem8
    }

    Model Type
    {
      dimension     : 2
      behaviour     : pstrain
      thickness     : 1
    }

    Domains
    {
      solid
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

    Element
    {
      id            : 1
      type          : ELEM_SOLID_2D_DISP
    }

    Body Force
    {
        value        :  0  0  0
        TimeFunction :  1
    }

    Boundary Conditions
    {
        leftedge
        {
            type   :  fixed
        }
    }

    Tractions
    {
        rightedge
        {
            type          : traction
            value         : 0.0  6.25  0.0
            timefunction  : 1
        }
    }

    Time Functions
    {

    ! f(t) = p1 + p2*t + p3*sin(p4*t+p5) + p6*cos(p7*t+p8)
    !
    ! id   t0    t1    p1   p2     p3    p4    p5    p6    p7    p8

       1   0.0   1000.0  0.0  1.0    0.0   0.0   0.0   0.0   0.0   0.0

    }

    Solver
    {

    solvertype       :  newton

    timescheme       :  STEADY

    spectralRadius   :  0.0

    finalTime        :  1.0

    timeStep         :  1.0

    maximumSteps     :  1

    maximumIterations :  2

    tolerance         :  1.0e-7

    debug             :  0

    }


The contour plot of displacement magnitude, along with the element edges, is shown in the figure below.

.. image:: cooksmembrane2d-P2-nelem8-LE-nu0p3-dispMagn.png
  :width: 400



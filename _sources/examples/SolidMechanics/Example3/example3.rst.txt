Example 3 - Cook's membrane in 2D with incompressible Neo-Hookean material
===========================================================================

In this example, we simulate the benchmark example of Cook's membrane in plane-strain with incompressible Neo-Hookean material.

The configuration file is quite similar to the ones used in Examples 1 and 2. Only the important blocks of the config file are shown below.

## Material block

::

    Material
    {
      id               : 1
      name             : Matl_NeoHookean
      density          : 1.0
      data_deviatoric :  80.188
      data_volumetric :  -1  0.0049883
    }

Similar to Example 2, we specifiy the value of shear modulus, which is 80.188, for the deviatoric part.

For the volumetric part, the first entry, `-1`, specifies that the model is truly incompressible. In this case, the second entry is ignored.

We apply the load in several steps by using a ramp profile, as done in Example 2. Since we use truly incompressible material, we must use the mixed displacement-pressure element. So, we change the element type in `Element` block.

::

    Element
    {
      id            : 1
      type          : ELEM_SOLID_2D_MIXED1
    }

`ELEM_SOLID_2D_MIXED1` is the Taylor-Hood element. For the P2 element, quadratic element is used for the displacement field and the linear continuous element for the pressure field.


The contour plot of displacement magnitude, along with the element edges, is shown in the figure below. The displacement decreases further. This is due to the fact that the material becomes stiffer as the Poisson's ratio approaches the incompressibility limit.

.. image:: cooksmembrane2d-P2-nelem8-NH-nu0p5-dispMagn.png
  :width: 400





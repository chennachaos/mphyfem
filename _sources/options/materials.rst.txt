
Material types
==============


Solid Mechanics
----------------

The list of available material types is provided below along with a brief descrption and the options for each material.

* **Matl_LinearElastic**
    * Linear elastic material
    * Young's modulus and Poisson's ratio are specified with the field `data` as ::

        data: 100.0  0.3

    where 100.0 is the Young's modulus and 0.3 is the Poisson's ratio.

* **Matl_SaintVenantKirchhoff**
    * Saint Venant-Kirchhoff material, which is an extension of the linear elastic material to finite strains.

* **Matl_NeoHookean**
    * The strain energy function is assumed to be split into deviatoric and volumetric parts.
    * Input data is specified as the parameters for deviatoric part of the strain energy function with `data_deviatoric` field and for the volumetric part with `data_volumetric` field.
    
    ::

        data_deviatoric: <shear-modulus>
        data_volumetric: <U-type>  <inverse-of-bulk-modulus>

    For example, ::

        data_deviatoric: 100.0
        data_volumetric: 3 0.000001
    
    where 100.0 is the shear modulus, 3 is the type of the volumetric function and 0.000001 is the inverse of bulk modulus.

    We can specify truly incompressible material by specifying `Utype` as -1. That is, ::

        data_deviatoric: 100.0
        data_volumetric: -1  0.000001
    
    In this case, the inverse of bulk modulus value is ignored as it is not relevant for truly incompressible materials.


* **Matl_MooneyRivlin**
    * Similar to Neo-Hookean material, input is specified with `data_deviatoric` and `data_volumetric` fields.
    * `data_deviatoric` takes the coefficients C10 and C01.

* **Matl_Gent**
    * Similar to Neo-Hookean material, input is specified with `data_deviatoric` and `data_volumetric` fields.
    * `data_deviatoric` takes the value of Im.

* **Matl_ArrudaBoyce**
    * Similar to Neo-Hookean material, input is specified with `data_deviatoric` and `data_volumetric` fields.
    * `data_deviatoric` takes the coefficients.

.. * **Matl_LinearElastic_Viscoelastic**

* **Matl_NeoHookean_Viscoelastic**
    * Similar to Neo-Hookean material, input is specified with `data_deviatoric` and `data_volumetric` fields.
    * Parameters of the viscoelastic part is specified with `data_viscoelastic` as::

        data_viscoelastic :  <N>  <mu1>  <tau1>  <mu2>  <tau2>  ... <muN>  <tauN>
    
    For example, ::

        data_viscoelastic :  1  5e5  0.01
    
    specifies one data point with shear modulus of 5e5 and relaxation time of 0.01 seconds.


.. * **Matl_MooneyRivlin_Viscoelastic**
.. * **Matl_Gent_Viscoelastic**
.. * **Matl_Longevin8chain_Viscoelastic**


Magneto Mechanics
------------------
* Hard magnetic (HM) materials: the magnetic field is specified as the residual and applied magnetic fields in the input.
* Soft magnetic (SM) materials: the magnetic field is solved as the solution along with the displacements.


Growth (Morphoelasticity)
--------------------------


Electro Mechanics
-----------------


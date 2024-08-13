
Element types
==============


Solid Mechanics
----------------

The list of available element types is provided below along with a brief descrption and the options for each element.

2D elements
^^^^^^^^^^^
* **ELEM_SOLID_2D_DISP**
    * 2D element with pure displacement formulation. Full integration is used.
* **ELEM_SOLID_2D_BBAR**
    * 2D element with Bbar formulation in small strains and Fbar formulation in finite strains.
* **ELEM_SOLID_2D_MIXED0**
    * Q1/P0 element.
    * 2D element with mixed displacement-pressure formulation with element-wise constant pressure.
    * Only available for 4-noded quadrilateral and 6-noded triangular elements.
* **ELEM_SOLID_2D_MIXED1**
    * 2D element with mixed displacement-pressure formulation with element-wise linear and continuous pressure.
    * Only available for 6-noded triangular and 9-noded quadrilateral elements.
* **ELEM_SOLID_2D_TRIA7_MIXED1**
    * 2D element with mixed displacement-pressure formulation.
    * 7-noded triangular element (P2b, P2 with a cubic bubble) for displacements and element-wise linear and discontinuous pressure (P1dc).
* **ELEM_SOLID_2D_TRIA10_MIXED1**
    * 2D element with mixed displacement-pressure formulation.
    * 10-noded triangular element (P3, cubic element) for displacements and element-wise linear and discontinuous pressure (P1dc).

3D elements
^^^^^^^^^^^
* **ELEM_SOLID_3D_DISP**
    * 3D element with pure displacement formulation. Full integration is used.
* **ELEM_SOLID_3D_BBAR**
    * 3D element with Bbar formulation in small strains and Fbar formulation in finite strains.
* **ELEM_SOLID_3D_MIXED0**
    * Q1/P0 element.
    * 3D element with mixed displacement-pressure formulation with element-wise constant pressure.
    * Only available for 8-noded hexahedral and 10-noded tetrahedral elements.
* **ELEM_SOLID_3D_MIXED1**
    * 3D element with mixed displacement-pressure formulation with element-wise linear and continuous pressure.
    * Only available for 10-noded tetrahedral (P2/P1), 18-noded wedge (W2/W1) and 27-noded hexahedral (Q2/Q1) elements.
* **ELEM_SOLID_3D_TETRA11_MIXED1**
    * 3D element with mixed displacement-pressure formulation.
    * 11-noded tetrahedral element (P2b, P2 with a cubic bubble) for displacements and element-wise linear and discontinuous pressure (P1dc).
* **ELEM_SOLID_3D_TETRA20_MIXED1**
    * 3D element with mixed displacement-pressure formulation.
    * 20-noded tetrahedral element (P3, cubic element) for displacements and element-wise linear and discontinuous pressure (P1dc).


Magneto Mechanics
------------------
* Hard magnetic (HM) materials: the magnetic field is specified as the residual and applied magnetic fields in the input.
* Soft magnetic (SM) materials: the magnetic field is solved as the solution along with the displacements.

2D elements
^^^^^^^^^^^
* **ELEM_MAGNMECH_2D_HM_10**
    * 2D Q1/P0 element for HM materials.
    * 4-noded quadrilateral for displacement and element-wise constant pressure.
* **ELEM_MAGNMECH_2D_HM_21**
    * 2D P2/P1 or Q2/Q1 element for HM materials.
    * Quadratic element for displacement and linear element for pressure.
* **ELEM_MAGNMECH_2D_SM_221**
    * 2D P2/P1/P2 or Q2/Q1/Q2 element for SM materials.
    * Quadratic element for displacement, linear element for pressure, and quadratic element for magnetic field.

3D elements
^^^^^^^^^^^
* **ELEM_MAGNMECH_3D_HM_10**
    * 3D Q1/P0 element for HM materials.
    * 8-noded quadrilateral for displacement and element-wise constant pressure.
* **ELEM_MAGNMECH_3D_HM_21**
    * 3D P2/P1 or Q2/Q1 element for HM materials.
    * Quadratic element for displacement and linear element for pressure.
* **ELEM_MAGNMECH_3D_SM_221**
    * 3D P2/P1/P2 (Tetrahedron) or W2/W1/W2 (Wedge) or Q2/Q1/Q2 (Hexahedron) element for SM materials.
    * Quadratic element for displacement, linear element for pressure, and quadratic element for magnetic field.


Growth (Morphoelasticity)
--------------------------

2D elements
^^^^^^^^^^^
* **ELEM_GROWTH_2D_MIXED21**
    * 2D P2/P1 or Q2/Q1 element for Morphoelasticity.
    * Quadratic element for displacement and linear element for pressure.

3D elements
^^^^^^^^^^^
* **ELEM_GROWTH_3D_MIXED21**
    * 3D P2/P1 (Tetrahedron) or W2/W1 (Wedge) or Q2/Q1 (Hexahedron) element for Morphoelasticity.
    * Quadratic element for displacement and linear element for pressure.


Electro Mechanics
-----------------


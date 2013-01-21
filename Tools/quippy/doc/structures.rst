:mod:`quippy.structures` Structure generation routines
=========================================================

.. currentmodule: quippy.structures

.. automodule:: quippy.structures
    :synopsis: Tools for generating crystals and other atomic structures

.. function:: diamond(a[, z])

   Return an 8-atom diamond-structure with cubic lattice constant `a` and
   atomic number `z`, e.g.::

      a = diamond(5.44, 14)  # Silicon unit cell


.. function:: bcc(a[, z])

   Return a 2-atom bcc structure with cubic lattice constant `a` and
   atomic number `z`.


.. function:: fcc(a[, z])

   Return a 4-atom fcc-structure with cubic lattice constant `a` and
   atomic number `z`.
 

.. function:: a15(a[, z])

   Return an 8-atom a15-structure with cubic lattice constant `a` and
   atomic number `z`.


.. function:: alpha_quartz(a, c, u, x, y, z)

   Return 9-atom primitve trigonal alpha-quartz unit cell, with lattice constants
   `a` and `c` and internal coordinates `u` (Si), `x`, `y` and `z` (O).


.. function:: alpha_quartz_cubic(a, c, u, x, y, z)

   Return non-primitve 18-atom cubic alpha-quartz unit cell.


.. function:: graphene_cubic(a)

   Return cubic graphene unit cell with lattice parameter `a`.

.. function:: graphene_sheet(a, n, m, rep_x, repy_y)

   Return a graphene sheet with index `(n,m)` and lattice constant
   `a`, with `rep_x` repeat units in the x-direction and `rep_y`
   in the y-direction.

.. function:: graphene_tube(a, n, m, nz)

   Construct a `(n,m)` nanotube with lattice parameter `a` and `nz`
   unit cells along the tube length. Return a tuple `(tube, radius)`
   where `tube` is the new :class:`Atoms` object and `radius` is the
   radius of the nanotube. Example usage::

      >>> tube, r = graphene_tube(1.42, 6, 6, 1) # (6,6) nanotube, 1 unit cell long in z direction

.. function:: tube_radius(at)

   Return average radius of the nanotube `at`.


.. function:: slab(axes, a, [nx, ny,  width, height, nz, atnum, lat_type])

   Return a slab of material with the x, y, and z axes desribed by the
   Miller indices in the array axes (with ``x = axes[:,1])``, ``y =
   axes[:,2]`` and ``z = axes[:,3]``).  The extent of the slab should
   be given either as ``(nx, ny, nz)`` unit cells or as ``(width,
   height, nz)`` where `width` and `height` are measured in Angstrom
   and `nz` is the number of cells in the `z` direction.
   
   `atnum` can be used to initialise the `z` and `species` properties.
   `lat_type` should be of ``"diamond"```, ``"fcc"``, or ``"bcc"``
   (default is ``"diamond"``)

.. function:: supercell(at, n1, n2, n3)

   Replicates the unit cell `at` `n1` x `n2` x `n3` times along its
   lattice vectors.

.. function:: transform(at, t)

   Transform cell and lattice coordinates by the 3 x 3 matrix `t` and
   return a new :class:`Atoms` object.


.. function:: water()

   Return an :class:`Atoms` object containing one TIP3P water molecule in a box
   giving the correct density at 300K




.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
.. HQ X
.. HQ X   quippy: Python interface to QUIP atomistic simulation library
.. HQ X
.. HQ X   Copyright James Kermode 2010
.. HQ X
.. HQ X   These portions of the source code are released under the GNU General
.. HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
.. HQ X
.. HQ X   If you would like to license the source code under different terms,
.. HQ X   please contact James Kermode, james.kermode@gmail.com
.. HQ X
.. HQ X   When using this software, please cite the following reference:
.. HQ X
.. HQ X   http://www.jrkermode.co.uk/quippy
.. HQ X
.. HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

:class:`quippy.Atoms` -- store and manipulate atoms
===================================================

.. currentmodule:: quippy

.. class:: Atoms([source, n, lattice, data, properties, params, *readargs, **readkwargs])

   An :class:`Atoms` object is used to store and manipulate a
   collection of atoms which together form a molecule or crystal.

   :class:`Atoms` objects store two main kinds of data:
   :attr:`properties` and :attr:`params`. Properties are used
   for data which are extensive with the size of the system, e.g.
   position, velocity, acceleration or local energy. Parameters
   are used for intensive variables e.g. temperature or total
   energy.

   When constructing an :class:`Atoms` object, if the `source`
   argument is present then :meth:`Atoms.read` is invoked to
   initialise the :class:`Atoms` object from `source` which would
   typically be a filename, for example::

      >>> at = Atoms('input.xyz')

   Otherwise, the arguments `n`, `lattice`, `data`,
   `properties` and `params` are used to initialise the
   object. `n` (default value 0) sets the number of atoms and
   `lattice` (default identity matrix) the lattice vectors. To create
   an :class:`Atoms` object to store 10 atoms with a cubic unit cell
   of dimensions 10 x 10 x 10 A, you would write::

       >>> at = Atoms(n=10, lattice=10.0*fidentity(3))

   (The :func:`fidentity` function creates a :class:`FortranArray`
   representing an identity matrix). The optional arguments `data`,
   `properties` and `parameters` can be used to initialise
   attributes of the object. If one of `data` or `properties` is
   present then so must the other.

   Important attributes and methods are listed below. In the example
   code `at` refers to an :class:`Atoms` instance.

   .. attribute:: n
   
	Number of atoms. Also available as ``len(at)``.

   .. attribute:: lattice

	3 x 3 :class:`FortranArray` matrix with lattice vectors as columns:

	.. math::
	    \left(
	    \begin{array}{ccc}
	     | & | & | \\
	     \mathbf{a} & \mathbf{b} & \mathbf{c} \\
	     | & | & | \\
	    \end{array}
	    \right)
	     = \left(
	    \begin{array}{ccc}
	     R_{11} & R_{12} & R_{13} \\
	     R_{21} & R_{22} & R_{23} \\
	     R_{31} & R_{32} & R_{33} \\
	    \end{array}
	    \right)

   .. attribute:: g
    
      3 x 3 :class:`FortranArray` matrix of reciprocal lattice vectors, 
      i.e. the inverse of the :attr:`lattice` matrix.


   .. attribute:: cutoff

      Should be set using :meth:`set_cutoff` or :meth:`set_cutoff_factor`.
      If :attr:`use_uniform_cutoff` is true, `cutoff` is the cutoff
      distance in Angstrom.  Otherwise, it's a multiplier for
      ``bond_length(Zi,Zj)``.


   .. attribute:: cutoff_break

      If connectivity is being calculated hysteretically, a bond will
      be considered to have broken when it exceeds `cutoff_break`.

   .. attribute:: nneightol

      Two atoms count as nearest neighbours if they are closer than
      the sum of their covalent radii multiplied by
      `nneightol`. Default is 1.2.
      

   .. attribute:: params

      :class:`Dictionary` of parameters. Useful for storing data about
      this :class:`Atoms` object, for example the temperature, total
      energy or applied strain. The data stored here is automatically
      saved to and loaded from XYZ and NetCDF files.

   .. attribute:: properties

      :class:`Dictionary` of atomic properties. A property is an array
      of shape (`m`,`n`) where `n` is the number of atoms and `m` is
      either one (for scalar properties) or three (vector
      properties). Properties can be integer, real, string or logical.
      String properties have a fixed length of 10 characters.

      Each property is automatically visible as a
      :class:`FortranArray` attribute of the :class:`Atoms` object,
      for example the atomic positions are stored in a real vector
      property called `pos`, and can be accessed as ``at.pos``.

      Properties can be added with the :meth:`add_property` method and
      removed with :meth:`remove_property`.

   .. attribute:: connect

      Connectivity information, stored in a :class:`Connection`
      object. The connectivity can be calculated with the
      :meth:`calc_connect` method, and is best acceseed via the
      :meth:`n_neighbours` and :meth:`neighbour` methods.

   .. attribute:: hysteretic_connect

      Secondary :class:`Connection` object used for storing hysteretic
      bonding information. Hysteretic connectivity can be calculated
      with :meth:`calc_connect_hysteretic` and accessed by calling
      :meth:`neighbour` with ``alt_connect=at.hysteretic_connect``.

   .. attribute:: data

      :class:`Table` object storing the data for all of the atomic
      properties.  No need to access directly; use the automatically
      created view arrays such as `pos`, `z`, `velo`, etc.


   .. attribute:: neighbours

      A :class:`Neighbours` object, which, when indexed with an
      integer from 1 to `at.n`, returns an array of
      :class:`NeighbourInfo` objects, each of which corresponds to a
      particular pair `(i,j)` and has attributes `j`, `distance`,
      `diff`, `cosines` and `shift`.

      If connectivity information has not already been calculated
      :meth:`calc_connect` will be called automatically. The code to
      loop over the neighbours of all atoms is quite idiomatic::
      
        for i in frange(at.n):
	    for neighb in at.neighbours[i]:
		print (neighb.j, neighb.distance, neighb.diff, neighb.cosines, neighb.shift)

      Note that this attribute provides a more Pythonic interface to
      the atomic connectivity information than the wrapped Fortran
      functions :meth:`n_neighbours` and :meth:`neighbour` described
      below.

   .. attribute:: hysteretic_neighbours

      The same as :attr:`neighbours`, but this one accesses the
      :attr:`hysteretic_connect` :class:`Connection` object rather than
      :attr:`connect`. :attr:`cutoff` and :attr:`cutoff_break` should be
      set appropriately and :meth:`calc_connect_hysteretic` called
      before accessing this attribute.

   .. method:: read(source[, format, *args, **kwargs])

      Class method to read an :class:`Atoms` object from `source`. If
      `format` is ``None`` then it is inferred automatically, either
      from the file extension if `source` is a filename, or from
      the class of `source`. 

      If `source` corresponds to a known format then it used
      to construct an appropriate iterator from the :attr:`AtomsReaders`
      dictionary. See :ref:`fileformats` for a list of supported
      file formats. 

      If `source` corresponds to an unknown format then it is
      expected to be an iterator returning :class:`Atoms` objects.


   .. method:: write(dest[, format, properties, *args, **kwargs])

      Write this :class:`Atoms` object to `dest`. If `format` is
      absent it is inferred from the file extension or type of `dest`,
      as described for the :meth:`read` method.  If `properties` is
      present, it should be a list of property names to include in the
      output file, e.g. `['species', 'pos']`.

      See :ref:`fileformats` for a list of supported file formats.

   .. method:: show([property, arrows])

      Show this :class:`Atoms` object in AtomEye using the quippy
      :mod:`~quippy.atomeye` module.  If `property` is present it
      should be the name of a scalar property
      (e.g. ``"local_energy"``) or a rank one array of length ``at.n``
      to be used to colour the atoms. If `arrows` is present it should
      be the name of a vector property (e.g. ``"force"``) or a rank
      two array with shape (3, ``at.n``) to be used to draw arrows
      originating from each atom.


   .. method:: select([mask, list])

      Return an :class:`Atoms` containing a subset of the atoms in
      this object. One of either `mask` or `list` should be
      present. If `mask` is given it should be a rank one array of
      length `self.n`. In this case atoms corresponding to true
      values in `mask` will be included in the result.  If `list`
      is present it should be an arry of list containing atom indices
      to include in the result.


   .. method:: copy()

      Return a copy of this :class:`Atoms` object.
      

   .. method:: add_atoms(pos, z[, mass, velo, acc, travel, data])

      Add one ore more atoms to an :class:`Atoms` object. To add a
      single atom, `pos` should be an array of size 3 and `z` a
      single integer. To add multiple atoms either arrays of length
      `n_new` should be passed, or `data` should be given as a
      :class:`Table` object with all new atomic data in the same
      format as the existing :attr:`data` table

      
   .. method:: has_property(name)

      Return true if property `name` exists. Property names
      are case-insensitive.


   .. method:: add_property(name, value[, n_cols])

      Add a new property to this Atoms object.

      `name` is the name of the new property and `value` should be
      either a scalar or an array representing the value, which should
      be either integer, real, logical or string.

      If a scalar is given for `value` it is copied to every element
      in the new property.  `n_cols` can be specified to create a 2D
      property from a scalar initial value - the default is 1 which
      creates a 1D property.

      If an array is given for `value` it should either have shape
      (self.n,) for a 1D property or (n_cols,self.n) for a 2D
      property.  In this case `n_cols` is inferred from the shape of
      the `value` and shouldn't be passed as an argument.

      If property with the same type is already present then no error
      occurs. A warning is printed if the verbosity level is VERBOSE
      or higher. The value will be overwritten with that given in
      `value`.

      Here are some examples::

	  a = Atoms(n=10, lattice=10.0*fidentity(3))

	  a.add_property('mark', 1)                  # Scalar integer
	  a.add_property('bool', False)              # Scalar logical
	  a.add_property('local_energy', 0.0)        # Scalar real
	  a.add_property('force', 0.0, n_cols=3)     # Vector real
	  a.add_property('label', '')                # Scalar string

	  a.add_property('count', [1,2,3,4,5,6,7,8,9,10])  # From list
	  a.add_property('norm_pos', a.pos.norm())         # From 1D array
	  a.add_property('pos', new_pos)                   # Overwrite positions with array new_pos
	                                                   # which should have shape (3,10)
      
   .. method:: remove_property(name)

      Remove the property `name`.

   .. method:: remove_atoms(remove_list)

      Remove one or more atoms. `remove_list` can be either a single
      integer or a list of atom indices to be removed.


   .. method:: distance(i, j, shift)

      Return distance between atoms with indices `i` and `j` if
      they are separated by `shift` periodic cells:

      .. math::

         r_{ij} = \left| \mathbf{r}_j - \mathbf{r}_i + \mathbf{R} \cdot  \mathbf{s} \right|

      where :math:`\mathbf{R}` is the :attr:`lattice` matrix and :math:`\mathbf{s}` the shift.


   .. method:: distance_min_image([i, j, v, w, shift])

      Return minimum image distance between two atoms or positions.
      End points can be specified by any combination of atoms indices
      `i` and `j` and absolute coordinates `u` and `w`.

      If `shift` is present it should be an farray of dtype `int32`, which
      on exit will contain the periodic shift between the two atoms or points.

     
   .. method:: diff(i, j, shift)

      Return difference vector between atoms `i` and `j` if they are separated
      by `shift` periodic cells:

      .. math::

         \mathbf{u}_{ij} = \mathbf{r}_j - \mathbf{r}_i + \mathbf{R} \cdot  \mathbf{s}

      where :math:`\mathbf{R}` is the :attr:`lattice` matrix and :math:`\mathbf{s}` the shift.


   .. method:: diff_min_image([i, j, v, w])

      Return the minimum image difference vector between two atoms or
      positions. End points can be specified by any combination of
      atoms indices `i` and `j` and absolute coordinates `u` and
      `w`.

   .. method:: cosine(i, j, k)

      Return cosine of the angle j--i--k


   .. method:: cosine_neighbour(i, n, m)

      Return cosine of the angle n--i--m where `n` and `m` are 
      the nth and  mth neighbours of atom `i`


   .. method:: direction_cosines(i, j, shift)

      Given two atoms `i` and `j` and a `shift`, return the direction
      cosines of the difference vector from `i` to `j`.


   .. method:: direction_cosines(i, j)

      Direction cosines of of the minimum image difference
      vector from `i` to `j`.
      

   .. method:: set_atoms(z[, mass])

      Set atomic numbers (in the `z` integer property), species names
      (in `species` string property) and optionally masses (if `mass`
      property exists in :class:`Atoms` object).
      

   .. method:: map_into_cell

      Map atomic positions into the unit cell so that lattice
      coordinates satisfy :math:`-0.5 \le t_x,t_y,t_z < 0.5`


   .. method:: calc_connect([own_neighbour])

      Fast :math:`O(N)` connectivity calculation routine. It divides the
      unit cell into similarly shaped subcells, of sufficient size
      that sphere of radius :attr:`cutoff` is contained in a subcell, at
      least in the directions in which the unit cell is big enough. For very
      small unit cells, there is only one subcell, so the routine is
      equivalent to the standard :math:`O(N^2)` method.

      If `own_neighbour` is true, atoms can be neighbours with their
      own periodic images.


   .. method:: calc_connect_hysteretic([alt_connect, origin, extent, own_neighbour])

      As for :meth:`calc_connect`, but perform the connectivity update
      hystertically: atoms must come within :attr:`cutoff` to be considered
      neighbours, and then will remain connect until them move apart
      further than :attr:`cutoff_break`.

      Typically `alt_connect` should be set to the
      :attr:`hysteretic_connect` attribute. `origin` and `extent`
      vectors can be used to restrict the hysteretic region to only
      part of the entire system -- the :func:`estimate_origin_extent`
      function can be used to guess suitable values.


   .. method:: calc_dists()

      Updates the stored distance tables using the stored connectivity
      and shifts. This should be called every time any atoms are moved
      (e.g. it is called by :meth:`DynamicalSystem.advance_verlet`).


   .. method:: n_neighbours(i[, max_dist, max_factor, alt_connect])

      Return the number of neighbours that atom `i` has. If
      `max_dist` or `max_factor` is specified, then only neighbour closer
      than this cutoff are included. `alt_connect` can be set to another
      :class:`Connection` object to use alternative connectivity information, for
      example :attr:`hysteretic_connect`.


   .. method:: neighbour(i, n[, distance, diff, cosines, shift, index_bn, max_dist, jn, alt_connect])

      Return the index of the nth neighbour of atom `i`. Together
      with :meth:`n_neighbours`, this facilites a loop over the
      neighbours of atom `i`. Optionally, we return other geometric
      information, such as distance, direction cosines and difference
      vector, and also a direct index into the neighbour tables. If
      :math:`i \le j`, this is an index into
      ``connect.neighbour1[i]``, otherwise it is an index into
      ``connect.neighbour1[j]``. 

      Here's a typical loop construct. Note how `r` and `u` are
      created before the loop: arguments which are both optional and
      ``intent(out)`` in Fortran are converted to ``intent(in,out)``
      for quippy.::

	r = farray(0.0)
	u = fzeros(3)    
	for i in frange(at.n):
	    for n in frange(at.n_neighbours(i)):
		j = at.neighbour(i, n, distance=r, diff=u)
    
      If ``distance > max_dist``, we return 0, and do not waste time
      calculating other quantities. `alt_connect` has the same meaning as
      in :meth:`n_neighbours`.

      
   .. method:: is_nearest_neighbour(i, n)

      Test if atom i's nth neighbour is one of its nearest neighbours.

   
   .. method:: cell_volume():

      Return the (unsigned) volume of the simulation cell.

   
   .. method:: centre_of_mass([index_list, origin])

      Calculate the centre of mass of an atoms object, using the
      closest images to the `origin` atom, or first atom if this is
      not specified.  If an `index_list` is present, just calculate it
      for that subset of atoms (then the `origin` atom is the first in
      this list unless it is specified separately).
    
      .. note:: 

	 Because the origin can be specified separately it need not
	 be one of the atoms in `index_list`.

    
   .. method:: realpos(i)

      Return real position of atom `i`, taking into account the
      stored travel across the periodic boundary conditions.

   .. method:: set_cutoff(cutoff[, cutoff_break])

      Specify a uniform cutoff throughout the system. Optionally set
      `cutoff_break` for :meth:`calc_connect_hysteretic` at the same time.


   .. method:: set_cutoff_factor(factor[, factor_break])

      Specify the neighbour cutoff to be a multiple of the bond length
      of the two atoms' types. Optional argument `factor_break` is
      used by :meth:`calc_connect_hysteretic`.

   .. method:: density()

      Return density in units of :math:`g/m^3`. If a `mass` property
      exists, use that, otherwise we use `z` and `ElementMass` to
      calculate the total mass of the cell.

      
Structure generation routines
-----------------------------

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



Cluster carving routines
------------------------

.. function:: create_hybrid_weights(at, args_str)

   Given an atoms structure with a `hybrid_mark` property, this
   routine creates a `weight_region1` property, whose values are
   between 0 and 1.

.. function:: create_cluster_info_from_mark(at, args_str[, cut_bonds])

   Create a cluster using the 'hybrid_mark' property and options in
   `args_str`.  All atoms that are marked with anything other than
   ``HYBRID_NO_MARK`` will be included in the cluster; this includes
   active, transition and buffer atoms.

   Optionally return a list of the bonds cut when making the cluster in the 
   :class:`Table` `cut_bonds`.


.. function:: carve_cluster(at, args_str, cluster_info)

   Carve a cluster from `at` using the information in the
   :class:`Table` `cluster_info`, which should be generated by a
   call to :func:`create_cluster_info_from_hybrid_mark`.

   The output cluster contains all properties of the initial atoms object, and
   some additional columns, which are:

    ``index``
	index of the cluster atoms into the initial atoms object.

    ``termindex`` 
       nonzero for termination atoms, and is an index into
       the cluster atoms specifiying which atom is being terminated,
       it is used in collecting the forces.

    ``rescale`` 
	a real number which for nontermination atoms is 1.0,
        for termination atoms records the scaling applied to
        termination bond lengths
    
    ``shift``
	the shift of each atom

.. function:: construct_buffer(at, core, radius[, use_avgpos, verbosity])
   
   Given an atoms object, and a list of core atoms in the first
   integer of the `core` table, fill the `buffer` table with all atoms
   within `radius` of any core atom (which can be reached by
   connectivity hopping).  Optionally use the time averaged positions.


.. function:: create_embed_and_fit_list(at, fit_hops,[nneighb_only,min_images_only])

   Given an Atoms structure with an active region marked in the
   `hybrid_mark` property using ``HYBRID_ACTIVE_MARK``, grow the embed
   region by `fit_hops` bond hops to form a fit region. 

   Returns the  embedlist and fitlist as :class:`Table` objects
   with correct periodic shifts.


.. function:: estimate_origin_extent(at, active, cluster_radius)

   Return estimated `(origin, extent)` for hysteretic connectivity 
   calculator to include all atoms within a distance `cluster_radius`
   of an atom in the `active` array.


Miscellaneous routines
----------------------
   
.. function:: elastic_fields(at, a, c11, c12, c44)

   Calculate atom resolved local strain and stress for all
   fourfold-coordinated atoms in `at`. Crystal structure should be
   diamond.

   Properties are added to `at`. `a` should be the relaxed lattice
   constant in A, and `c11`, `c12` and `c44` are the independent
   elastic constants for the diamond structure, in GPa.

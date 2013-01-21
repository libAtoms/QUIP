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

:mod:`quippy.atoms` -- Representation of atomic configurations
==============================================================

.. automodule:: quippy.atoms
   :synopsis: Representation of an atomic configuration

.. autoclass:: Atoms

   Important :class:`Atoms` attributes include:

   .. attribute:: n
   
       Number of atoms. Also available as ``len(at)``.

   .. attribute:: lattice

       3 x 3 matrix with lattice vectors as columns:

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
    
      3 x 3 matrix of reciprocal lattice vectors,
      i.e. the inverse of the :attr:`lattice` matrix.

  
   .. attribute:: pos

      3 x `N` array of atomic positions, stored inside :attr:`properties` dictionary.


   .. attribute:: Z

      Array of atomic numbers, stored inside :attr:`properties` dictionary.


   .. attribute:: species

      Array of atomic symbols, stored inside :attr:`properties` dictionary.

   
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
		print (neighb.j, neighb.distance, neighb.diff,
		       neighb.cosines, neighb.shift)

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

   .. autoattribute:: info

   .. autoattribute:: arrays

   .. autoattribute:: indices


   Important :class:`Atoms` methods include:

   .. automethod:: iteratoms

   .. automethod:: equivalent

   .. automethod:: read

   .. automethod:: read_from

   .. automethod:: copy

   .. automethod:: copy_from

   .. automethod:: write

   .. automethod:: select

   .. automethod:: md5_hash

   .. automethod:: get_atom

   .. automethod:: print_atom

   .. automethod:: density

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

   .. automethod:: add_property

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

   .. automethod:: mem_estimate


   Important attributes and methods inherited
   from :class:`ase.atoms.Atoms` include the following. See the ASE
   documentation of the :class:`ase.atoms.Atoms` for full details.

   .. autoattribute:: positions

   .. autoattribute:: numbers

   .. autoattribute:: cell

   .. autoattribute:: constraints

   .. automethod:: get_positions

   .. automethod:: set_positions

   .. automethod:: get_atomic_numbers

   .. automethod:: set_atomic_numbers

   .. automethod:: get_cell

   .. automethod:: set_cell

   .. automethod:: get_momenta

   .. automethod:: set_momenta

   .. automethod:: get_calculator

   .. automethod:: set_calculator

   .. automethod:: get_potential_energy

   .. automethod:: get_potential_energies

   .. automethod:: get_forces

   .. automethod:: get_stress

   .. automethod:: get_stresses


.. autoclass:: Neighbours
   :members:

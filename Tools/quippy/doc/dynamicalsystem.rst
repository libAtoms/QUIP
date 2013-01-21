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

:mod:`quippy.dynamicalsystem` -- run molecular dynamics simulation
====================================================================

.. currentmodule:: quippy.dynamicalsystem

.. class:: DynamicalSystem(atoms[, velocity, acceleration, constraints, rigidbodies])
   
   A DynamicalSystem object stores the dynamical state of a
   molecular dynamics simulation. It has methods to integrate the
   equations of motion using the Velocity-Verlet algorithm.

   A DynamicalSystem is constructed from an :class:`Atoms` object
   `atoms`. The initial velocities and accelerations can optionally be
   specificed as `(3,atoms.n)` arrays. The `constraints` and
   `rigidbodies` arguments can be used to specify the number of
   constraints and rigid bodies respectively (the default is zero in
   both cases).

   Attributes:

   .. attribute:: atoms

      :class:`Atoms` object used to store all dynamical variables. Important :attr:`~Atoms.properties`
      include:
	
	`pos`
		Real vector. Positions of atoms.
	`velo`
		Real vector. Velocity of each atom.
	`acc`
		Real vector. Acceleration of each atom.
	`avgpos`
		Real vector. Time-averaged position of each atom.
	`oldpos`
		Real vector. Position of each atom at previous time step.
	`mass`
		Real scalar. Atomic masses.
	`travel`
		Integer vector. Travel across the periodic boundary conditions.
	`move_mask`
		Integer scalar. Mask used to fix atoms. Defaults to 1 for all atoms.
	`damp_mask`
		Integer salar. Mask to select atoms for damping. Defaults to 1 for all atoms.
	`thermostat_mask`
		Integer scalar. Mask to select which thermostat to use for each
		atom. Defaults to 1 for all atoms.
	`avg_ke`
		Real scalar. Time-averaged kinetic energy of each atom.

   .. autoattribute:: avg_temp
   .. autoattribute:: avg_time
   .. autoattribute:: cur_temp
   .. autoattribute:: dw
   .. autoattribute:: epot
   .. autoattribute:: ext_energy
   .. autoattribute:: initialised
   .. autoattribute:: n
   .. autoattribute:: nconstraints
   .. autoattribute:: ndof
   .. autoattribute:: nrigid
   .. autoattribute:: nsteps
   .. autoattribute:: random_seed
   .. autoattribute:: t
   .. autoattribute:: thermostat_dw
   .. autoattribute:: thermostat_work
   .. autoattribute:: work

   Important DynamicalSystem methods are listed below:

   .. method:: add_atoms(z[mass, p, v, a, t, data])

      Equivalent of :meth:`Atoms.add_atoms`, but also amends the number of degrees
      of freedom correctly.

   .. method:: remove_atoms(remove_list)
 
      Equivalent of :meth:`Atoms.remove_atoms`, but also amends the number of degrees
      of freedom correctly.
     
   .. method:: add_heat(heat, ekin)

      Add or remove heat from the system by scaling the atomic
      velocities.

   .. method:: add_thermostat(type_bn, t[, gamma, q, tau, p])

      Add a new thermostat to this DynamicalSystem. `type_bn` should
      be one of 0 (no thermostat), 1 (Langevin), 2 (Nose-Hoover), 3
      (Nose-Hoover-Langevin) or 4 (Langevin NPT). `t` is the target
      temperature. `q` is the Nose-Hoover coupling constant. Only one
      of `tau` or `gamma` should be given. `p` is the external
      pressure for the case of Langevin NPT.


   .. method:: advance_verlet1(dt[, virial])

      Advance the velocities by half the time-step `dt` and the
      positions by a full time-step. A typical MD loop should
      resemble the following code. See the :ref:`moleculardynamics` 
      section of the quippy tutorial for more details. ::

         ds.atoms.calc_connect() 
	 for n in range(n_steps):
	    ds.advance_verlet1(dt)
	    pot.calc(ds.atoms, force=True, energy=True)
	    ds.advance_verlet2(dt, ds.atoms.force)
	    ds.print_status(epot=ds.atoms.energy)
	    if n % connect_interval == 0:
	       ds.atoms.calc_connect()
      

   .. method:: advance_verlet2(dt, f[, virial])

      Advances the velocities by the second half time-step.

   .. method:: advance_verlet(dt, f[, virial])

      Equivalent to :meth:`advance_verlet2` followed by
      :meth:`advance_verlet1`. This allows a simple MD loop::

	  for n in range(n_steps):
	      pot.calc(ds.atoms, force=True, energy=True)
	      ds.advance_verlet(dt, ds.atoms.force)
	      ds.print_status(epot=ds.atoms.energy)
	      if n % connect_interval == 0:
	         ds.atoms.calc_connect()

   .. method:: angular_momentum([origin, indices])

      Return the angular momentum of all the atoms in this DynamicalSystem, defined by

      .. math::
          
	  \mathbf{L} = \sum_{i} \mathbf{r_i} \times \mathbf{v_i}

      Optionally the `origin` can be specified, or only specific atoms in the list `indices`
      can be included.


   .. method:: centre_of_mass_velo([index_list])

      Return the velocity of the center of mass, optionally including
      only the atoms in `index_list`.

   .. method:: centre_of_mass_acc([index_list])

      Return the acceleration of the center of mass, optionally including
      only the atoms in `index_list`.
      

   .. method:: enable_damping(damp_time)

      Enable damping, with damping time set to `damp_time`. Only atoms
      flagged in the `damp_mask` property will be affected.

   .. method:: disable_damping()

      Disable damping for this DynamicalSystem.

   .. method:: kinetic_energy()

      Return the total kinetic energy :math:`E_k = \sum_{i} \frac{1}{2} m v^2`

   .. method:: momentum([indices])

      Return the total momentum :math:`\mathbf{p} = \sum_i \mathbf{m_i} \mathbf{v_i}`.
      Optionally only include the contribution of a subset of atoms in the array `indices`.

   .. method:: print_status([label, epot, instantaneous])

      Print a status line showing the current time, temperature, the
      mean temperature the total energy and the total momentum for this
      DynamicalSystem. If present, the optional `label` parameter should
      be a one character label for the log lines and is printed in the
      first column of the output. `epot` should be the potential energy
      if this is available. `instantaneous` has the same meaning as
      in :meth:`temperature`.

   .. method:: rescale_velo(temp[, mass_weighted, zero_l])

      Rescale the atomic velocities to temperature `temp`. If the
      current temperature is zero, we first randomise the
      velocities. If `mass_weighted` is true, then the velocites are
      weighted by :math:`1/sqrt{m}`. Linear momentum is zeroed
      automatically.  If `zero_l` is true then the angular momentum is
      also zeroed.

   .. method:: run(pot, [dt, n_steps, save_interval, connect_interval, args_str])

      Generator to return snapshots from a trajectory. For each step,
      forces are evaluated using the :class:`Potential` `pot` and
      the DynamicalSystem is advanced by a time `dt` (default 1 fs).
      `n_steps` (default 10 steps) are carried out in total, with
      the generator yielding a result every `save_interval`
      steps. The connectivity is recalculated every
      `connect_interval` steps.  `args_str` can be used to supply
      extra arguments to :meth:`Potential.calc`.

   .. method:: temperature(this[, region, include_all, instantaneous])

      Return the temperature, assuming each degree of freedom
      contributes :math:`\frac{1}{2}kT`. By default only moving and
      thermostatted atoms are included --- this can be overriden by
      setting `include_all` to true. `region` can be used to restrict
      the calculation to a particular thermostat
      region. `instantaneous` controls whether the calculation should
      be carried out using the current values of the velocities and
      masses, or whether to return the value at the last Verlet step
      (the latter is the default).

   .. method:: zero_momentum([indices])

      Change velocities to those that the system would have in the
      zero momentum frame.  Optionalally zero the total momentum of a
      subset of atoms, specified by `indices`.


.. autoclass:: Dynamics
   :members:


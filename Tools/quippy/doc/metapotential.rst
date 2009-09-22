:class:`quippy.MetaPotential`
=============================

.. currentmodule:: quippy

.. class:: MetaPotential(args_str, pot[, pot2, bulk_scale])
   
   A MetaPotential is a generalised :class:`Potential`, where the
   requirement for the forces to be the derivative of a single
   Hamiltonian is relaxed - this is needed for QM/MM force-mixing,
   amongst other things.

   It is constructed with an `args_str` which specifies the
   type of MetaPotential and one or more Potential objects.
   `Simple` MetaPotentials contains only a single Potential::

      pot = Potential('IP SW', xml_string)
      mp = MetaPotential(simple, pot)

   An example of a more complex MetaPotential is a `ForceMixing`
   MetaPotential, which combines forces from two different
   potentials, commonly one QM and one classical::

      mm_pot = Potential('IP SW', sw_xml_string)
      qm_pot = Potential('TB DFTB', dftb_xml_string)
      metapot = MetaPotential('ForceMixing', mm_pot, qm_pot)

   The `bulk_scale` argument can be used if the two Potentials being
   combined predict different bulk properties. Space can be rescaled
   to match the lattice constants, and energy to match the bulk
   moduli. If you want ot do this, `bulk_scale` should be a suitable
   Atoms object to calibrate these rescalings.

   The third type of MetaPotential is a `Sum` MetaPotential. This
   combines two Potentials simply by adding (or subtracting) the
   contributions.

   Attributes:

   .. attribute:: is_forcemixing
      
      True if this is a `ForceMixing` MetaPotential.

   .. attribute:: is_sum
      
      True if this is a `Sum` MetaPotential.

   .. attribute:: is_simple
      
      True if this is a `Simpl` MetaPotential.


   Methods:

   .. method:: cutoff()

      As :meth:`Potential.cutoff`, but combine the cutoffs of the
      component Potentials appropriately.


   .. method:: calc(at[, e, local_e, f, df, virial, args_str, err]))

      As :meth:`Potential.calc`. 


   .. method:: minim(at,method,convergence_tol,max_steps,[linminroutine,do_print,print_inoutput,print_cinoutput,do_pos,do_lat,args_str,eps_guess,use_n_MetaPotential.minim,use_fire,lattice_fix,hook_print_interval,status])

      :param at: Atoms object
      :param method: Minimisation method - one of ``"cg"`` or ``"sd"``
      :param convergence_tol: Value of :math:`|\mathbf{f}|^2` for convergence
      :param max_steps: Maximum number of minimisation steps
      :param linminroutine: Line minimisation routine. Should be one of ``"NR_LINMIN"``, ``"FAST_LINMIN"`` or ``"LINMIN_DERIV"``
      :param do_print: if true, print configurations using minim's hook()
      :param print_inoutput: InOutput object to print configs to, needed if `do_print` is true
      :param print_cinoutput: CInOutput object to print configs to, needed if `do_print` is true
      :param do_pos: do relaxation w.r.t. positions and/or lattice (is neither is included, do both)
      :param do_lat: do relaxation w.r.t. positions and/or lattice (is neither is included, do both)
      :param args_str: extra arguments to pass to :meth:`calc`
      :param eps_guess: `eps_guess` argument to pass to minim
      :param use_n_minim: if true, use :func:`n_minim` instead of :func:`minim`
      :param use_fire: if true, use :func:`fire_minim` instead of minim
      :param lattice_fix: 3x3 array mask to fix some components of lattice. Defaults to all false.
      :param hook_print_interval: how often to print from hook function
      :param status: set to 1 if an error occurred during minimisation

      Minimise the configuration `at` under the action of this
      MetaPotential.  Returns number of minimisation steps taken. If
      an error occurs or convergence is not reached within `max_steps`
      steps, `status` will be set to 1 on exit.

      Example usage (see :ref:`geomopt` section of tutorial for full explanation)::
      
         at0 = diamond(5.44, 14)
	 at0.calc_connect()
	 pot = Potential('IP SW', """<SW_params n_types="1">
	         <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984)</comment>
		 <per_type_data type="1" atomic_num="14" />

		 <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
		   p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

		 <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
		   lambda="21.0" gamma="1.20" eps="2.1675" />
		</SW_params>""")
         metapot = MetaPotential('Simple', pot)
	 metapot.minim(at0, 'cg', 1e-7, 100, do_pos=True, do_lat=True)
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

:class:`quippy.Potential` -- evaluate interatomic potentials
============================================================

.. currentmodule:: quippy

.. class:: Potential(args_str,[pot1,pot2,param_str,bulk_scale,mpi_obj,error])

   A :class:`Potential` object represents an interatomic potential, a
   tight binding model or an interface to an external code used to
   perform calculations. It is initialised from an `args_str` describing
   the type of potential, and an XML formatted string `param_str`
   giving the parameters.

   Types of Potential:

     ====================  ==========================================================
     `arg_str` prefix      Description
     ====================  ==========================================================   
     ``IP``                Interatomic Potential
     ``TB``                Tight Binding Model
     ``FilePot``           File potential, used to communicate with external program
     ``CallbackPot``       Callback potential, computation done by Python function
     ====================  ==========================================================

   
   Types of interatomic potential available:

     ======================= ==========================================================
     ``arg_str`` prefix      Description
     ======================= ==========================================================   
     ``IP ASAP``             ASAP quartz potential
     ``IP BOP``              Bond order potential for metails
     ``IP Brenenr``          Brenner (1990) potential for carbon
     ``IP Brenner_2002``     Brenner (2002) reactive potential for carbon
     ``IP Brenner_Screened`` Pastewka modified Brenner reactive potential for carbon
     ``IP EAM_ErcolAd``      Embedded atom potential of Ercolessi and Adams
     ``IP FB``               Flikkema and Bromley potential for SiO2
     ``IP FS``	  	     Finnis=Sinclair potential for metals
     ``IP GAP``              Gaussian approximation potential
     ``IP LJ``               Lennard-Jones potential
     ``IP Si_MEAM``          Silicon modified embedded attom potential
     ``IP SW``               Stillinger-Weber potential for silicon
     ``IP Tersoff``          Tersoff potential for silicon
     ======================= ==========================================================   

   Types of tight binding potential available:

     ======================= ==========================================================
     ``arg_str`` prefix      Description
     ======================= ==========================================================   
     ``TB Bowler``           Bowler tight binding model 
     ``TB DFTB``             Density functional tight binding
     ``TB GSP``              Goodwin-Skinner-Pettifor tight binding model
     ``TB NRL_TB``           Naval Research Laboratory tight binding model
     ======================= ==========================================================   

   
   Examples of the XML parameters for each of these potential can be
   found in the `QUIP_Core/parameters
   <http://src.tcm.phy.cam.ac.uk/viewvc/jrk33/repo/trunk/QUIP/QUIP_Core/parameters>`_
   directory of the QUIP svn repository.

   Important methods of Potential objects are listed below.

   .. method:: cutoff

      Return the cutoff of this Potential, in Angstrom. This is the
      minimum neighbour connectivity cutoff that should be used: if
      you're doing MD you'll want to use a slightly larger cutoff so
      that new neighbours don't drift in between connectivity updates.

   .. method:: calc(at[, energy, local_energy, force, virial, local_virial, args_str])

      Apply this Potential to the Atoms object `at`, which must have
      connectivity information (i.e. :meth:`Atoms.calc_connect` should
      have been called). The optional arguments determine what should
      be calculated and how it will be returned. Each physical
      quantity has a corresponding optional argument, which can either
      be an `True` to store the result inside the Atoms object
      (i.e. in :attr:`Atoms.params` or in :attr:`Atoms.properties`) with
      the default name, a string to specify a different property or
      parameter name, or an array of the the correct shape to receive
      the quantity in question, as set out in the table below.

	================  ============= =============== ================================
	Array argument    Quantity      Shape           Default	storage location
	================  ============= =============== ================================
	`energy`          Energy        ``()``  	`energy` param
	`local_energy`    Local energy  ``(at.n,)``     `local_energy` property
	`force`           Force         ``(3,at.n)``    `force` property
	`virial`          Virial tensor ``(3,3)``       `virial` param
	`local_virial`    Local virial   ``(3,3,at.n)`` `local_virial` property
	================  ============= =============== ================================

      For the energy, a rank-0 array is required to store the results
      of this calculation. This is necessary since the underlying
      Fortran argument is both ``optional`` and ``intent(out)``. We
      can avoid this inconvenience by passing `energy=True` argument
      which returns the result in the `energy` parameter of the Atoms
      object `at`.

      The `args_str` argument can be a string or dictionary containing
      additional arguments which depend on the particular Potential
      being used. As mentioned in the documentation of
      :meth:`oo_fortran.FortranDerivedType._runroutine`, as a special
      case unexpected keyword arguments will be converted to
      additional `args_str` options.

      Not all Potentials support all of these quantities: a
      :exc:`RuntimeError` will be raised if you ask for something that
      is not supported.


   .. method:: minim(at,method,convergence_tol,max_steps,[linminroutine,do_print,print_inoutput,print_cinoutput,do_pos,do_lat,args_str,eps_guess,use_n_Potential.minim,use_fire,lattice_fix,hook_print_interval,status])

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
      Potential.  Returns number of minimisation steps taken. If
      an error occurs or convergence is not reached within `max_steps`
      steps, `status` will be set to 1 on exit.

      Example usage (see :ref:`geomopt` section of tutorial for full explanation)::
      
         at0 = diamond(5.44, 14)
	 at0.calc_connect()
	 pot = Potential('IP SW', param_str="""<SW_params n_types="1">
	         <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984)</comment>
		 <per_type_data type="1" atomic_num="14" />

		 <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
		   p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

		 <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
		   lambda="21.0" gamma="1.20" eps="2.1675" />
		</SW_params>""")
	 pot.minim(at0, 'cg', 1e-7, 100, do_pos=True, do_lat=True)

   .. method:: set_callback(at, callback)

      For a :class:`Potential` of type `CallbackPot`, this method is
      used to set the callback function. `callback` should be a Python
      function (or other callable, such as a bound method or class
      instance) which takes a single argument, of type
      :class:`Atoms`. Information about which quantities should be
      computed can be obtained from the `calc_energy`, `calc_local_e`,
      `calc_force`, and `calc_virial` keys in `at.params`. Results
      should be returned either as `at.params` entries (for energy and
      virial) or by adding new atomic properties (for forces and local
      energy).

      Here's an example implementation of a simple callback::

        def example_callback(at):
	    if at.calc_energy:
	       at.params['energy'] = ...
	      
	    if at.calc_force:
	       at.add_property('force', 0.0, n_cols=3)
	       at.force[:,:] = ...

	p = Potential('CallbackPot')
	p.set_callback(example_callback)
	p.calc(at, energy=True)
	print at.energy
	...
 

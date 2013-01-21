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

:mod:`quippy.potential` -- evaluate interatomic potentials
==========================================================

.. currentmodule:: quippy.potential

.. autoclass:: Potential
   :members:

   .. method:: calc(at, energy, force, virial, local_energy, local_virial, args_str, error)
    
        Apply this :class:`Potential` to the :class:`~quippy.atoms.Atoms` object
        `at`, which must have connectivity information
        (i.e. :meth:`Atoms.calc_connect` should have been called). The
        optional arguments determine what should be calculated and how
        it will be returned. Each physical quantity has a
        corresponding optional argument, which can either be an `True`
        to store the result inside the Atoms object (i.e. in
        :attr:`Atoms.params` or in :attr:`Atoms.properties`) with the
        default name, a string to specify a different property or
        parameter name, or an array of the the correct shape to
        receive the quantity in question, as set out in the table
        below.

          ================  ============= =============== ================================
          Array argument    Quantity      Shape           Default storage location
          ================  ============= =============== ================================
          `energy`          Energy        ``()``  	  `energy` param
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


   .. method:: minim(at, method, convergence_tol, max_steps, linminroutine[, do_print, print_inpoutput, print_cinoutput, do_pos, do_lat, args_str, eps_guess, fire_minim_dt0, fire_minim_dt_max, external_pressure, use_precond, hook_print_interval, error])

        :param at: Atoms object
        :param method: Minimisation method - one of ``"cg"``,  ``"cg_n"``,  ``"pcg"``, ``"sd"``,  ``"lbfgs"``,  ``"fire"``
        :param convergence_tol: Value of :math:`|\mathbf{f}|^2` for convergence
        :param max_steps: Maximum number of minimisation steps
        :param linminroutine: Line minimisation routine. Should be one of ``"NR_LINMIN"``, ``"FAST_LINMIN"`` or ``"LINMIN_DERIV"``
        :param do_print: if true, print configurations using minim's hook()
        :param print_inoutput: :class:`InOutput` object to print configs to, needed if `do_print` is true
        :param print_cinoutput: :class:`CInOutput` object to print configs to, needed if `do_print` is true
        :param do_pos: do relaxation w.r.t. positions and/or lattice (is neither is included, do both)
        :param do_lat: do relaxation w.r.t. positions and/or lattice (is neither is included, do both)
        :param args_str: extra arguments to pass to :meth:`calc`
        :param eps_guess: `eps_guess` argument to pass to minim
        :param fire_minim_dt0: Initial time step for FIRE minimiser
        :param fire_minim_dt_max: Maximum time step for FIRE minimiser
        :param external_pressure: 3x3 array giving external pressure tensor
        :param use_precond: If true, use preconditioner for minimisation (only implemented for ``method="cg_n"``)
        :param hook_print_interval: how often to print from hook function
        :param error: set to 1 if an error occurs during minimisation


        Minimise the configuration `at` under the action of this
        :class:`Potential`.  Returns number of minimisation steps taken. If
        an error occurs or convergence is not reached within `max_steps`
        steps, `status` will be set to 1 on exit.

        Example usage (see :ref:`geomopt` section of the quippy tutorial for full explanation)::

             at0 = diamond(5.44, 14)
             at0.calc_connect()
             pot = Potential('IP SW', param_str='''<SW_params n_types="1">
                     <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984)</comment>
                     <per_type_data type="1" atomic_num="14" />

                     <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
                       p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

                     <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
                       lambda="21.0" gamma="1.20" eps="2.1675" />
                    </SW_params>''')
             pot.minim(at0, 'cg', 1e-7, 100, do_pos=True, do_lat=True)
      

    .. method:: cutoff()

        Return the cutoff of this :class:`Potential`, in Angstrom. This is the
        minimum neighbour connectivity cutoff that should be used: if
        you're doing MD you'll want to use a slightly larger cutoff so
        that new neighbours don't drift in between connectivity updates
        (see the :attr:`cutoff_skin` attribute)

 
.. autoclass:: Minim
   :members:

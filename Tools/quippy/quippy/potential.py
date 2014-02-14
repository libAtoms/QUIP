# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import logging
import weakref
from math import sqrt

import numpy as np

from quippy import _potential
from quippy._potential import *
from quippy import get_fortran_indexing
from quippy.clusters import HYBRID_NO_MARK, HYBRID_ACTIVE_MARK
from quippy.oo_fortran import update_doc_string
from quippy.atoms import Atoms
from quippy.util import quip_xml_parameters, dict_to_args_str
from quippy.elasticity import stress_matrix
from quippy.farray import fzeros

__doc__ = _potential.__doc__
__all__ = _potential.__all__ + ['force_test', 'Minim', 'ForceMixingPotential']

potlog = logging.getLogger('quippy.potential')

def calculator_callback_factory(calculator):
    """Return a Python function which can be used as a quippy
       :class:`Potential` callback routine for the given ASE calculator."""

    def callback(at):
        at.set_calculator(calculator)
        if at.calc_energy:
            at.params['energy'] = at.get_potential_energy()
        if at.calc_local_e:
            at.add_property('local_e', 0.0, overwrite=True)
            at.local_e[:] = at.get_potential_energies()
        if at.calc_force:
            at.add_property('force', 0.0, n_cols=3, overwrite=True)
            at.force[:] = at.get_forces().T
        if at.calc_virial:
            stress = at.get_stress()
            virial = fzeros((3,3))
            virial[:,:] = stress_matrix(-stress*at.get_volume())
            at.params['virial'] = virial
        if at.calc_local_virial:
            stresses = at.get_stresses()
            at.add_property('local_virial', 0.0, n_cols=9, overwrite=True)
            lv = at.local_virial.view(np.ndarray)
            vol_per_atom = at.get_volume()/len(at)
            lv[...] = -stresses.T.reshape((len(at),9))*vol_per_atom

    return callback


class Potential(_potential.Potential):
    __doc__ = update_doc_string(_potential.Potential.__doc__, r"""
The :class:`Potential` class also implements the ASE
:class:`ase.calculators.interface.Calculator` interface via the
the :meth:`get_forces`, :meth:`get_stress`, :meth:`get_stresses`,
:meth:`get_potential_energy`, :meth:`get_potential_energies`
methods. This simplifies calculation since there is no need
to set the cutoff or to call :meth:`~quippy.atoms.Atoms.calc_connect`,
as this is done internally. The example above reduces to::

    atoms = diamond(5.44, 14)
    atoms.rattle(0.01)
    atoms.set_calculator(pot)
    forces = atoms.get_forces()
    print forces

Note that the ASE force array is the transpose of the QUIP force
array, so has shape (len(atoms), 3) rather than (3, len(atoms)).

The optional arguments `pot1`, `pot2` and `bulk_scale` are
used by ``Sum`` and ``ForceMixing`` potentials (see also
wrapper class :class:`ForceMixingPotential`)

An :class:`quippy.mpi_context.MPI_context` object can be
passed as the `mpi_obj` argument to restrict the
parallelisation of this potential to a subset of the

The `callback` argument is used to implement the calculation of
the :class:`Potential` in a Python function: see :meth:`set_callback` for
an example.

In addition to the builtin QUIP potentials, it is possible to
use any ASE calculator as a QUIP potential by passing it as
the `calculator` argument to the :class:`Potential` constructor, e.g.::

   from ase.calculators.morse import MorsePotential
   pot = Potential(calculator=MorsePotential)

`cutoff_skin` is used to set the :attr:`cutoff_skin` attribute.

`atoms` if given, is used to set the calculator associated
with `atoms` to the new :class:`Potential` instance, by calling
:meth:'.Atoms.set_calculator`.

.. note::

    QUIP potentials do not compute stress and per-atom stresses
    directly, but rather the virial tensor which has units of stress
    :math:`\times` volume, i.e. energy. If the total stress is
    requested, it is computed by dividing the virial by the atomic
    volume, obtained by calling :meth:`.Atoms.get_volume`. If per-atom
    stresses are requested, a per-atom volume is needed. By default
    this is taken to be the total volume divided by the number of
    atoms. In some cases, e.g. for systems containing large amounts of
    vacuum, this is not reasonable. The ``vol_per_atom`` calc_arg can
    be used either to give a single per-atom volume, or the name of an
    array in :attr:`.Atoms.arrays` containing volumes for each atom.

""",
    signature='Potential(init_args[, pot1, pot2, param_str, param_filename, bulk_scale, mpi_obj, callback, calculator, cutoff_skin, atoms])')

    callback_map = {}

    def __init__(self, init_args=None, pot1=None, pot2=None, param_str=None,
                 param_filename=None, bulk_scale=None, mpi_obj=None,
                 callback=None, calculator=None, cutoff_skin=1.0, atoms=None,
                 fpointer=None, finalise=True,
                 error=None, **kwargs):

        self.atoms = None
        self._prev_atoms = None
        self.energy = None
        self.energies = None
        self.forces = None
        self.stress = None
        self.stresses = None
        self.elastic_constants = None
        self.unrelaxed_elastic_constants = None
        self.numeric_forces = None
        self._calc_args = {}
        self._default_quantities = []
        self.cutoff_skin = cutoff_skin

        if callback is not None or calculator is not None:
            if init_args is None:
                init_args = 'callbackpot'

        if param_filename is not None:
            param_str = open(param_filename).read()

        if init_args is None and param_str is None:
            raise ValueError('Need one of init_args,param_str,param_filename')

        if init_args is not None:
            if init_args.lower().startswith('callbackpot'):
                if not 'label' in init_args:
                    init_args = init_args + ' label=%d' % id(self)
            else:
                # if param_str missing, try to find default set of QUIP params
                if param_str is None and pot1 is None and pot2 is None:
                    param_str = quip_xml_parameters(init_args)

        if kwargs != {}:
            if init_args is not None:
                init_args = init_args + ' ' + dict_to_args_str(kwargs)
            else:
                init_args = dict_to_args_str(kwargs)

        _potential.Potential.__init__(self, init_args, pot1=pot1, pot2=pot2,
                                         param_str=param_str,
                                         bulk_scale=bulk_scale,
                                         mpi_obj=mpi_obj,
                                         fpointer=fpointer, finalise=finalise,
                                         error=error)
        
        if init_args is not None and init_args.lower().startswith('callbackpot'):
            _potential.Potential.set_callback(self, Potential.callback)

            if callback is not None:
                self.set_callback(callback)

            if calculator is not None:
                self.set_callback(calculator_callback_factory(calculator))

        if atoms is not None:
            atoms.set_calculator(self)


    __init__.__doc__ = _potential.Potential.__init__.__doc__


    def calc(self, at, energy=None, force=None, virial=None,
             local_energy=None, local_virial=None,
             args_str=None, error=None, **kwargs):

        if not isinstance(args_str, basestring):
           args_str = dict_to_args_str(args_str)

        kw_args_str = dict_to_args_str(kwargs)

        args_str = ' '.join((self.get_calc_args_str(), kw_args_str, args_str))

        if isinstance(energy, basestring):
            args_str = args_str + ' energy=%s' % energy
            energy = None
        if isinstance(energy, bool) and energy:
            args_str = args_str + ' energy'
            energy = None

        if isinstance(force, basestring):
            args_str = args_str + ' force=%s' % force
            force = None
        if isinstance(force, bool) and force:
            args_str = args_str + ' force'
            force = None

        if isinstance(virial, basestring):
            args_str = args_str + ' virial=%s' % virial
            virial = None
        if isinstance(virial, bool) and virial:
            args_str = args_str + ' virial'
            virial = None

        if isinstance(local_energy, basestring):
            args_str = args_str + ' local_energy=%s' % local_energy
            local_energy = None
        if isinstance(local_energy, bool) and local_energy:
            args_str = args_str + ' local_energy'
            local_energy = None

        if isinstance(local_virial, basestring):
            args_str = args_str + ' local_virial=%s' % local_virial
            local_virial = None
        if isinstance(local_virial, bool) and local_virial:
            args_str = args_str + ' local_virial'
            local_virial = None

        potlog.debug('Potential invoking calc() on n=%d atoms with args_str "%s"' %
                     (len(at), args_str))
        _potential.Potential.calc(self, at, energy, force,
                                  virial, local_energy,
                                  local_virial,
                                  args_str, error)

    calc.__doc__ = update_doc_string(_potential.Potential.calc.__doc__,
       """In Python, this method is overloaded to set the final args_str to
          :meth:`get_calc_args_str`, followed by any keyword arguments,
          followed by an explicit `args_str` argument if present. This ordering
          ensures arguments explicitly passed to :meth:`calc` will override any
          default arguments.""")


    @staticmethod
    def callback(at_ptr):
        from quippy import Atoms
        at = Atoms(fpointer=at_ptr, finalise=False)
        if 'label' not in at.params or at.params['label'] not in Potential.callback_map:
            raise ValueError('Unknown Callback label %s' % at.params['label'])
        Potential.callback_map[at.params['label']](at)


    def set_callback(self, callback):
        """
        For a :class:`Potential` of type `CallbackPot`, this method is
        used to set the callback function. `callback` should be a Python
        function (or other callable, such as a bound method or class
        instance) which takes a single argument, of type
        :class:`~quippy.atoms.Atoms`. Information about which quantities should be
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
        """
        Potential.callback_map[str(id(self))] = callback


    def wipe(self):
        """
        Mark all quantities as needing to be recalculated
        """
        self.energy = None
        self.energies = None
        self.forces = None
        self.stress = None
        self.stresses = None
        self.numeric_forces = None
        self.elastic_constants = None
        self.unrelaxed_elastic_constants = None


    def update(self, atoms):
        """
        Set the :class:`~quippy.atoms.Atoms` object associated with this :class:`Potential` to `atoms`.

        Called internally by :meth:`get_potential_energy`,
        :meth:`get_forces`, etc.  Only a weak reference to `atoms` is
        kept, to prevent circular references.  If `atoms` is not a
        :class:`quippy.atoms.Atoms` instance, then a copy is made and a
        warning will be printed.
        """
        # we will do the calculation in place, to minimise number of copies,
        # unless atoms is not a quippy Atoms
        if isinstance(atoms, Atoms):
            self.atoms = weakref.proxy(atoms)
        else:
            potlog.debug('Potential atoms is not quippy.Atoms instance, copy forced!')
            self.atoms = Atoms(atoms)

        # check if atoms has changed since last call
        if self._prev_atoms is not None and self._prev_atoms.equivalent(self.atoms):
            return

        # Mark all quantities as needing to be recalculated
        self.wipe()

        # do we need to reinitialise _prev_atoms?
        if self._prev_atoms is None or len(self._prev_atoms) != len(self.atoms) or not self.atoms.connect.initialised:
            self._prev_atoms = Atoms()
            self._prev_atoms.copy_without_connect(self.atoms)
            self._prev_atoms.add_property('orig_pos', self.atoms.pos)
        else:
            # _prev_atoms is OK, update it in place
            self._prev_atoms.z[...] = self.atoms.z
            self._prev_atoms.pos[...] = self.atoms.pos
            self._prev_atoms.lattice[...] = self.atoms.lattice

        # do a calc_connect(), setting cutoff_skin so full reconnect will only be done when necessary
        self.atoms.set_cutoff(self.cutoff(), cutoff_skin=self.cutoff_skin)
        potlog.debug('Potential doing calc_connect() with cutoff %f cutoff_skin %r' % (self.atoms.cutoff, self.cutoff_skin))
        self.atoms.calc_connect()

    # Synonyms for `update` for compatibility with ASE calculator interface
    def initialize(self, atoms):
        self.update(atoms)

    def set_atoms(self, atoms):
        self.update(atoms)

    def calculation_required(self, atoms, quantities):
        self.update(atoms)
        for quantity in quantities:
            if getattr(self, quantity) is None:
                return True
        return False


    def calculate(self, atoms, quantities=None):
        """
        Perform a calculation of `quantities` for `atoms` using this Potential.

        Automatically determines if a new calculation is required or if previous
        results are still appliciable (i.e. if the atoms haven't moved since last call)
        Called internally by :meth:`get_potential_energy`, :meth:`get_forces`, etc.
        """
        if quantities is None:
            quantities = ['energy', 'forces', 'stress']

        # Add any default quantities
        quantities = set(self.get_default_quantities() + quantities)

        if len(quantities) == 0:
            raise RuntimeError('Nothing to calculate')

        if not self.calculation_required(atoms, quantities):
            return

        args_map = {
            'energy':          {'energy': None},
            'energies':        {'local_energy': None},
            'forces':          {'force': None},
            'stress':          {'virial': None},
            'numeric_forces':  {'force': 'numeric_force',
                                'force_using_fd': True,
                                'force_fd_delta': 1.0e-5},
            'stresses':        {'local_virial': None},
            'elastic_constants': {},
            'unrelaxed_elastic_constants': {}
            }

        # list of quantities that require a call to Potential.calc()
        calc_quantities = ['energy', 'energies', 'forces', 'numeric_forces', 'stress', 'stresses']

        # list of other quantities we know how to calculate
        other_quantities = ['elastic_constants', 'unrelaxed_elastic_constants']

        calc_args = {}
        calc_required = False
        for quantity in quantities:
            if quantity in calc_quantities:
                calc_required = True
                calc_args.update(args_map[quantity])
            elif quantity not in other_quantities:
                raise RuntimeError("Don't know how to calculate quantity '%s'" % quantity)

        if calc_required:
            self.calc(self.atoms, args_str=dict_to_args_str(calc_args))

        if 'energy' in quantities:
            self.energy = float(self.atoms.energy)
        if 'energies' in quantities:
            self.energies = self.atoms.local_energy.view(np.ndarray)
        if 'forces' in quantities:
            self.forces = self.atoms.force.view(np.ndarray).T
        if 'numeric_forces' in quantities:
            self.numeric_forces = self.atoms.numeric_force.view(np.ndarray).T
        if 'stress' in quantities:
            stress = -self.atoms.virial.view(np.ndarray)/self.atoms.get_volume()
            # convert to 6-element array in Voigt order
            self.stress = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                                    stress[1, 2], stress[0, 2], stress[0, 1]])
        if 'stresses' in quantities:
            lv = np.array(self.atoms.local_virial) # make a copy
            vol_per_atom = self.get('vol_per_atom', self.atoms.get_volume()/len(atoms))
            if isinstance(vol_per_atom, basestring):
                vol_per_atom = self.atoms.arrays[vol_per_atom]
            self.stresses = -lv.T.reshape((len(atoms), 3, 3), order='F')/vol_per_atom

        if 'elastic_constants' in quantities:
            cij_dx = self.get('cij_dx', 1e-2)
            cij = fzeros((6,6))
            self.calc_elastic_constants(self.atoms, fd=cij_dx,
                                        args_str=self.get_calc_args_str(),
                                        c=cij, relax_initial=False, return_relaxed=False)
            if not get_fortran_indexing():
                cij = cij.view(np.ndarray)
            self.elastic_constants = cij

        if 'unrelaxed_elastic_constants' in quantities:
            cij_dx = self.get('cij_dx', 1e-2)
            c0ij = fzeros((6,6))
            self.calc_elastic_constants(self.atoms, fd=cij_dx,
                                        args_str=self.get_calc_args_str(),
                                        c0=c0ij, relax_initial=False, return_relaxed=False)
            if not get_fortran_indexing():
                c0ij = c0ij.view(np.ndarray)
            self.unrelaxed_elastic_constants = c0ij


    def get_potential_energy(self, atoms):
        """
        Return potential energy of `atoms` calculated with this Potential
        """
        self.calculate(atoms, ['energy'])
        return self.energy


    def get_potential_energies(self, atoms):
        """
        Return array of atomic energies calculated with this Potential
        """
        self.calculate(atoms, ['energies'])
        return self.energies.copy()


    def get_forces(self, atoms):
        """
        Return forces on `atoms` calculated with this Potential
        """
        self.calculate(atoms, ['forces'])
        return self.forces.copy()


    def get_numeric_forces(self, atoms):
        """
        Return forces on `atoms` computed with finite differences of the energy
        """
        self.calculate(atoms, ['numeric_forces'])
        return self.numeric_forces.copy()


    def get_stress(self, atoms):
        """
        Return stress tensor for `atoms` computed with this Potential

        Result is a 6-element array in Voigt notation:
           [sigma_xx, sigma_yy, sigma_zz, sigma_yz, sigma_xz, sigma_xy]
        """
        self.calculate(atoms, ['stress'])
        return self.stress.copy()


    def get_stresses(self, atoms):
        """
        Return the per-atoms virial stress tensors for `atoms` computed with this Potential
        """
        self.calculate(atoms, ['stresses'])
        return self.stresses.copy()


    def get_elastic_constants(self, atoms):
        """
        Calculate elastic constants of `atoms` using this Potential.

        Returns  6x6 matrix :math:`C_{ij}` of elastic constants.

        The elastic contants are calculated as finite difference
        derivatives of the virial stress tensor using positive and
        negative strains of magnitude the `cij_dx` entry in
        ``calc_args``.
        """
        self.calculate(atoms, ['elastic_constants'])
        return self.elastic_constants.copy()


    def get_unrelaxed_elastic_constants(self, atoms):
        """
        Calculate unrelaxed elastic constants of `atoms` using this Potential

        Returns 6x6 matrix :math:`C^0_{ij}` of unrelaxed elastic constants.

        The elastic contants are calculated as finite difference
        derivatives of the virial stress tensor using positive and
        negative strains of magnitude the `cij_dx` entry in
        :attr:`calc_args`.
        """
        self.calculate(atoms, ['unrelaxed_elastic_constants'])
        return self.unrelaxed_elastic_constants.copy()


    def get_default_quantities(self):
        "Get the list of quantities to be calculated by default"
        return self._default_quantities[:]


    def set_default_quantities(self, quantities):
        "Set the list of quantities to be calculated by default"
        self._default_quantities = quantities[:]


    def get(self, param, default=None):
        """
        Get the value of a ``calc_args`` parameter for this :class:`Potential`

        Returns ``None`` if `param` is not in the current ``calc_args`` dictionary.

        All calc_args are passed to :meth:`calc` whenever energies,
        forces or stresses need to be re-computed.
        """
        return self._calc_args.get(param, default)


    def set(self, **kwargs):
        """
        Set one or more calc_args parameters for this Potential

        All calc_args are passed to :meth:`calc` whenever energies,
        forces or stresses need to be computed.

        After updating the calc_args, :meth:`set` calls :meth:`wipe`
        to mark all quantities as needing to be recaculated.
        """
        self._calc_args.update(kwargs)
        self.wipe()


    def get_calc_args(self):
        """
        Get the current ``calc_args``
        """
        return self._calc_args.copy()


    def set_calc_args(self, calc_args):
        """
        Set the ``calc_args`` to be used subsequent :meth:`calc` calls
        """
        self._calc_args = calc_args.copy()


    def get_calc_args_str(self):
        """
        Get the ``calc_args`` to be passed to :meth:`calc` as a string
        """
        return dict_to_args_str(self._calc_args)


    def get_cutoff_skin(self):
        return self._cutoff_skin

    def set_cutoff_skin(self, cutoff_skin):
        self._cutoff_skin = cutoff_skin
        self._prev_atoms = None # force a recalculation

    cutoff_skin = property(get_cutoff_skin, set_cutoff_skin,
                           doc="""
                           The `cutoff_skin` attribute is only relevant when the ASE-style
                           interface to the Potential is used, via the :meth:`get_forces`,
                           :meth:`get_potential_energy` etc. methods. In this case the
                           connectivity of the :class:`~quippy.atoms.Atoms` object for which
                           the calculation is requested is automatically kept up to date by
                           using a neighbour cutoff of :meth:`cutoff` + `cutoff_skin`, and
                           recalculating the neighbour lists whenever the maximum displacement
                           since the last :meth:`Atoms.calc_connect` exceeds `cutoff_skin`.
                           """)


from quippy import FortranDerivedTypes
FortranDerivedTypes['type(potential)'] = Potential


try:
    from ase.optimize import Optimizer
except ImportError:
    Optimizer = object


class Minim(Optimizer):
    """
    Minimise forces and/or virial tensor components wrt atomic
    positions and/or cell.

    This class is a wrapper around the
    :meth:`quippy.potential.Potential.minim` routine, compatible with
    the ASE `optimizer inferface <https://wiki.fysik.dtu.dk/ase/ase/optimize.html>`_.

    `method` should be one of ``"sd"``, (steepest descent), ``"cg"``
    (conjugate gradients, the default), ``"cg_n"`` (Noam Bernstein's conjugate gradients
    implementation), ``"pcg"``, (preconditioned conjugate gradients),
    ``"lbfgs"`` (L-BFGS), or ``"fire"`` (FIRE).

    Example usage::

       from quippy.structures import diamond
       from quippy.potential import Potential, Minim

       orig_atoms = diamond(5.44, 14)
       atoms = orig_atoms.copy()
       atoms.rattle(0.01)  # randomise the atomic positions a little

       potential = Potential('IP SW')
       atoms.set_calculator(potential)

       minimiser = Minim(atoms, relax_positions=True, relax_cell=False)
       minimiser.run(fmax=0.01)

       print orig_atoms.positions - atoms.positions # should be approx zero

    Note that if the :class:`~quippy.atoms.Atoms` object passed to the
    :class:`Minim` constructor is a :class:`quippy.atoms.Atoms`
    instance, the minimisation is done in place.  Otherwise, a copy is
    made, but the relaxed positions and cell vectors are copied back
    at the end of the :meth`run` method.

    """
    def __init__(self, atoms, restart=None,
                 relax_positions=True, relax_cell=True,
                 logfile='-', trajectory=None,
                 method='cg', linminroutine=None,
                 eps_guess=None,
                 fire_dt0=None, fire_dt_max=None,
                 external_pressure=None,
                 use_precond=None, cutoff_skin=1.0):

        calc = atoms.get_calculator()
        if calc is None:
            raise RuntimeError('Atoms object has no calculator')

        if not isinstance(calc, Potential):
            calc = Potential(calculator=calc)
        self.potential = calc

        # we will do the calculation in place, to minimise number of copies,
        # unless atoms is not a quippy Atoms
        self._atoms = atoms
        if not isinstance(atoms, Atoms):
            potlog.warn('Minim atoms is not quippy.Atoms instance, copy forced!')
            atoms = Atoms(atoms)
        self.atoms = atoms

        self.nsteps = 0
        if restart is not None:
            raise NotImplementedError
        if logfile != '-':
            raise NotImplementedError
        if trajectory is not None:
            raise NotImplementedError

        self.method = method
        self.linminroutine = linminroutine
        self.do_pos = relax_positions
        self.do_lat = relax_cell
        self.eps_guess = eps_guess
        self.fire_dt0 = fire_dt0
        self.fire_dt_max = fire_dt_max
        self.external_pressure = external_pressure
        self.use_precond = use_precond
        self.cutoff_skin = cutoff_skin


    def run(self, fmax=0.05, steps=100000000, convergence_tol=None):
        """
        Run the minimiser until maximum force is below `fmax`.

        Maximum number of minimisation steps is given by `steps`.

        Note that QUIP minim convergence criteria is actually based on
        :math:`|\mathbf{f}|^2` rather than on the maximum force. Here
        we convert using the relation::

                convergence_tol = 3*len(atoms)*fmax**2

        which is only approximately equivalent. To specify the converengence
        tolerance exactly, pass a value for the `convergence_tol` argument.
        """

        if convergence_tol is None:
            convergence_tol = 3*len(self.atoms)*fmax**2 # FIXME only approximately equivalent

        if self.atoms is not self._atoms:
            if self.do_pos:
                self.atoms.set_positions(self._atoms.get_positions())
            if self.do_lat:
                self.atoms.set_cell(self._atoms.get_cell())

        self.atoms.set_cutoff(self.potential.cutoff()+self.cutoff_skin)
        self.atoms.calc_connect()

        # check for constraints, only FixAtoms is supported for now
        if self.atoms.constraints:
            try:
                import ase.constraints
            except ImportError:
                raise RuntimeError('atoms has constraints but cannot import ase.constraints module')

        for constraint in self.atoms.constraints:
            if isinstance(constraint, ase.constraints.FixAtoms):
                self.atoms.add_property('move_mask', 1, overwrite=True)
                self.atoms.move_mask.view(np.ndarray)[constraint.index] = 0 # 0-based indices
            else:
                raise NotImplementedError('cannot convert ASE constraint %r to QUIP constraint' % constraint)

        potlog.debug('Minim.run() calling QUIP minim() with convergence_tol %f and args_str %s' %
                     (convergence_tol, self.potential.get_calc_args_str()))
        self.nsteps = self.potential.minim(self.atoms, self.method, convergence_tol, steps,
                                           self.linminroutine, do_pos=self.do_pos, do_lat=self.do_lat,
                                           args_str=self.potential.get_calc_args_str(), eps_guess=self.eps_guess,
                                           fire_minim_dt0=self.fire_dt0, fire_minim_dt_max=self.fire_dt_max,
                                           external_pressure=self.external_pressure,
                                           use_precond=self.use_precond)

        if self.atoms is not self._atoms:
            # update with results of minimisation
            if self.do_pos:
                self._atoms.set_positions(self.atoms.get_positions())
            if self.do_lat:
                self._atoms.set_cell(self.atoms.get_cell())



    def get_number_of_steps(self):
        """
        Return number of steps taken during minimisation
        """
        return self.nsteps


class ForceMixingPotential(Potential):
    """
    Subclass of :class:`Potential` for mixing forces from two Potentials
    """

    def __init__(self, pot1, pot2, bulk_scale=None, mpi_obj=None,
                 callback=None, calculator=None, cutoff_skin=1.0, atoms=None,
                 qm_list=None, fpointer=None, finalise=True, error=None, **kwargs):

        args_str = 'ForceMixing'
        Potential.__init__(self, args_str,
                           pot1=pot1, pot2=pot2, bulk_scale=bulk_scale,
                           mpi_obj=mpi_obj, cutoff_skin=cutoff_skin,
                           atoms=atoms, fpointer=fpointer, finalise=finalise,
                           error=error)
        if qm_list is not None:
            self.set_qm_atoms(qm_list)
        self.set(**kwargs)


    def get_qm_atoms(self):
        """
        Return the current list of QM atom indices as a list
        """
        if self.atoms is None:
            raise RuntimeError('No atoms assocated with this ForceMixingPotential!')
        return list((self.atoms.hybrid == HYBRID_ACTIVE_MARK).nonzero()[0])


    def set_qm_atoms(self, qm_list):
        """
        Set the QM atoms, given as a list of atom indices
        """
        if self.atoms is None:
            raise RuntimeError('No atoms assocated with this ForceMixingPotential!')
        if not self.atoms.has_property('hybrid'):
            self.atoms.add_property('hybrid', HYBRID_NO_MARK)
        self.atoms.hybrid[:] = HYBRID_NO_MARK 
        self.atoms.hybrid[qm_list] = HYBRID_ACTIVE_MARK
         


def force_test(at, p, dx=1e-4):
    """
    Compare analyric and numeric forces for the Potential `p` with Atoms `at`

    Finite difference derivates are calculated by moving each atom by `dx`.
    """
    analytic_f = fzeros((3,at.n))
    p.calc(at, force=analytic_f)
    num_f = fzeros((3,at.n))
    ep, em = farray(0.0), farray(0.0)

    for i in frange(at.n):
        for j in (1,2,3):
            ap = at.copy()
            ap.pos[j,i] += dx
            p.calc(ap, energy=ep)
            print 'e+', j,i,ep
            ap.pos[j,i] -= 2.0*dx
            p.calc(ap, energy=em)
            print 'e-', j,i,em
            num_f[j,i] = -(ep - em)/(2*dx)

    return analytic_f, num_f, analytic_f - num_f

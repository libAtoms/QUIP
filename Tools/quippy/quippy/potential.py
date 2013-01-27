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
from quippy.oo_fortran import update_doc_string
from quippy.atoms import Atoms
from quippy.units import GPA
from quippy.util import quip_xml_parameters, dict_to_args_str
from quippy.elasticity import stress_matrix
from quippy.farray import fzeros

__doc__ = _potential.__doc__
__all__ = _potential.__all__ + ['force_test', 'Minim']

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
            virial[:,:] = stress_matrix(-stress*at.get_volume()/GPA)
            at.params['virial'] = virial
        if at.calc_local_virial:
            stresses = at.get_stresses()
            at.add_property('local_virial', 0.0, n_cols=9, overwrite=True)
            lv = at.local_virial.view(np.ndarray)
            vol = at.get_volume()
            for i, stress in enumerate(stresses):
                lv[:,i] = -(stress*vol/GPA).reshape((9,), order='F')

    return callback


class Potential(_potential.Potential):
    __doc__ = update_doc_string(_potential.Potential.__doc__, """
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
wrapper classes :class:`ForceMixingPotential` and
:class:`LOTFPotential`).

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
""", signature='Potential(args_str[, pot1, pot2, param_str, param_filename, bulk_scale, mpi_obj, callback, calculator, cutoff_skin, fortran_indexing])')

    callback_map = {}

    def __init__(self, args_str=None, pot1=None, pot2=None, param_str=None,
                 param_filename=None, bulk_scale=None, mpi_obj=None,
                 callback=None, calculator=None, cutoff_skin=1.0,
                 fortran_indexing=True, fpointer=None, finalise=True,
                 error=None, **kwargs):

        self.atoms = None
        self._prev_atoms = None
        self.energy = None
        self.energies = None
        self.forces = None
        self.stress = None
        self.stresses = None
        self.numeric_forces = None
        self._calc_args = {}
        self._default_quantities = []
        self.cutoff_skin = cutoff_skin

        if callback is not None or calculator is not None:
            if args_str is None:
                args_str = 'callbackpot'

        if param_filename is not None:
            param_str = open(param_filename).read()

        if args_str is None and param_str is None:
            raise ValueError('Need one of args_str,param_str,param_filename')

        if args_str is None:
            xml_label = param_str[:param_str.index('\n')].translate(None,'<>')
            args_str = 'xml_label=%s' % xml_label

        if args_str.lower().startswith('callbackpot'):
            if not 'label' in args_str:
                args_str = args_str + ' label=%d' % id(self)
        else:
            # if param_str missing, try to find default set of QUIP params
            if param_str is None and pot1 is None and pot2 is None:
                param_str = quip_xml_parameters(args_str)

        if kwargs != {}:
            args_str = args_str + ' ' + dict_to_args_str(kwargs)

        _potential.Potential.__init__(self, args_str, pot1=pot1, pot2=pot2,
                                         param_str=param_str,
                                         bulk_scale=bulk_scale,
                                         mpi_obj=mpi_obj,
                                         fortran_indexing=fortran_indexing,
                                         fpointer=fpointer, finalise=finalise,
                                         error=error)

        if args_str.lower().startswith('callbackpot'):
            _potential.Potential.set_callback(self, Potential.callback)

            if callback is not None:
                self.set_callback(callback)

            if calculator is not None:
                self.set_callback(calculator_callback_factory(calculator))


    __init__.__doc__ = _potential.Potential.__init__.__doc__

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
            potlog.warn('Potential atoms is not quippy.Atoms instance, copy forced!')
            self.atoms = Atoms(atoms)

        # check if atoms has changed since last call
        if self._prev_atoms is not None and self._prev_atoms.equivalent(self.atoms):
            return

        # mark all quantities as needing to be recalculated
        self.energy = None
        self.energies = None
        self.forces = None
        self.stress = None
        self.stresses = None
        self.numeric_forces = None

        # do we need to reinitialise _prev_atoms?
        if self._prev_atoms is None or len(self._prev_atoms) != len(self.atoms) or not self.atoms.connect.initialised:
            self._prev_atoms = Atoms()
            self._prev_atoms.copy_without_connect(self.atoms)
            self._prev_atoms.add_property('orig_pos', self.atoms.pos)
            self.atoms.set_cutoff(self.cutoff() + self.cutoff_skin)
            potlog.debug('Potential doing initial calc_connect() with cutoff %f' % self.atoms.cutoff)
            self.atoms.calc_connect()
        else:
            # _prev_atoms is OK, update it in place
            self._prev_atoms.z[...] = self.atoms.z
            self._prev_atoms.pos[...] = self.atoms.pos
            self._prev_atoms.lattice[...] = self.atoms.lattice

            self._prev_atoms.undo_pbc_jumps(persistent=True)
            self._prev_atoms.calc_msd(persistent=True)
            max_disp_sq = self._prev_atoms.msd_displ_mag.max()
            potlog.debug('Potential max_disp %f' % sqrt(max_disp_sq))

            if max_disp_sq > self.cutoff_skin**2:
                # do a calc_connect() if any atoms have moved more than cutoff_skin
                potlog.debug('Potential calling atoms.calc_connect()')
                self.atoms.calc_connect()
                self._prev_atoms.orig_pos[...] = self.atoms.pos
            else:
                # just update the neighbour distance tables in place
                potlog.debug('Potential calling atoms.calc_dists()')
                self.atoms.calc_dists()

           

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

        Automatically determines if a new calculation is required of if previous
        results are still appliciable (i.e. if the atoms haven't moved since last call)
        Called internally by :meth:`get_potential_energy`, :meth:`get_forces`, etc.
        """
        if quantities is None:
            quantities = ['energy', 'forces', 'stress']

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
            'stresses':        {'local_virial': None}
            }

        calc_args = self.get_calc_args()
        for quantity in self.get_default_quantities() + quantities:
            calc_args.update(args_map[quantity])

        potlog.debug('Potential calling calc() on atoms with n=%d with args %s' % (len(atoms), calc_args))
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
            self.stress = -(self.atoms.virial.view(np.ndarray)*
                            GPA/self.atoms.get_volume())
        if 'stresses' in quantities:
            self.stresses = np.zeros((len(self.atoms), 3, 3))

            lv = self.atoms.local_virial.view(np.ndarray)
            vol = self.atoms.get_volume()
            for i in range(len(self.atoms)):
                self.stresses[i,:,:] = -(lv[:,i].reshape((3,3),order='F')*
                                         GPA/vol)


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
        Return virial stress tensor for `atoms` computed with this Potential
        """
        self.calculate(atoms, ['stress'])
        return self.stress.copy()


    def get_stresses(self, atoms):
        """
        Return the per-atoms virial stress tensors for `atoms` computed with this Potential
        """
        self.calculate(atoms, ['stresses'])
        return self.stresses.copy()
    

    def get_elastic_constants(self, atoms, fd=0.01, relaxed=True,
                              relax_initial=True, relax_tol=None,
                              relax_method=None):
        """
        Calculate elastic constants of `atoms` using this Potential.

        Returns 6x6 matrix :math:`C_{ij}`, or if `relaxed` is False,
        the 6x6 matrix :math:`C^0_{ij}` of unrelaxed elastic constants.

        The elastic contants are calculated as finite difference
        derivatives of the virial stress tensor using positive and
        negative strains of magnitude `fd`.
        
        If `relax_initial` is True, the positions and cell of `atoms`
        are relaxed before calculating the elastic constants (note
        that this is done on a copy, so `atoms` is not modified).

        `relax_tol` and `relax_method` are the minimisation tolerance
        and method to use both for the relaxation of the input
        structure, and for the relaxations of the internal coordinates
        at each applied strain (if `relaxed` is True, which is the
        default). See :meth:`minim` for full details.
        """
        c = c0 = None
        if relaxed:
            c = fzeros((6,6))
        else:
            c0 = fzeros((6,6))
        self.calc_elastic_constants(atoms, fd=fd, args_str=self.get_calc_args_str(),
                                    c=c, c0=c0, relax_initial=relax_initial, return_relaxed=False,
                                    relax_tol=relax_tol, relax_method=relax_method)
        if relaxed:
            if not self.fortran_indexing:
                c = c.view(np.ndarray)
            return c
        else:
            if not self.fortran_indexing:
                c0 = c0.view(cp.ndarray)
            return c0

    def get_default_quantities(self):
        return self._default_quantities[:]

    def set_default_quantities(self, quantities):
        self._default_quantities = quantities[:]

    default_quantities = property(get_default_quantities,
                                  set_default_quantities,
                                  doc="Get or set the list of quantities to be calculated by default")
    

    def get_calc_args(self):
        return self._calc_args.copy()

    def set_calc_args(self, calc_args):
        self._calc_args = calc_args

    calc_args = property(get_calc_args, set_calc_args,
                         doc="Dictionary of extra arguments to be passed to Potential.calc() routine")

    def get_calc_args_str(self):
        """
        Return the extra arguments to pass to Potential.calc() as a string
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
                           connectivity of the :class:`~quippy.atoms.Atoms` object for which the calculation
                           is requested is automatically kept up to date by using a neighbour
                           cutoff of :meth:`cutoff` + `cutoff_skin`, and recalculating the 
                           neighbour lists whenever the maximum displacement since the last
                           :meth:`Atoms.calc_connect` exceeds `cutoff_skin`.
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
    the ASE `ASE Optimizer inferface <https://wiki.fysik.dtu.dk/ase/ase/optimize.html>`_.

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
        if not isinstance(calc, Potential):
            calc = Potential(calculator=calc)
        self.potential = calc

        # we will do the calculation in place, to minimise number of copies,
        # unless atoms is not a quippy Atoms
        if not isinstance(atoms, Atoms):
            self._atoms = atoms
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
                _save = self.fortran_indexing
                self.fortran_indexing = False
                self.move_mask[constraint.index] = 0 # 0-based indices
                self.fortran_indexing = _save
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
    Subclass of Potential specifically for mixing forces from two Potentials
    """

    def __init__(self, pot1, pot2, args_str=None, param_str=None,
                 param_filename=None, bulk_scale=None, mpi_obj=None,
                 callback=None, calculator=None, 
                 fortran_indexing=True, fpointer=None, finalise=True,
                 cutoff_skin=1.0, error=None, **kwargs):
        
        if args_str is None:
            args_str = 'ForceMixing'
        Potential.__init__(self, args_str,
                           pot1=pot1, pot2=pot2, param_str=param_str,
                           mpi_obj=mpi_obj, fortran_indexing=fortran_indexing,
                           fpointer=fpointer, finalise=finalise,
                           cutoff_skin=cutoff_skin,
                           error=error, **kwargs)


class LOTFPotential(ForceMixingPotential):
    """
    Subclass of ForceMixing specificically for 'Learn on the Fly' scheme
    """

    def __init__(self, pot1, pot2, args_str=None, param_str=None,
                 param_filename=None, bulk_scale=None, mpi_obj=None,
                 callback=None, calculator=None, 
                 fortran_indexing=True, fpointer=None, finalise=True,
                 cutoff_skin=None, error=None, **kwargs):

        if args_str is None:
            args_str = 'ForceMixing method=lotf_adj_pot_svd %s'
        ForceMixingPotential.__init__(self, args_str=args_str,
                                      pot1=pot1, pot2=pot2, param_str=param_str,
                                      mpi_obj=mpi_obj, fortran_indexing=fortran_indexing,
                                      fpointer=fpointer, finalise=finalise,
                                      cutoff_skin=cutoff_skin,
                                      error=error, **kwargs)
    


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



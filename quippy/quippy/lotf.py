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

"""
Module containing :class:`LOTFDynamics` integrator for doing 'Learn on the Fly'
predictor-corrector dynamics within the ASE molecular dynamics framework.
"""

import numpy as np

from quippy.clusters import (HYBRID_ACTIVE_MARK, HYBRID_NO_MARK,
                             construct_hysteretic_region,
                             create_hybrid_weights,
                             create_cluster_simple)
from quippy.potential import Potential, ForceMixingPotential
from quippy.table import Table
from quippy.system import system_get_random_seed, system_set_random_seeds, system_timer
from quippy.util import args_str


__all__ = ['LOTFDynamics', 'update_hysteretic_qm_region', 'iter_atom_centered_clusters']


from ase.md.md import MolecularDynamics

class LOTFDynamics(MolecularDynamics):
    """
    'Learn on the Fly' dynamics with extrapolation and interpolation

    Subclass of :class:`ase.md.md.MolecularDynamics` which
    implements LOTF predictor/corrector dynamics. See the :meth:`step`
    method for a description of the algorithm.

    The calculator associated with `Atoms` should be a
    :class:`quippy.potential.Potential` object.
    """

    Extrapolation = 0
    Interpolation = 1

    def __init__(self, atoms, timestep, extrapolate_steps,
                 trajectory=None, logfile=None, loginterval=1,
                 check_force_error=False, qm_update_func=None):

        MolecularDynamics.__init__(self, atoms, timestep,
                                   trajectory, logfile, loginterval)
        self.extrapolate_steps = extrapolate_steps
        self.nsteps = 0
        self.saved_state = None
        self.state = LOTFDynamics.Extrapolation
        self.check_force_error = check_force_error
        self.calc = self.atoms.get_calculator()
        self.set_qm_update_func(qm_update_func)
        if not isinstance(self.calc, Potential):
           raise RuntimeError('LOTFDynamics only works with quippy.potential.Potential calculators!')

        self.initialise_adjustable_potential(self.atoms)


    def set_qm_update_func(self, qm_update_func):
        """
        Specify a function that will be called up update the QM region

        `qm_update_func` should take a single argument, which is the
        :class:`~.Atoms` object representing the current state of the
        dynamical system.
        """
        self.qm_update_func = qm_update_func


    def get_time(self):
        """
        Return current time (in ASE time units)
        """
        return self.nsteps*self.dt


    def get_number_of_steps(self):
        """
        Return number of steps taken since start of dynamics
        """
        return self.nsteps


    def advance_verlet(self, f):
        """
        Single advance using velocity Verlet with forces `f`

        Do not call this method directly: it is invoked by :meth:`step`
        once for each extrapolation and interpolation step.

        Note that order of the two half-steps is inverted, to fit in
        with the predictor-corrector loops. This means the momenta
        are half-stepped outside of this routine.
        """
        # on entry we have p(t-dt/2), r(t), F(t)

        #  p(t) = p(t-dt/2) + F(t)/2
        self.atoms.set_momenta(self.atoms.get_momenta() + 0.5*self.dt*f)

        p = self.atoms.get_momenta()
        p += 0.5 * self.dt * f
        
        # r(t+dt) = r(t) + p(t+dt/2)/m dt
        self.atoms.set_positions(self.atoms.get_positions() +
                            self.dt * p / self.atoms.get_masses()[:,np.newaxis])

        # p(t+dt/2) = p(t) + F(t) dt/2
        self.atoms.set_momenta(p)

        # now we have p(t+dt/2), r(t+dt)


    @property
    def state_label(self):
        """
        One-character label for current dynamical state: either "E" or "I"
        for extrapolation or interpolation, respectively.
        """
        if self.state == LOTFDynamics.Extrapolation:
            return 'E'
        elif self.state == LOTFDynamics.Interpolation:
            return 'I'
        else:
            raise RuntimeError('Unknown state %s' % self.state)


    def get_extrapolate_steps(self):
        return self._extrapolate_steps

    def set_extrapolate_steps(self, extrapolate_steps):
        self._extrapolate_steps = extrapolate_steps

    extrapolate_steps = property(get_extrapolate_steps, set_extrapolate_steps,
                                 doc="""
                                 Number of small time-steps used when extrapolating
                                 or interpolating in the LOTF predictor-corrector loop.
                                 See :meth:`step` for a more detailed description.
                                 """)

    def get_state(self):
        """
        Save dynamical state, including random number seeds
        """
            
        saved_positions = self.atoms.get_positions()
        saved_momenta = self.atoms.get_momenta()
        saved_nsteps = self.nsteps
        saved_seed = system_get_random_seed()
        saved_np_random_state = np.random.get_state()
        saved_state = (saved_positions, saved_momenta, saved_nsteps,
                       saved_seed, saved_np_random_state)
        return saved_state
                       

    def set_state(self, saved_state):
        """
        Restore dynamics to a previously saved state.
        Sets arrays directly, ignoring any constraints.
        """

        (saved_positions, saved_momenta, saved_nsteps,
         saved_seed, saved_np_random_state) = saved_state
        self.atoms.set_array('positions', saved_positions)
        self.atoms.set_array('momenta', saved_momenta)
        self.nsteps = saved_nsteps
        system_set_random_seeds(saved_seed)
        np.random.set_state(saved_np_random_state)


    def step(self):
        """
        Advance the dynamics by one LOTF cycle

        The *Learn on the Fly* algorithm proceeds as follows:

           1. **QM selection** using the function registered
              with the :meth:`set_qm_update_func` method.

           2. **Save state** of the dynamical system using the
              :meth:`get_state` method.

           3. **Extrapolation** of the dynamics with fixed parameters for
              :attr:`extrapolate_steps`, each of size :attr:`dt`.
              This is the *predictor* part of the predictor-corrector cycle.
              The adjustable potential parameters are remapped if
              the set of QM atoms has changed since the most recent fit.

           4. **QM force calculation** and **optimisation** of the adjustable
              potential parameters to reproduce the QM forces at the fit point.
              With the linear 'spings' method, the fit can be reposed as a linear
              algebra problem and solved via singular value decomposition (SVD).
              See the adjustable potential code for full details in
              :svn:`QUIP_Core/AdjustablePotential.f95`.

           5. **Restore state** of the dynamical system saved at step 2,
              returning to before the extrapolation phase.

           6. **Interpolation** of the dynamics, with the adjustable potential
              parameters interpolated between their values at the two fit points.

        """

        system_timer('lotf_step')
        system_timer('lotf_extrapolate')

        if self.qm_update_func is not None:
            self.qm_update_func(self.atoms)

        saved_state = self.get_state()
        self.state = LOTFDynamics.Extrapolation

        for i in range(self.extrapolate_steps):
            # Compute extrapolation forces
            f = self.get_extrapolation_forces(self.atoms, i)

            # Possibly check force error
            if self.check_force_error:
                self._calc_force_error(f, i)

            # Apply any constraints to the forces
            for constraint in self.atoms.constraints:
                constraint.adjust_forces(self.atoms.arrays['positions'], f)            

            # Do the velocity Verlet step
            self.advance_verlet(f)
            self.nsteps += 1

            self.call_observers()
        system_timer('lotf_extrapolate')

        # Compute QM forces and fit the LOTF adjustable potential to them
        system_timer('lotf_qm_fit')
        self.fit_forces(self.atoms)
        system_timer('lotf_qm_fit')

        system_timer('lotf_interpolate')
        self.set_state(saved_state)
        self.state = LOTFDynamics.Interpolation

        for i in range(self.extrapolate_steps):
            # Compute interpolation forces
            f = self.get_interpolation_forces(self.atoms, i)

            # Possibly check force error
            if self.check_force_error:
               self._calc_force_error(f, i)

            # Apply any constraints to the forces
            for constraint in self.atoms.constraints:
                constraint.adjust_forces(self.atoms.arrays['positions'], f)

            # Do the velocity Verlet step
            self.advance_verlet(f)
            self.nsteps += 1            

            self.call_observers()
        system_timer('lotf_interpolate')
        system_timer('lotf_step')


    def run(self, steps):
        """
        Run LOTF dynamics forwards for given number of small time steps

        Number of LOTF cycles is given by ``step/extrapolate_steps``.
        """
        lotf_steps = max(1, steps/self.extrapolate_steps)
        for i in range(lotf_steps):
            self.step()
   

    def initialise_adjustable_potential(self, atoms):
        """
        Bootstrap the potential by doing a QM calculation and fitting the adj. pot.

        Returns the initial QM/MM forces on the atoms
        """
        self.calc.set(method='lotf_adj_pot_svd',
                      calc_weights=True,
                      lotf_do_qm=True,
                      lotf_do_init=True,
                      lotf_do_fit=True,
                      lotf_do_interp=False)
        return self.calc.get_forces(atoms)


    def get_extrapolation_forces(self, atoms, i):
        """
        Get LOTF forces at extrapolation step `i`
        """
        if i == 0:
            self.calc.set(calc_weights=True,
                          lotf_do_init=True,
                          lotf_do_map=True)
        else:
            self.calc.set(calc_weights=False,
                          lotf_do_init=False,
                          lotf_do_map=False)

        self.calc.set(method='lotf_adj_pot_svd',
                      lotf_do_qm=False,
                      lotf_do_fit=False,                 
                      lotf_do_interp=False)
        return self.calc.get_forces(atoms)


    def get_interpolation_forces(self, atoms, i):
        """
        Get LOTF forces at interpolation step `i`
        """
        if i == 0:
            self.calc.set(calc_weights=True)
        else:
            self.calc.set(calc_weights=False)
        interp_frac = float(i)/float(self.extrapolate_steps)
        self.calc.set(method='lotf_adj_pot_svd',
                      lotf_do_qm=False,
                      lotf_do_init=False,
                      lotf_do_interp=True,
                      lotf_do_fit=False,
                      lotf_interp=interp_frac)
        return self.calc.get_forces(atoms)


    def fit_forces(self, atoms):
        """
        Evaluate QM forces and optimise the adjustable potential to reproduce them

        Returns the QM forces on the atoms at the fit point
        """
        self.calc.set(method='lotf_adj_pot_svd',
                      calc_weights=True,
                      lotf_do_qm=True,
                      lotf_do_init=False,
                      lotf_do_map=False,
                      lotf_do_fit=True,
                      lotf_do_interp=False)
        return self.calc.get_forces(atoms)


    def get_reference_forces(self, atoms, i):
        """
        Do a reference QM calculation and return forces resulting from standard force mixing.

        Useful for checking errors during predictor/corrector extrapolation and interpolation.
        Does not change the current set of adjustable potential parameters.
        """
        if i == 0:
            self.calc.set(calc_weights=True)
        else:
            self.calc.set(calc_weights=False)
        self.calc.set(method='force_mixing',
                      lotf_do_qm=True,
                      lotf_do_init=False,
                      lotf_do_map=False,
                      lotf_do_fit=False,
                      lotf_do_interp=False)
        return self.calc.get_forces(atoms)


    def _calc_force_error(self, f, i, qm_only=True):
       """
       Compute reference forces and compare with current forces `f`
       """
       if qm_only:
           mask = self.calc.get_qm_atoms()
       else:
           mask = Ellipsis

       ref_forces = self.get_reference_forces(self.atoms, i)
       force_error = ref_forces[mask,:] - f[mask,:]

       self.atoms.set_array('ref_forces', ref_forces)
       self.atoms.set_array('force_error', ref_forces - f)
       
       self.rms_force_error = np.sqrt((force_error**2).mean())
       self.max_force_error = abs(force_error).max()



def update_hysteretic_qm_region(atoms, old_qm_list, qm_centre, qm_inner_radius,
                                qm_outer_radius, use_avgpos=False, update_marks=True):
    """
    Update the QM region in `atoms`

    Calls :meth:`~.ForceMixingPotential.set_qm_atoms` to mark the atoms
    if the Atoms' calculator is an instance of :class:`~.ForceMixingPotential`.

    Parameters
    ----------
    old_qm_list : list
       List of QM atoms at previous update
    qm_centre : array with shape (3,)
       Position of the new centre of the QM region
    qm_inner_radius : float
       Atoms which come within this distance of any core atom will become
       selected for QM treatment.
    qm_outer_radius : float
       Atoms stay QM until they are more than this distance from any core atom
    update_marks : bool
       If true (default), call :meth:`~.ForceMixingPotential.set_qm_atoms` to mark the atoms

    Returns
    -------
    qm_list : list
       List of the new QM atom indices
    """

    qm_table = Table.from_atom_list(atoms, old_qm_list)
    construct_hysteretic_region(qm_table, atoms, inner_radius=qm_inner_radius,
                                outer_radius=qm_outer_radius,
                                centre=qm_centre, use_avgpos=use_avgpos,
                                loop_atoms_no_connectivity=True)
    qm_list = qm_table.to_atom_list()
    
    print('update_qm_region: QM region with %d atoms centred on %s' %
          (len(qm_list), qm_centre))

    if update_marks:
        calc = atoms.get_calculator()
        if isinstance(calc, ForceMixingPotential):
            calc.set_qm_atoms(qm_list)

    return qm_list


def iter_atom_centered_clusters(at, mark_name='hybrid_mark', **cluster_args):
    """
    Iterate over all atom-centered (little) clusters in `at`.

    If `at` has a property with name `mark_name` (default `"hybrid_mark"`),
    only those atoms where ``hybrid_mark == HYBRID_ACTIVE_MARK`` are
    included. Otherwise, all atoms are included.

    Clusters are constructed by calling
    :func:`~quippy.clusters.create_hybrid_weights` and
    :func:`~quippy.clusters.create_cluster_simple` once for each atom
    of interest, after setting ``hybrid_mark=HYBRID_ACTIVE_MARK`` for that
    atom only.

    Any keyword arguments given are passed along to both cluster
    creation functions.
    """

    saved_hybrid_mark = None
    if hasattr(at, 'hybrid_mark'):
        saved_hybrid_mark = at.hybrid_mark.copy()

    if hasattr(at, mark_name):
        indices = (getattr(at, mark_name) == HYBRID_ACTIVE_MARK).nonzero()[0]
    else:
        indices = at.indices
    
    at.add_property('hybrid_mark', 0, overwrite=True)

    for i in indices:
        at.hybrid_mark[:] = 0
        at.hybrid_mark[i] = True
        create_hybrid_weights(at, args_str=args_str(cluster_args))
        c = create_cluster_simple(at, args_str=args_str(cluster_args))
        yield c

    if saved_hybrid_mark is not None:
        at.add_property('hybrid_mark', saved_hybrid_mark, overwrite=True)
    else:
        del at.properties['hybrid_mark']

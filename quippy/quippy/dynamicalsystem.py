# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Portions of this code were written by
# HQ X     Tamas K. Stenczel, James Kermode
# HQ X
# HQ X   Copyright 2019
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   https://warwick.ac.uk/fac/sci/eng/staff/jrk
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import warnings
from math import sqrt

import numpy as np

from quippy import dynamicalsystem_module
import quippy.convert
import quippy.atoms_types_module

import quippy._quippy as _quippy

import f90wrap.runtime

# from quippy._dynamicalsystem import *
# from quippy.units import MASSCONVERT
# taking it explicitly now:
# Units.f95:68 MASSCONVERT from old quippy
MASSCONVERT = 103.6426957074462
# from quippy.atoms import Atoms
# from quippy.farray import farray
# from quippy.io import AtomsWriter
# from quippy.system import InOutput, OUTPUT

import ase.io
import ase.md
from ase.io.trajectory import Trajectory
from ase.optimize import optimize

# time conversion units from ase.units
from ase.units import fs

# thermostat and barostat type constants are defined in Thermostats.f95,
# which is not directly wrapped by quippy, so we reproduce them here

THERMOSTAT_NONE = 0
THERMOSTAT_LANGEVIN = 1
THERMOSTAT_NOSE_HOOVER = 2
THERMOSTAT_NOSE_HOOVER_LANGEVIN = 3
THERMOSTAT_LANGEVIN_NPT = 4
THERMOSTAT_LANGEVIN_PR = 5
THERMOSTAT_NPH_ANDERSEN = 6
THERMOSTAT_NPH_PR = 7
THERMOSTAT_LANGEVIN_OU = 8
THERMOSTAT_LANGEVIN_NPT_NB = 9
THERMOSTAT_ALL_PURPOSE = 10

_thermostat_types = {
    'THERMOSTAT_NONE': 0,
    'THERMOSTAT_LANGEVIN': 1,
    'THERMOSTAT_NOSE_HOOVER': 2,
    'THERMOSTAT_NOSE_HOOVER_LANGEVIN': 3,
    'THERMOSTAT_LANGEVIN_NPT': 4,
    'THERMOSTAT_LANGEVIN_PR': 5,
    'THERMOSTAT_NPH_ANDERSEN': 6,
    'THERMOSTAT_NPH_PR': 7,
    'THERMOSTAT_LANGEVIN_OU': 8,
    'THERMOSTAT_LANGEVIN_NPT_NB': 9,
    'THERMOSTAT_ALL_PURPOSE': 10
}

BAROSTAT_NONE = 0
BAROSTAT_HOOVER_LANGEVIN = 1

_barostat_types = {
    'BAROSTAT_NONE': 0,
    'BAROSTAT_HOOVER_LANGEVIN': 1,
}

__all__ = ['Dynamics', 'DynamicalSystem']


class DynamicalSystem(dynamicalsystem_module.DynamicalSystem):
    __doc__ = dynamicalsystem_module.DynamicalSystem.__doc__

    # def __init__(self, atoms_in, velocity=None, acceleration=None, constraints=None,
    #     restraints=None, rigidbodies=None, error=None, handle=None):
    #
    #
    #     quip_atoms_in = quippy.convert.ase_to_quip(atoms_in)
    #
    #     dynamicalsystem_module.DynamicalSystem.__init__(self, atoms_in=quip_atoms_in, velocity=velocity,
    #                                                     acceleration=acceleration, constraints=constraints,
    #                                                     restraints=restraints, rigidbodies=rigidbodies, error=error,
    #                                                     handle=handle)

    def run(self, pot, dt, n_steps, summary_interval=None, hook_interval=None, write_interval=None,
            trajectory=None, args_str=None, hook=None,
            save_interval=None):

        if hook is None and hook_interval is not None:
            raise ValueError('hook_interval not permitted when hook is not present')

        if hook is None:
            traj = []
            save_hook = lambda: traj.append(self.atoms.copy())
            dynamicalsystem_module.DynamicalSystem.run(self, pot, dt, n_steps,
                                                       save_hook, hook_interval=save_interval,
                                                       summary_interval=summary_interval,
                                                       write_interval=write_interval,
                                                       trajectory=trajectory,
                                                       args_str=args_str)
            return traj
        else:
            dynamicalsystem_module.DynamicalSystem.run(self, pot, dt, n_steps, hook, hook_interval=hook_interval,
                                                       summary_interval=summary_interval, write_interval=write_interval,
                                                       trajectory=trajectory, args_str=args_str)

    #run.__doc__ = dynamicalsystem_module.DynamicalSystem.run.__doc__


class Dynamics(optimize.Dynamics):
    """
    Wrapper around :class:`DynamicalSystem` integrator compatible with
    ASE :class:`ase.md.MolecularDynamics` interface.

    .. note::
         This class uses ASE units for time (and hence velocities
         and momenta) which differ from the femtoseconds used in QUIP.
         When using this class, all times should be given in ASE time units.

    :param atoms: quippy or ASE Atoms instance. Note that if atoms is
                 not a quippy Atoms object, a copy will be made.

    :param timestep: in ASE time units (multiply by :const:`ase.units.fs`
                    to convert from femtoseconds)

    :param trajectory: output file to which to write the trajectory.
                      Can be a string to create a new output object.

    :param trajectoryinterval: interval at which to write frames

    :param initialtemperature: if not ``None``, rescale initial velocities
                              to a specified temperature, given in K.

    :param logfile: filename or open InOutput object to write log lines to.
              Use special filename ``"-"`` for stdout (default).

    :param loginterval: interval at which to write log lines
    """

    def __init__(self, atoms, timestep, trajectory,
                 trajectoryinterval=10, initialtemperature=None,
                 logfile='-', loginterval=1, loglabel='D'):

        # check for type first
        if not isinstance(atoms, ase.Atoms):
            atoms = ase.Atoms(atoms)

        # construct the quip atoms object which we will use to calculate on
        self.ase_atoms = atoms

        if not atoms.has('momenta'):  # so that there is a velocity initialisation on the quip object
            atoms.set_momenta(np.zeros((len(atoms), 3)))
        self._quip_atoms = quippy.convert.ase_to_quip(self.ase_atoms)

        # add the mass separately, because converter is not doing it
        _quippy.f90wrap_atoms_add_property_real_a(this=self._quip_atoms._handle, name='mass',
                                                  value=self.ase_atoms.get_masses() * MASSCONVERT)

        # initialise accelerations as zero, so that we have the objects in QUIP
        _quippy.f90wrap_atoms_add_property_real_2da(this=self._quip_atoms._handle, name='acc',
                                                    value=np.zeros((len(atoms), 3)))

        self._ds = DynamicalSystem(self._quip_atoms)

        # checking initial temperature and velocities
        if initialtemperature is not None:
            if np.max(np.abs(self._ds.atoms.velo)) > 1e-3:
                msg = 'initialtemperature given but Atoms already has non-zero velocities!'
                raise RuntimeError(msg)
            self._ds.rescale_velo(initialtemperature)

        # setting the time
        if 'time' in self.ase_atoms.info:
            self.set_time(self.ase_atoms.info['time'])  # from ASE units to fs

        self.observers = []
        self.set_timestep(timestep)

        if trajectory is not None:
            # fixme: this is not been worked on yet
            # # We also want to save the positions of all atoms after every 100th time step. --- ASE DOCS
            # traj = Trajectory('moldyn3.traj', 'w', atoms)
            # dyn.attach(traj.write, interval=50)

            # if isinstance(trajectory, str):
            #     trajectory = AtomsWriter(trajectory)
            # self.attach(trajectory, trajectoryinterval, self._ds.atoms)
            raise NotImplementedError

        # fixme: this is not been worked on yet
        self.loglabel = loglabel
        if logfile is not None:
            # if isinstance(logfile, str):
            #     logfile = InOutput(logfile, OUTPUT)
            # self.attach(self.print_status, loginterval, logfile)
            raise NotImplementedError

        # fixme: this is not been worked on yet
        self._calc_virial = False
        self._virial = np.zeros((3, 3))

    def get_time(self):
        return float(self._ds.t * fs)

    def converged(self):
        """ MD is 'converged' when number of maximum steps is reached. """
        return self.nsteps >= self.max_steps

    def get_timestep(self):
        return self._dt * fs

    def set_timestep(self, timestep):
        self._dt = timestep / fs

    timestep = property(get_timestep, set_timestep, doc="Set timestep, in ASE time units")

    def get_number_of_steps(self):
        return int(self._ds.nsteps)

    nsteps = property(get_number_of_steps)

    def insert_observer(self, function, position=0, interval=1,
                        *args, **kwargs):
        """Insert an observer."""
        if not callable(function):
            function = function.write
        self.observers.insert(position, (function, interval, args, kwargs))

    def attach(self, function, interval=1, *args, **kwargs):
        """Attach callback function.

        At every *interval* steps, call *function* with arguments
        *args* and keyword arguments *kwargs*."""

        if not hasattr(function, '__call__'):
            function = function.write
        self.observers.append((function, interval, args, kwargs))

    def call_observers(self):
        for function, interval, args, kwargs in self.observers:
            if self._ds.nsteps % interval == 0:
                function(*args, **kwargs)

    def step(self, forces):
        """
        Advance dynamics by one time-step.

        Returns forces at the end of the time-step, suitable for passing to next step()
        """
        # assert (self._ds.atoms.is_same_fortran_object(self._quip_atoms))

        # on entry we have r(t), v(t), p(t), f(t) in atoms and ds.atoms

        # keep a copy of r(t)
        r_of_t = self.ase_atoms.get_positions()

        # set current accelerations a(t) using incoming f(t)
        # only necessary for first step since advance_verlet2() does it
        # (and may add thermostat forces to acceleration, which we don't
        # want to overwrite)

        if self._ds.nsteps == 0:
            self._quip_atoms.acc[:] = forces.T / self._quip_atoms.mass

        # first half of the Velocity Verlet step for ds.atoms (which is self._quip_atoms):
        #    p(t+dt/2) = p(t) + F(t) dt/2        ->   v(t+dt/2)  = v(t) + a(t) dt/2
        #    r(t+dt)   = r(t) + p(t+dt/2)/m dt   ->   r(t+dt)    = r(t) + v(t+dt/2) dt
        # tks32: this is done with QUIP, so ase_atoms needs an update of pos and momenta
        # momenta updated directly, to avoid call on adjust forces
        self._ds.advance_verlet1(self._dt, virial=self._virial)
        self.ase_atoms.arrays['momenta'] = self.ase_atoms.get_masses()[:, np.newaxis] * \
                                           quippy.convert.velocities_quip_to_ase(self._quip_atoms.velo)
        self.ase_atoms.arrays['positions'] = self._quip_atoms.pos.T

        # now we have r(t+dt), v(t+1/2dt), p(t+1/2dt), a(t) in ds.atoms

        if self.ase_atoms.constraints:
            # fixme: there are only ASE constraints, right? I assume so here
            # keep a copy of new positions r(t+dt) so we don't lose them
            r_of_t_plus_dt = self._quip_atoms.get_positions()

            # manually revert positions of atoms to r(t)
            # NB: do not call set_positions() as this applies constraints
            # tks32: set ase_atoms only, for
            self.ase_atoms.arrays['positions'] = r_of_t

            # If there are any ASE constraints, we need to call
            # set_positions(r(t+dt)) and set_momenta(p(t+dt/2)) to invoke
            # constraints' adjust_positions() and adjust_forces() routines

            # set_positions() calls adjust_positions(), which
            # will recieve correct old and new positions
            self.ase_atoms.set_positions(r_of_t_plus_dt)
            self._quip_atoms.pos[:] = r_of_t_plus_dt.T.copy()
            self._quip_atoms.calc_dists()  # uodates dstance tables in quip, call it on mevement of atoms

            # set_momenta() calls adjust_forces() with only the new momenta
            self.ase_atoms.set_momenta(self.ase_atoms.get_momenta())

            # Update velocities v(t+dt/2) for changes in momenta made by the constraints
            self._quip_atoms.velo[:] = quippy.convert.velocities_ase_to_quip(self.ase_atoms.get_velocities())

        # Now we have r(t+dt), v(t+dt/2), p(t+dt/2)

        # Calculate f(t+dt) on atoms at the new positions r(t+dt)
        # This is invoked on the original atoms object, so that any
        # constraints' adjust_forces() routines will be called to
        # modify the forces
        forces = self.ase_atoms.get_forces()

        # compute the virial if necessary, i.e. if we have a barostat or are doing NPT
        if self._calc_virial:
            stress = self.ase_atoms.get_stress()
            virial = -stress * self.ase_atoms.get_volume()
            self._virial = np.zeros((3, 3))
            if stress.shape == (3, 3):
                self._virial[:, :] = virial
            else:
                self._virial[0, 0] = virial[0]
                self._virial[1, 1] = virial[1]
                self._virial[2, 2] = virial[2]
                self._virial[1, 2] = self._virial[2, 1] = virial[3]
                self._virial[0, 2] = self._virial[2, 0] = virial[4]
                self._virial[0, 1] = self._virial[1, 0] = virial[5]
            # print 'Computed virial:', self._virial

        # Second half of the Velocity Verlet step
        #   p(t+dt) = p(t+dt/2) + F(t+dt)/2    ->    v(t+dt) = v(t+dt/2) + a(t+dt)/2
        self._ds.advance_verlet2(self._dt, forces.T, virial=self._virial)
        self.ase_atoms.arrays['momenta'] = self.ase_atoms.get_masses()[:, np.newaxis] * \
                                           quippy.convert.velocities_quip_to_ase(self._quip_atoms.velo)

        if self.ase_atoms.constraints:
            # Update momenta, honouring constraints
            self.ase_atoms.set_momenta(self.ase_atoms.get_momenta())
            self._quip_atoms.velo[:] = quippy.convert.velocities_ase_to_quip(self.ase_atoms.get_velocities())

        # Now we have r(t+dt), v(t+dt), p(t+dt), a(t+dt) in atoms

        # Copy status into atoms.params
        # TODO: this would nee to go to the new observer, rigth?
        self._quip_atoms.params['time'] = self._ds.t * fs  # from fs to ASE time units
        self._quip_atoms.params['nsteps'] = self._ds.nsteps
        self._quip_atoms.params['cur_temp'] = self._ds.cur_temp
        self._quip_atoms.params['avg_temp'] = self._ds.avg_temp
        self._quip_atoms.params['dW'] = self._ds.dw
        self._quip_atoms.params['work'] = self._ds.work
        self._quip_atoms.params['Epot'] = self._ds.epot
        self._quip_atoms.params['Ekin'] = self._ds.ekin
        self._quip_atoms.params['Wkin'] = self._ds.wkin
        self._quip_atoms.params['thermostat_dW'] = self._ds.thermostat_dw
        self._quip_atoms.params['thermostat_work'] = self._ds.thermostat_work

        # return f(t+dt)
        return forces

    def run(self, steps=50):
        """
        Run dynamics forwards for `steps` steps.
        """
        f = self._quip_atoms.get_forces()
        for step in range(steps):
            f = self.step(f)
            self.call_observers()

    def print_status(self, file=None):
        self._ds.print_status(self.loglabel, file=file)

    def set_time(self, time):
        self._ds.t = time / fs

    time = property(get_time, set_time, doc="Time in ASE units (Note: NOT the same as femtoseconds)")

    def get_number_of_degrees_of_freedom(self):
        return self._ds.ndof

    number_of_degrees_of_freedom = property(get_number_of_degrees_of_freedom,
                                            doc="Get number of degrees of freedom in system, including QUIP constraints")

    def get_number_of_constraints(self):
        return self._ds.nconstraints

    number_of_constraints = property(get_number_of_constraints,
                                     doc="Get number of constraints")

    def get_number_of_restraints(self):
        return self._ds.nrestraints

    number_of_restraints = property(get_number_of_restraints,
                                    doc="Get number of restraints")

    def get_number_of_rigid_bodies(self):
        return self._ds.nrigid

    number_of_rigid_bodies = property(get_number_of_rigid_bodies,
                                      doc="Get number of rigid_bodies")

    def get_temperature(self):
        return self._ds.cur_temp

    def set_temperature(self, temperature):
        """
        Randomise velocities to a target temperature (given in K)
        """
        self._ds.rescale_velo(temperature)
        self.ase_atoms.arrays['momenta'] = self.ase_atoms.get_masses()[:, np.newaxis] * \
                                           quippy.convert.velocities_quip_to_ase(self._quip_atoms.velo)

    temperature = property(get_temperature, set_temperature,
                           doc="Get or set the current temperature")

    def get_average_temperature(self):
        return self._ds.avg_temp

    average_temperature = property(get_average_temperature,
                                   doc="Average temperature")

    def get_averaging_time(self):
        return self._ds.avg_time

    def set_averaging_time(self, avg_time):
        self._ds.avg_time = avg_time

    averaging_time = property(get_averaging_time, set_averaging_time,
                              doc="Averaging time used for average temperature and kinetic energy")

    def get_damping(self):
        if not self._ds.is_damping_enabled():
            return None
        return self._ds.get_damping_time()

    def set_damping(self, damp_time):
        if damp_time is None:
            self._ds.disable_damping()
        else:
            self._ds.enable_damping(damp_time)

    damping = property(get_damping, set_damping, doc=
    """Get or set the damping time constant in fs. Set to
                       None to disable damping.  By default damping applies
                       to all atoms, but can be selectively enabled with the
                       'damp_mask' property.""")

    def get_number_of_thermostats(self):
        """
        Return the number of active thermostats

        Excludes thermostat with index zero, which is used for damping.
        """
        n = self._ds.n_thermostat()
        return int(n - 1)

    def add_thermostat(self, type, T, gamma=None,
                       Q=None, tau=None, tau_cell=None,
                       p=None, NHL_tau=None, NHL_mu=None,
                       massive=None):
        """
        Add a new thermostat to this Dynamics, and return its index.

        By default, all thermostats apply to the whole system. The
        'thermostat_region' property can be used to control which
        thermostats affect which atoms. It should be set to the index
        of the required thermostat for each atom.

        :param type: string or integer giving thermostat type. The following types are supported::
           ``THERMOSTAT_NONE``, ``THERMOSTAT_LANGEVIN``, ``THERMOSTAT_NOSE_HOOVER``, ``THERMOSTAT_NOSE_HOOVER_LANGEVIN``,
           ``THERMOSTAT_LANGEVIN_NPT``, ``THERMOSTAT_LANGEVIN_PR``, ``THERMOSTAT_NPH_ANDERSEN``, ``THERMOSTAT_NPH_PR``,
           ``THERMOSTAT_LANGEVIN_OU``, ``THERMOSTAT_LANGEVIN_NPT_NB``, ``THERMOSTAT_ALL_PURPOSE``.

        :param T:  target temperature for this thermostat, in K.

        :param gamma: decay constant, in units  of 1/fs

        :param tau:  time constant, in units of fs. tau == 1/gamma, and
                     only one of the two should be specified.

        :param tau_cell: time constant for the cell degrees of freedom
                         (for variable cell thermostats only).

        :param p:   target pressure (for variable cell thermostats)

        :param NHL_tau:  time constant for Nose-Hoover-Langevein thermostats

        :param NHL_mu:   thermostat mass Nose-Hoover-Langevin thermostats

        :param massive:  set to True to enable massive Nose-Hoover-Langevin
        """

        type = _thermostat_types.get(type, type)
        if type in (THERMOSTAT_LANGEVIN_NPT, THERMOSTAT_LANGEVIN_NPT_NB):
            self._calc_virial = True
        new_index = self._ds.n_thermostat()
        region_i = np.zeros(0, dtype=np.int32)  # fixme: will this work instead of farray?

        self._ds.add_thermostat(type, T, gamma,
                                Q, tau, tau_cell,
                                p, NHL_tau,
                                NHL_mu, massive,
                                region_i=region_i)
        assert (new_index == int(region_i))
        return new_index

    def update_thermostat(self, T=None, p=None, index=1):
        """
        Change the target temperature `T` or presure `p` for thermostat `index`.
        """
        self._ds.update_thermostat(T, p, index)

    def remove_thermostat(self, index):
        """
        Remove thermostat number `index`.

        `index` should be in range 1 <= index <= dyn.get_number_of_thermostats()
        """
        n = self.get_number_of_thermostats()
        if index < 1 or index > n:
            raise ValueError('index outside range 1 <= index <= (n_thermostat=%d)' % n)
        self._ds.remove_thermostat(index)

    def print_thermostats(self):
        """
        Print a table of the current thermostats in this system
        """
        self._ds.print_thermostats()

    def get_thermostat_temperatures(self):
        """
        Return the current temperatures of all thermostated regions
        """
        temperatures = fzeros(self.get_number_of_thermostats() + 1)
        if self.damping:
            return np.array(temperatures)
        else:
            return np.array(temperatures)[1:]

    def set_barostat(self, type, p_ext, hydrostatic_strain, diagonal_strain,
                     finite_strain_formulation, tau_epsilon, W_epsilon=None,
                     T=None, W_epsilon_factor=None):
        """
        Add a barostat to this dynamical system

        Langevin NPT barostat based on fluctuation/dissipation, as
        originally implemented for hydrostatic strain using
        EOM from Quigley, D. and Probert, M.I.J., J. Chem. Phys.
        120 11432, subsequently rediscretized by Noam Bernstein
        based on ideas from and discussion with B. Leimkuhler.

        Modified for arbitrary strain based on (but not exactly using)
        formulation in E. Tadmor and R. Miller Modeling Materials:
        Continuum, Atomistic and Multiscale Techniques
        (Cambridge University Press, 2011). Chap 9.5.  Their
        definition of stress (factors of F, F^T, J, and F^-T) for finite
        strain is optionally used, but certain terms in EOM are excluded,
        and their discretization is entirely ignored.

        :param type:   The type of barostat to be used. Currently only
                       ``BAROSTAT_HOOVER_LANGEVIN`` is supported.

        :param p_ext:  External pressure in QUIP units eV/A^3

        :param hystrostatic_strain:  Set to True to constrain a
                                     hydrostatic strain tensor

        :param diagonal_strain:  Set to True to constrain a
                                 diagonal strain tensor

        :param finite_strain_formulation:  If True, use finite instead of
                      infinitessimal strain formulation

        :param tau_epsilon:  time constant for Langevin part (in fs)

        :param W_epsilon:    fictious cell mass used in Langevin NPT

        :param T:            target temperature for Langevin part

        :param W_epsilon_factor:  Scaling factor for fictious cell mass
        """

        self._calc_virial = True
        type = _barostat_types.get(type, type)
        self._ds.set_barostat(type, p_ext, hydrostatic_strain, diagonal_strain,
                              finite_strain_formulation, tau_epsilon, W_epsilon,
                              T, W_epsilon_factor)

    def update_barostat(self, p, T):
        """
        Update target pressure `p` or temperature `T` for NPT barostat
        """
        self._ds.update_barostat(self, p, T)

    def get_state(self):
        saved_ds = DynamicalSystem(self._quip_atoms)
        saved_ds.save_state(self._ds)
        return saved_ds

    def set_state(self, saved_state):
        self.restore_state(saved_state)

    state = property(get_state, set_state,
                     doc="""Save or restore current state of this dynamical system""")

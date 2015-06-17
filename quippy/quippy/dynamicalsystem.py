# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
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
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import warnings
from math import sqrt

import numpy as np

from quippy import _dynamicalsystem
from quippy._dynamicalsystem import *
from quippy.units import MASSCONVERT
from quippy.atoms import Atoms
from quippy.farray import farray
from quippy.io import AtomsWriter
from quippy.system import InOutput, OUTPUT

# time conversion units from ase.units
try:
    from ase.units import fs
except ImportError:
    _e = 1.60217733e-19          # elementary charge
    _amu = 1.6605402e-27         # atomic mass unit, kg
    second = 1e10 * sqrt(_e / _amu)
    fs = 1e-15 * second

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
    'BAROSTAT_NONE' : 0,
    'BAROSTAT_HOOVER_LANGEVIN': 1,
    }

__doc__ = _dynamicalsystem.__doc__
__all__ = ( _dynamicalsystem.__all__+ ['Dynamics', 'fs'] +
            _thermostat_types.keys() + _barostat_types.keys())


class DynamicalSystem(_dynamicalsystem.DynamicalSystem):

    __doc__ = _dynamicalsystem.DynamicalSystem.__doc__

    def run(self, pot, dt, n_steps, summary_interval=None, hook_interval=None, write_interval=None,
            trajectory=None, args_str=None, hook=None,
            save_interval=None):

        if hook is None and hook_interval is not None:
            raise ValueError('hook_interval not permitted when hook is not present')

        if hook is None:
            traj = []
            save_hook = lambda: traj.append(self.atoms.copy())
            _dynamicalsystem.DynamicalSystem.run(self, pot, dt, n_steps,
                                                 save_hook, hook_interval=save_interval,
						 summary_interval=summary_interval,
                                                 write_interval=write_interval,
                                                 trajectory=trajectory,
                                                 args_str=args_str)
            return traj
        else:
            _dynamicalsystem.DynamicalSystem.run(self, pot, dt, n_steps, hook, hook_interval=hook_interval, 
						 summary_interval=summary_interval, write_interval=write_interval,
                                                 trajectory=trajectory, args_str=args_str)

    run.__doc__ = _dynamicalsystem.DynamicalSystem.run.__doc__


from quippy import FortranDerivedTypes
FortranDerivedTypes['type(dynamicalsystem)'] = DynamicalSystem


class Dynamics(object):
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
       
        # we will do the calculation in place, to minimise number of copies,
        # unless atoms is not a quippy Atoms
        if not isinstance(atoms, Atoms):
            warnings.warn('Dynamics atoms is not quippy.Atoms instance, copy forced!')
            atoms = Atoms(atoms)
        self.atoms = atoms

        if self.atoms.has('masses'):
            if self.atoms.has_property('mass'):
                if np.max(np.abs(self.atoms.mass/MASSCONVERT - self.atoms.get_masses())) > 1e-3:
                    raise RuntimeError('Dynamics confused as atoms has inconsistent "mass" and "masses" arrays')
            else:
                self.atoms.add_property('mass', self.atoms.get_masses()*MASSCONVERT)
        else:
            if self.atoms.has_property('mass'):
                self.atoms.set_masses(self.atoms.mass/MASSCONVERT)
            else:
                self.atoms.set_masses('defaults')
        
        if self.atoms.has('momenta'):
            if self.atoms.has_property('velo'):
                if np.max(np.abs(self.atoms.velo*sqrt(MASSCONVERT) - self.atoms.get_velocities().T)) > 1e-3:
                    raise RuntimeError('Dynamics confused as atoms has inconsistent "velo" and "momenta" arrays')
            else:
                self.atoms.add_property('velo', (self.atoms.get_velocities()/sqrt(MASSCONVERT)).T)
        else:
            if self.atoms.has_property('velo'):
                self.atoms.set_velocities(self.atoms.velo.T*sqrt(MASSCONVERT))
            else:
                # start with zero momenta
                self.atoms.set_momenta(np.zeros_like(self.atoms.positions))
                self.atoms.add_property('velo', 0., n_cols=3)            

        self._ds = DynamicalSystem(self.atoms)

        if initialtemperature is not None:
            if np.max(np.abs(self._ds.atoms.velo)) > 1e-3:
                msg = ('initialtemperature given but Atoms already '+
                       'has non-zero velocities!')
                raise RuntimeError(msg)
            self._ds.rescale_velo(initialtemperature)

        # now self._ds.atoms is either a Fortran shallowcopy of atoms,
        # or a copy if input atoms was not an instance of quippy.Atoms
        
        if 'time' in atoms.info:
            self.set_time(atoms.info['time']) # from ASE units to fs

        self.observers = []
        self.set_timestep(timestep)

        if trajectory is not None:
            if isinstance(trajectory, basestring):
                trajectory = AtomsWriter(trajectory)
            self.attach(trajectory, trajectoryinterval, self._ds.atoms)

        self.loglabel = loglabel
        if logfile is not None:
            if isinstance(logfile, basestring):
                logfile = InOutput(logfile, OUTPUT)
            self.attach(self.print_status, loginterval, logfile)

        self._calc_virial = False
        self._virial = np.zeros((3,3))


    def get_timestep(self):
        return self._dt*fs

    def set_timestep(self, timestep):
        self._dt = timestep/fs

    timestep = property(get_timestep, set_timestep,
                        doc="Set timestep, in ASE time units")


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
        assert(self._ds.atoms.is_same_fortran_object(self.atoms))

        # on entry we have r(t), v(t), p(t), f(t) in atoms and ds.atoms

        # keep a copy of r(t)
        r_of_t = self.atoms.get_positions()

        # set current accelerations a(t) using incoming f(t)
        # only necessary for first step since advance_verlet2() does it
        # (and may add thermostat forces to acceleration, which we don't
        # want to overwrite)

        if self._ds.nsteps == 0:
            self.atoms.acc[...] = (forces.T)/self.atoms.mass

        # first half of the Velocity Verlet step for ds.atoms:
        #    p(t+dt/2) = p(t) + F(t) dt/2        ->   v(t+dt/2)  = v(t) + a(t) dt/2
        #    r(t+dt)   = r(t) + p(t+dt/2)/m dt   ->   r(t+dt)    = r(t) + v(t+dt/2) dt
        self._ds.advance_verlet1(self._dt, virial=self._virial)
        self.atoms.arrays['momenta'][...] = (self.atoms.velo*sqrt(MASSCONVERT)*self.atoms.mass/MASSCONVERT).T

        # now we have r(t+dt), v(t+1/2dt), p(t+1/2dt), a(t) in ds.atoms

        if self.atoms.constraints:
            # keep a copy of new positions r(t+dt) so we don't lose them
            r_of_t_plus_dt = self.atoms.get_positions()

            # manually revert positions of atoms to r(t)
            # NB: do not call set_positions() as this applies constraints
            self.atoms.arrays['positions'][...] = r_of_t

            # If there are any ASE constraints, we need to call
            # set_positions(r(t+dt)) and set_momenta(p(t+dt/2)) to invoke
            # constraints' adjust_positions() and adjust_forces() routines

            # set_positions() calls adjust_positions(), which
            # will recieve correct old and new positions
            self.atoms.set_positions(r_of_t_plus_dt)
            self.atoms.calc_dists()

            # set_momenta() calls adjust_forces() with only the new momenta
            self.atoms.set_momenta(self.atoms.get_momenta())

            # Update velocities v(t+dt/2) for changes in momenta made by the constraints
            self.atoms.velo[...] = (self.atoms.arrays['momenta'].T/self.atoms.mass*MASSCONVERT)/sqrt(MASSCONVERT)

        # Now we have r(t+dt), v(t+dt/2), p(t+dt/2)

        # Calculate f(t+dt) on atoms at the new positions r(t+dt)
        # This is invoked on the original atoms object, so that any
        # constraints' adjust_forces() routines will be called to
        # modify the forces
        forces = self.atoms.get_forces()

        # compute the virial if necessary, i.e. if we have a barostat or are doing NPT
        if self._calc_virial:
            stress = self.atoms.get_stress()
            virial = -stress*self.atoms.get_volume()
            self._virial = np.zeros((3,3))
            if stress.shape == (3,3):
                self._virial[:,:] = virial
            else:
                self._virial[0,0] = virial[0]
                self._virial[1,1] = virial[1]
                self._virial[2,2] = virial[2]
                self._virial[1,2] = self._virial[2,1] = virial[3]
                self._virial[0,2] = self._virial[2,0] = virial[4]
                self._virial[0,1] = self._virial[1,0] = virial[5]
            print 'Computed virial:', self._virial

        # Second half of the Velocity Verlet step
        #   p(t+dt) = p(t+dt/2) + F(t+dt)/2    ->    v(t+dt) = v(t+dt/2) + a(t+dt)/2
        self._ds.advance_verlet2(self._dt, forces.T, virial=self._virial)
        self.atoms.arrays['momenta'][...] = (self.atoms.velo*sqrt(MASSCONVERT)*self.atoms.mass/MASSCONVERT).T


        if self.atoms.constraints:
            # Update momenta, honouring constraints
            self.atoms.set_momenta(self.atoms.get_momenta())
            self.atoms.velo[...] = (self.atoms.arrays['momenta'].T/self.atoms.mass*MASSCONVERT)/sqrt(MASSCONVERT)

        # Now we have r(t+dt), v(t+dt), p(t+dt), a(t+dt) in atoms

        # Copy status into atoms.params
        self.atoms.params['time'] = self._ds.t*fs # from fs to ASE time units
        self.atoms.params['nsteps'] = self._ds.nsteps
        self.atoms.params['cur_temp'] = self._ds.cur_temp
        self.atoms.params['avg_temp'] = self._ds.avg_temp
        self.atoms.params['dW'] = self._ds.dw
        self.atoms.params['work'] = self._ds.work
        self.atoms.params['Epot'] = self._ds.epot
        self.atoms.params['Ekin'] = self._ds.ekin
        self.atoms.params['Wkin'] = self._ds.wkin
        self.atoms.params['thermostat_dW'] = self._ds.thermostat_dw
        self.atoms.params['thermostat_work'] = self._ds.thermostat_work
 
        # return f(t+dt)
        return forces


    def run(self, steps=50):
        """
        Run dynamics forwards for `steps` steps. 
        """
        f = self.atoms.get_forces()
        for step in xrange(steps):
            f = self.step(f)
            self.call_observers()            


    def print_status(self, file=None):
        self._ds.print_status(self.loglabel, file=file)
        

    def get_time(self):
        return float(self._ds.t*fs)

    def set_time(self, time):
        self._ds.t = time/fs

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
        self.atoms.set_momenta((self.atoms.velo*sqrt(MASSCONVERT)*self.atoms.mass/MASSCONVERT).T)

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
        return int(n-1)


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
        region_i = farray(0, dtype=np.int32)

        self._ds.add_thermostat(type, T, gamma,
                                Q, tau, tau_cell,
                                p, NHL_tau,
                                NHL_mu, massive,
                                region_i=region_i)
        assert(new_index == int(region_i))
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
        temperatures = fzeros(self.get_number_of_thermostats()+1)
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
        saved_ds = DynamicalSystem(self.atoms)
        saved_ds.save_state(self._ds)
        return saved_ds

    def set_state(self, saved_state):
        self.restore_state(saved_state)

    state = property(get_state, set_state,
                     doc="""Save or restore current state of this dynamical system""")



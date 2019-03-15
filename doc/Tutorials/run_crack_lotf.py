"""
run_crack_lotf.py
Script to run LOTF molecular dynamics for a crack slab,
incrementing the load in small steps until fracture starts.
James Kermode <james.kermode@kcl.ac.uk>
February 2013
"""
import numpy as np

from ase.constraints import FixAtoms
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
import ase.units as units

from quippy import set_fortran_indexing
from quippy.atoms import Atoms
from quippy.potential import Potential
from quippy.io import AtomsWriter

from quippy.crack import (get_strain,
                          get_energy_release_rate,
                          ConstantStrainRate,
                          find_crack_tip_stress_field)

# additional requirements for the QM/MM simulation:
from quippy.potential import ForceMixingPotential
from quippy.lotf import LOTFDynamics, update_hysteretic_qm_region


# ******* Start of parameters ***********

input_file = 'crack.xyz'         # File from which to read crack slab structure
sim_T = 300.0*units.kB           # Simulation temperature
nsteps = 1000000                   # Total number of timesteps to run for
timestep = 1.0*units.fs          # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang      # Amount by which potential cutoff is increased
                                 # for neighbour calculations
tip_move_tol = 10.0              # Distance tip has to move before crack
                                 # is taken to be running
strain_rate = 1e-5*(1/units.fs)  # Strain rate
traj_file = 'traj.xyz'           # Trajectory output file (NetCDF or XYZ format)
traj_interval = 10               # Number of time steps between
                                 # writing output frames
param_file = 'params.xml'        # Filename of XML file containing
                                 # potential parameters
mm_init_args = 'IP SW'           # Initialisation arguments for
                                 # classical potential

# additional parameters for the QM/MM simulation:
qm_init_args = 'TB DFTB'         # Initialisation arguments for QM potential
qm_inner_radius = 6.0*units.Ang  # Inner hysteretic radius for QM region
qm_outer_radius = 8.0*units.Ang  # Inner hysteretic radius for QM region
extrapolate_steps = 10           # Number of steps for predictor-corrector
                                 # interpolation and extrapolation

# ******* End of parameters *************

set_fortran_indexing(False)

# ********** Read input file ************

print 'Loading atoms from file %s' % input_file
atoms = Atoms(input_file)

orig_height = atoms.info['OrigHeight']
orig_crack_pos = atoms.info['CrackPos'].copy()

# ***** Setup constraints *******

top = atoms.positions[:, 1].max()
bottom = atoms.positions[:, 1].min()
left = atoms.positions[:, 0].min()
right = atoms.positions[:, 0].max()

# fix atoms in the top and bottom rows
fixed_mask = ((abs(atoms.positions[:, 1] - top) < 1.0) |
              (abs(atoms.positions[:, 1] - bottom) < 1.0))
fix_atoms = FixAtoms(mask=fixed_mask)
print('Fixed %d atoms\n' % fixed_mask.sum())
atoms.set_constraint([fix_atoms])

# Increase epsilon_yy applied to all atoms at constant strain rate
strain_atoms = ConstantStrainRate(orig_height, strain_rate*timestep)

# ******* Set up potentials and calculators ********

mm_pot = Potential(mm_init_args,
                   param_filename=param_file,
                   cutoff_skin=cutoff_skin)

# Density functional tight binding (DFTB) potential
qm_pot = Potential(qm_init_args,
                   param_filename=param_file)

# Construct the QM/MM potential, which mixes QM and MM forces.
# The qm_args_str parameters control how the QM calculation is carried out:
# we use a single cluster, periodic in the z direction and terminated
# with hydrogen atoms. The positions of the outer layer of buffer atoms
# are not randomised.
qmmm_pot = ForceMixingPotential(pot1=mm_pot,
                                pot2=qm_pot,
                                qm_args_str='single_cluster cluster_periodic_z carve_cluster '+
                                            'terminate cluster_hopping=F randomise_buffer=F',
                                fit_hops=4,
                                lotf_spring_hops=3,
                                hysteretic_buffer=True,
                                hysteretic_buffer_inner_radius=7.0,
                                hysteretic_buffer_outer_radius=9.0,
                                cluster_hopping_nneighb_only=False,
                                min_images_only=True)

# Use the force mixing potential as the Atoms' calculator
atoms.set_calculator(qmmm_pot)
qmmm_pot.atoms = atoms

# *** Set up the initial QM region ****

qm_list = update_hysteretic_qm_region(atoms, [], orig_crack_pos,
                                      qm_inner_radius, qm_outer_radius)

# ********* Setup and run MD ***********

# Set the initial temperature to 2*simT: it will then equilibriate to
# simT, by the virial theorem
MaxwellBoltzmannDistribution(atoms, 2.0*sim_T)

# Initialise the dynamical system
dynamics = LOTFDynamics(atoms, timestep, extrapolate_steps)

# array to store time averaged stress field
avg_sigma = np.zeros((len(atoms), 3, 3))

# Print some information every time step
def printstatus():
    if dynamics.nsteps == 1:
        print """
State      Time/fs    Temp/K     Strain      G/(J/m^2)  CrackPos/A D(CrackPos)/A
---------------------------------------------------------------------------------"""

    log_format = ('%(label)-4s%(time)12.1f%(temperature)12.6f'+
                  '%(strain)12.5f%(G)12.4f%(crack_pos_x)12.2f    (%(d_crack_pos_x)+5.2f)')

    atoms.info['label'] = dynamics.state_label  # Label for the status line
    atoms.info['time'] = dynamics.get_time()/units.fs
    atoms.info['temperature'] = (atoms.get_kinetic_energy() /
                                 (1.5*units.kB*len(atoms)))
    atoms.info['strain'] = get_strain(atoms)
    atoms.info['G'] = get_energy_release_rate(atoms)/(units.J/units.m**2)

    crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot,
                                            avg_sigma=avg_sigma)
    atoms.info['crack_pos_x'] = crack_pos[0]
    atoms.info['d_crack_pos_x'] = crack_pos[0] - orig_crack_pos[0]

    print log_format % atoms.info

dynamics.attach(printstatus)

def atom_straining(atoms):
    crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot,
                                            avg_sigma=avg_sigma)
    # keep straining until the crack tip has advanced to tip_move_tol
    if not atoms.info['is_cracked'] and (crack_pos[0] - orig_crack_pos[0]) < tip_move_tol:
      strain_atoms.apply_strain(atoms)
    elif not atoms.info['is_cracked']:
      atoms.info['is_cracked'] = True

dynamics.attach(atom_straining, 1, dynamics.atoms)

# Function to update the QM region at the beginning of each extrapolation cycle
def update_qm_region(atoms):
   crack_pos = find_crack_tip_stress_field(atoms, calc=mm_pot,
                                           avg_sigma=avg_sigma)
   qm_list = qmmm_pot.get_qm_atoms(atoms)
   qm_list = update_hysteretic_qm_region(atoms, qm_list, crack_pos,
                                         qm_inner_radius, qm_outer_radius)
   qmmm_pot.set_qm_atoms(qm_list, atoms)
   #assert (atoms.hybrid == 1).sum() == len(qm_list)

dynamics.set_qm_update_func(update_qm_region)


# Save frames to the trajectory every `traj_interval` time steps
# but only when interpolating
trajectory = AtomsWriter(traj_file)

def traj_writer(dynamics):
   if dynamics.state == LOTFDynamics.Interpolation:
      # copy time-averaged stress into Atoms so it will be written to trajectory file
      dynamics.atoms.set_array('avg_sigma', avg_sigma.reshape((len(atoms), 9)))
      trajectory.write(dynamics.atoms)

dynamics.attach(traj_writer, traj_interval, dynamics)

# Start running!
dynamics.run(nsteps)

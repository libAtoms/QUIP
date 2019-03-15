"""
run_crack_classical_1.py

Script to run classical molecular dynamics for a crack slab,
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

# ******* Start of parameters ***********

input_file = 'crack.xyz'         # File from which to read crack slab structure
sim_T = 300.0*units.kB           # Simulation temperature
nsteps = 10000                   # Total number of timesteps to run for
timestep = 1.0*units.fs          # Timestep (NB: time base units are not fs!)
cutoff_skin = 2.0*units.Ang      # Amount by which potential cutoff is increased
                                 # for neighbour calculations
tip_move_tol = 10.0              # Distance tip has to move before crack
                                 # is taken to be running
strain_rate = 1e-5*(1/units.fs)  # Strain rate
traj_file = 'traj.nc'            # Trajectory output file (NetCDF format)
traj_interval = 10               # Number of time steps between
                                 # writing output frames
param_file = 'params.xml'        # Filename of XML file containing
                                 # potential parameters
mm_init_args = 'IP SW'           # Initialisation arguments for
                                 # classical potential

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

# Request Potential to compute per-atom stresses whenever we
# compute forces, to save time when locating the crack tip
mm_pot.set_default_quantities(['stresses'])

atoms.set_calculator(mm_pot)


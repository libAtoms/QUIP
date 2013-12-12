"""
make_crack.py

Script to generate a crack slab, and apply initial strain ramp

James Kermode <james.kermode@kcl.ac.uk>
January 2013
"""
from ase.structure import bulk
from ase.lattice.cubic import Diamond
from ase.constraints import FixAtoms
import ase.units as units

from quippy import set_fortran_indexing
from quippy.atoms import Atoms
from quippy.potential import Potential, Minim
from quippy.elasticity import youngs_modulus, poisson_ratio
from quippy.io import write
from quippy.crack import (print_crack_system,
                          G_to_strain,
                          thin_strip_displacement_y,
                          find_crack_tip_stress_field)

# ******* Start of parameters ***********

# There are three possible crack systems, choose one and uncomment it

# System 1. (111)[0-11]
crack_direction = (-2, 1, 1)      # Miller index of x-axis
cleavage_plane = (1, 1, 1)        # Miller index of y-axis
crack_front = (0, 1, -1)          # Miller index of z-axis

# # System 2. (110)[001]
# crack_direction = (1,-1,0)
# cleavage_plane = (1,1,0)
# crack_front = (0,0,1)

# # System 3. (110)[1-10]
# crack_direction = (0,0,-1)
# cleavage_plane = (1,1,0)
# crack_front = (1,-1,0)

width = 200.0*units.Ang              # Width of crack slab
height = 100.0*units.Ang             # Height of crack slab
vacuum = 100.0*units.Ang             # Amount of vacuum around slab
crack_seed_length = 40.0*units.Ang   # Length of seed crack
strain_ramp_length = 30.0*units.Ang  # Distance over which strain is ramped up
initial_G = 5.0*(units.J/units.m**2) # Initial energy flow to crack tip

relax_fmax = 0.025*units.eV/units.Ang  # Maximum force criteria for relaxation

param_file = 'params.xml'            # XML file containing interatomic potential parameters
mm_init_args = 'IP SW'               # Initialisation arguments for the classical potential

output_file = 'crack.xyz'            # File to which structure will be written

# ******* End of parameters *************

set_fortran_indexing(False)

# ********** Build unit cell ************

# 8-atom diamond cubic unit cell for silicon, with guess at lattice
# constant of 5.44 A
si_bulk = bulk('Si', 'diamond', a=5.44, cubic=True)


# ********** Setup potential ************

# Stillinger-Weber (SW) classical interatomic potential, from QUIP
mm_pot = Potential(mm_init_args,
                   param_filename=param_file)


# ***** Find eqm. lattice constant ******

# find the equilibrium lattice constant by minimising atoms wrt virial
# tensor given by SW pot (possibly replace this with parabola fit in
# another script and hardcoded a0 here)
si_bulk.set_calculator(mm_pot)

print('Minimising bulk unit cell...')
minim = Minim(si_bulk, relax_positions=True, relax_cell=True)
minim.run(fmax=1e-4)
    
a0 = si_bulk.cell[0, 0]
print('Lattice constant %.3f A\n' % a0)

# make a new bulk cell with correct a0 (so that off-diagonal lattice values are exactly zero)
si_bulk = bulk('Si', 'diamond', a=a0, cubic=True)
si_bulk.set_calculator(mm_pot)

# ******* Find elastic constants *******

# Get 6x6 matrix of elastic constants C_ij
c = mm_pot.get_elastic_constants(si_bulk)
print('Elastic constants (GPa):')
print((c / units.GPa).round(0))
print('')

E = youngs_modulus(c, cleavage_plane)
print('Young\'s modulus %.1f GPa' % (E / units.GPa))
nu = poisson_ratio(c, cleavage_plane, crack_direction)
print('Poisson ratio % .3f\n' % nu)

# **** Setup crack slab unit cell ******

print_crack_system(crack_direction, cleavage_plane, crack_front)

# now, we build system aligned with requested crystallographic orientation
unit_slab = Diamond(directions=[crack_direction,
                                cleavage_plane,
                                crack_front],
                    size=(1, 1, 1),
                    symbol='Si',
                    pbc=True,
                    latticeconstant=a0)

# You could check elastic constants of this rotated system:
# should lead to same Young's modulus and Poisson ratio

print('Unit slab with %d atoms per unit cell:' % len(unit_slab))
print(unit_slab.cell)
print('')

# center vertically half way along the vertical bond between atoms 0 and 1
unit_slab.positions[:, 1] += (unit_slab.positions[1, 1] -
                              unit_slab.positions[0, 1]) / 2.0

# map positions back into unit cell
unit_slab.set_scaled_positions(unit_slab.get_scaled_positions())

# Make a surface unit cell by adding some vaccum along y
surface = unit_slab.copy()
surface.center(vacuum, axis=1)


# ********** Surface energy ************

# Calculate surface energy per unit area
surface.set_calculator(mm_pot)
E_surf = surface.get_potential_energy()
E_per_atom_bulk = si_bulk.get_potential_energy() / len(si_bulk)
area = surface.get_volume() / surface.cell[1, 1]
gamma = ((E_surf - E_per_atom_bulk * len(surface)) /
         (2.0 * area))

print('Surface energy of %s surface %.4f J/m^2\n' %
      (cleavage_plane, gamma / (units.J / units.m ** 2)))


# ***** Setup crack slab supercell *****

# Now we will build the full crack slab system,
# approximately matching requested width and height
nx = int(width / unit_slab.cell[0, 0])
ny = int(height / unit_slab.cell[1, 1])

# make sure ny is even so slab is centered on a bond
if ny % 2 == 1:
    ny += 1

# make a supercell of unit_slab
crack_slab = unit_slab * (nx, ny, 1)

# open up the cell along x and y by introducing some vaccum
crack_slab.center(vacuum, axis=0)
crack_slab.center(vacuum, axis=1)

# centre the slab on the origin
crack_slab.positions[:, 0] -= crack_slab.positions[:, 0].mean()
crack_slab.positions[:, 1] -= crack_slab.positions[:, 1].mean()

orig_width = (crack_slab.positions[:, 0].max() -
              crack_slab.positions[:, 0].min())
orig_height = (crack_slab.positions[:, 1].max() -
               crack_slab.positions[:, 1].min())

print(('Made slab with %d atoms, original width and height: %.1f x %.1f A^2' %
       (len(crack_slab), orig_width, orig_height)))

top = crack_slab.positions[:, 1].max()
bottom = crack_slab.positions[:, 1].min()
left = crack_slab.positions[:, 0].min()
right = crack_slab.positions[:, 0].max()

# fix atoms in the top and bottom rows
fixed_mask = ((abs(crack_slab.positions[:, 1] - top) < 1.0) |
              (abs(crack_slab.positions[:, 1] - bottom) < 1.0))
const = FixAtoms(mask=fixed_mask)
crack_slab.set_constraint(const)
print('Fixed %d atoms\n' % fixed_mask.sum())


# ****** Apply initial strain ramp *****

strain = G_to_strain(initial_G, E, nu, orig_height)

crack_slab.positions[:, 1] += thin_strip_displacement_y(
                                 crack_slab.positions[:, 0],
                                 crack_slab.positions[:, 1],
                                 strain,
                                 left + crack_seed_length,
                                 left + crack_seed_length + strain_ramp_length)

print('Applied initial load: strain=%.4f, G=%.2f J/m^2' %
      (strain, initial_G / (units.J / units.m**2)))


# ***** Relaxation of crack slab  *****

# optionally, relax the slab, keeping top and bottom rows fixed
print('Relaxing slab...')
crack_slab.set_calculator(mm_pot)
minim = Minim(crack_slab, relax_positions=True, relax_cell=False)
minim.run(fmax=relax_fmax)

# Find initial position of crack tip
crack_pos = find_crack_tip_stress_field(crack_slab, calc=mm_pot)
print 'Found crack tip at position %s' % crack_pos

# Save all calculated materials properties inside the Atoms object
crack_slab.info['nneightol'] = 1.3 # nearest neighbour tolerance
crack_slab.info['LatticeConstant'] = a0
crack_slab.info['C11'] = c[0, 0]
crack_slab.info['C12'] = c[0, 1]
crack_slab.info['C44'] = c[3, 3]
crack_slab.info['YoungsModulus'] = E
crack_slab.info['PoissonRatio_yx'] = nu
crack_slab.info['SurfaceEnergy'] = gamma
crack_slab.info['OrigWidth'] = orig_width
crack_slab.info['OrigHeight'] = orig_height
crack_slab.info['CrackDirection'] = crack_direction
crack_slab.info['CleavagePlane'] = cleavage_plane
crack_slab.info['CrackFront'] = crack_front
crack_slab.info['strain'] = strain
crack_slab.info['G'] = initial_G
crack_slab.info['CrackPos'] = crack_pos
crack_slab.info['is_cracked'] = False


# ******** Save output file **********

# save results in extended XYZ format, including extra properties and info
print('Writing crack slab to file %s' % output_file)
write(output_file, crack_slab)

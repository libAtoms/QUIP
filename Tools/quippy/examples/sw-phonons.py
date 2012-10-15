from phonopy import Phonopy
from phonopy.structure.atoms import Atoms as PhonopyAtoms
from phonopy.units import VaspToTHz
import numpy as np
from phonopy.phonon.band_structure import BandStructure
from pylab import clf, ylim, vlines, savefig, draw
from quippy import diamond, Potential

# Define unit cell and potential to be used for calculating forces
a = 5.44
bulk = diamond(a, 14)
pot = Potential('IP SW')
bulk.set_calculator(pot)

# Phonopy pre-process
print "------"
print "Phonon"
print "------"
# 1st arg. is the input unit cell.
# 2nd arg. is the supercell lattice relative to the unit cell.
# 'distance' is the distance of displacements.
# Default symmetry tolerance is 1e-5 in fractional coordinates.
phonon=Phonopy(bulk, [[1,0,0],[0,1,0],[0,0,1]], distance=0.01, factor=VaspToTHz)
phonon.print_displacements()
supercells = phonon.get_supercells_with_displacements()

# Do force calculations using QUIP Potential
set_of_forces = []
for cell in supercells:
    forces = pot.get_forces(cell)
    drift_force = forces.sum(axis=0)
    print "        ---------------------------------"
    print "     ", "%11.5f"*3 % tuple(drift_force)
    # Simple translational invariance
    for force in forces:
        force -= drift_force / forces.shape[0]
    set_of_forces.append(forces)

# Phonopy post-process
# 1st arg. is transformation from input unit cell to the primitive lattice
# 2nd arg. is list of the calculated forces.
phonon.set_post_process([[0,0.5,0.5],[0.5,0,0.5],[0.5,0.5,0]], set_of_forces)

# Gamma, X(1/2, 1/2,0),  W(3/4,1/4,1/2), X, k(3/4,3/8,3/8), Gamma, L(1/2,1/2,1/2)
bands = []         # Gamma to X
q_start  = np.array([0.0,0.0,0.0])
q_end    = np.array([0.5,0.5,0.0])
band = []
for i in range(51):
    band.append(q_start + (q_end - q_start) / 50 * i)
bands.append(band)

q_start  = np.array([0.5,0.5,0.0])  # X to W
q_end   =  np.array([0.75,0.25,0.5])
band = []
for i in range(51):
    band.append(q_start + (q_end - q_start) / 50 * i)
bands.append(band)

q_start  = np.array([0.75,0.25,0.5]) # W to X
q_end    = np.array([0.5,0.5,0.0])
band = []
for i in range(51):
    band.append(q_start + (q_end - q_start) / 50 * i)
bands.append(band)

q_start  = np.array([0.5,0.5,0.0])       # X to K
q_end    = np.array([0.75,0.375,0.375])
band = []
for i in range(51):
    band.append(q_start + (q_end - q_start) / 50 * i)
bands.append(band)

q_start  = np.array([0.75,0.375,0.375]) # K to Gamma
q_end = np.array([0.0, 0.0, 0.0])
band = []
for i in range(51):
    band.append(q_start + (q_end - q_start) / 50 * i)
bands.append(band)

q_start  = np.array([0.0,0.0,0.0])   # Gamma to L 
q_end    = np.array([0.5,0.5,0.5])
band = []
for i in range(51):
    band.append(q_start + (q_end - q_start) / 50 * i)
bands.append(band)

# plot the band structure
clf()
phonon.set_band_structure(bands)
phonon.plot_band_structure([r'$\Gamma$', '$X$', '$W$', '$X$', '$K$', r'$\Gamma$', '$L$'])
ylim(0,20.)

# save results to file band.yaml
bs=BandStructure(bands, phonon.dynamical_matrix, phonon.primitive, factor=VaspToTHz)
bs.write_yaml()

# add vertical  lines to plot
vlines(bs.special_point, 0, 20.)
draw()
savefig('sw-phonons.pdf')

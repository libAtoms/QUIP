#!/usr/bin/env python

from quippy import *
from quippy.crack import *
from numpy import *

stem = sys.argv[1]

a0 = Atoms(stem+'.xyz')
bulk = Atoms(stem+'_bulk.xyz')

traj = AtomsReader(sys.argv[2])

for at in traj:
    at.add_property('move_mask', a0.move_mask)
    G = crack_strain_energy_release_rate(at, bulk, stem=stem)

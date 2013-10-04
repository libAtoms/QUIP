#!/usr/bin/env python

import sys
import os
from math import pi
import numpy as np

from quippy.atoms import Atoms
from quippy.structures import void_analysis
from quippy.farray import fzeros

infile = sys.argv[1]
basename = os.path.splitext(infile)[0]

a = Atoms(infile)

grid_size = 0.25
min_void_size = 2.0

nx = int((a.pos[1,:].max() - a.pos[1,:].min())/grid_size)
ny = int((a.pos[2,:].max() - a.pos[2,:].min())/grid_size)
nz = int((a.pos[3,:].max() - a.pos[3,:].min())/grid_size)
n = nx*ny*nz

grid = fzeros((3, n))
radii = fzeros(n)

void_analysis(a, grid_size, min_void_size, grid, radii)

extent = (grid[:,-1] - grid[:,1])

grid = np.array(grid.reshape(3, nx, ny, nz))
radii = np.array(radii.reshape((nx, ny, nz)))

a.write(basename+'.cube', data=radii, extent=extent, 
        origin=[a.pos[1,:].min(), a.pos[2,:].min(), a.pos[3,:].min()])

voids = (radii > min_void_size).nonzero()

void_positions = []
void_radii = []

while len(voids[0]) > 0:
    cvoids = np.c_[voids]

    rs = radii[voids]
    biggest_void = rs.argmax()

    bx, by, bz = cvoids[biggest_void,:]
    p = grid[:, bx, by, bz]
    r = rs[biggest_void]
    v = 4./3*pi*r**3

    print '%d voids remaining' % len(rs)
    print 'biggest void %r, pos %r radius %.2f' % (cvoids[biggest_void,:], p, r)
    
    void_pos = np.array([grid[:, i, j, k] for (i,j,k) in cvoids])

    # keep only the non-overlapping voids
    void_sep = np.array([a.distance_min_image(p, vp) for vp in void_pos])
    non_overlap_voids = cvoids[void_sep >= r + rs]
    voids = (non_overlap_voids[:,0], non_overlap_voids[:,1], non_overlap_voids[:,2])

    void_positions.append(p)
    void_radii.append(r)

void_positions = np.array(void_positions)
void_radii = np.array(void_radii)

#a.add_atoms(pos=void_positions.T, z=[5]*len(void_radii))

def write_vmd_script(f, positions, radii):
    f = open(f, 'w')
    f.write('draw delete all\n')
    for (pos, radius) in zip(positions, radii):
        f.write('draw sphere {%.3f %.3f %.3f} radius %.3f\n' % tuple(list(pos) + [radius]))
    f.close()

np.savetxt(basename+'.dat', np.c_[void_radii, void_positions])
write_vmd_script(basename+'.tcl', void_positions, void_radii)


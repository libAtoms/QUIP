#!/usr/bin/env python
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

from enthought.mayavi import mlab
import quippy
import numpy as np

atom_colours_lut = np.array([list(quippy.ElementColours[k])+[1] for k in range(110)])*255

def add_balls(at, colours=None, colorbar=False, atom_scale_factor=1.0, vmin=None, vmax=None):

    r = np.array([quippy.ElementCovRad[z] for z in at.z])

    if colours is None:
        scalars = at.z
        vmin = 0
        vmax = 110
    else:
        scalars = colours
    
    pts = mlab.quiver3d(at.pos[1,:], at.pos[2,:], at.pos[3,:], r, [0]*len(r), [0]*len(r), scale_factor=atom_scale_factor,
                        scale_mode = 'vector',
                        scalars=scalars, mode='sphere', name='balls', vmin=vmin, vmax=vmax)
    
    pts.glyph.color_mode = 'color_by_scalar'
    #pts.glyph.glyph.scale_mode = 'scale_by_vector'
    pts.glyph.glyph_source.glyph_source.center = [0, 0, 0]

    if colours is None:
        pts.module_manager.scalar_lut_manager.lut.table = atom_colours_lut
    
    if colorbar: mlab.colorbar()
        
    return pts


def add_bonds(at, pts, cutoff_factor=1.2, bond_radius=0.2, bond_colour=(.55, .55, .55)):
    at.set_cutoff_factor(cutoff_factor)
    at.calc_connect()

    bonds = []
    for i in quippy.frange(at.n):
        i_is_min_image = at.connect.is_min_image[i]
        for n in at.neighbours[i]:
            # only include minimum image bonds with i < j
            j_is_min_image = at.connect.is_min_image[n.j]
            if i_is_min_image and j_is_min_image and i > n.j: continue

            # only include bonds within this periodic image
            if any(n.shift != (0, 0, 0)): continue
            bonds.append([i-1,n.j-1])

    pts.mlab_source.dataset.lines = np.array(bonds)
    return mlab.pipeline.surface(mlab.pipeline.tube(pts, tube_radius=bond_radius), color=bond_colour)
    
    

def add_cell(at, pts, origin=[0.,0.,0.], shift=[.5, .5, .5], supercell=[1, 1, 1]):
    if not at.is_orthorhombic:
        raise ValueError('Cannot draw cell outline for non-orthorhombic cells')

    na, nb, nc = supercell
    extent = np.array([0., at.lattice[1,1], 0., at.lattice[2,2], 0., at.lattice[3,3]])
    shift = np.dot(np.array(at.lattice), shift)
    shift = np.array([shift[0], shift[0], shift[1], shift[1], shift[2], shift[2]])
    origin = np.array(origin)
    origin = np.array([origin[0], origin[0], origin[1], origin[1], origin[2], origin[2]])

    print 'origin', origin
    print 'shift', shift
    
    for i in range(na):
        for j in range(nb):
            for k in range(nc):
                disp = np.array([i*at.lattice[1,1], i*at.lattice[1,1],
                                 j*at.lattice[2,2], j*at.lattice[2,2],
                                 k*at.lattice[3,3], k*at.lattice[3,3]])
                return mlab.outline(pts, extent=origin+disp+extent-shift, name='outline_%d_%d_%d' % (i,j,k))


def add_arrows(at, arrows, arrow_colour=(0, 0, 0), arrow_scale_factor=1.0, clf=True):
    pos = np.array(at.pos)
    arrows = np.array(arrows)
    return mlab.quiver3d(pos[0,:], pos[1,:], pos[2,:],
                  arrows[0,:], arrows[1,:], arrows[2,:],
                  scale_factor=arrow_scale_factor, color=arrow_colour, name='arrows')

def scalar_field(at, x=None, y=None, z=None, v=None):
    # Default to attributes read from a .cube file
    if x is None: x = at.grid_x
    if y is None: y = at.grid_y
    if z is None: z = at.grid_z
    if v is None: v = at.data

    return mlab.pipeline.scalar_field(x, y, z, v)

def draw_atoms(at, colours=None, colorbar=None, atom_scale_factor=1.0,
               cell=True, origin=[0., 0., 0.], shift=[.5, .5, .5], supercell=[1, 1, 1],
               bonds=True, cutoff_factor=1.2, bond_radius=0.2, bond_colour=(.55, .55, .55),
               arrows=None, arrow_colour=(0,0,0), arrow_scale_factor=1.0, clf=True):

    if clf: mlab.clf()

    fig = mlab.gcf()
    fig.scene.disable_render = True

    balls = add_balls(at, colours, colorbar, atom_scale_factor)

    if cell:
        add_cell(at, balls, origin, shift, supercell)
        
    if bonds:
        add_bonds(at, balls, cutoff_factor, bond_radius, bond_colour)

    if arrows:
        arrows = add_arrows(at, arrows, arrow_colour, arrow_scale_factor)

    fig.scene.disable_render = False

    if arrows:
        return (balls, arrows)
    else:
        return balls
        

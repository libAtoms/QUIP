#!/usr/bin/env python

from quippy import *
from quippy.surface import *
from quippy.structures import *
from quippy.elasticity import *
import optparse
import numpy as np

p = optparse.OptionParser(usage='%prog [options]')

p.add_option('-s', '--struct', action='store', help='Structure name (mandatory)')
p.add_option('-v', '--vol-per-atom', action='store', type='float', help='Volume per atom (default -1)', default=-1.0)
p.add_option('-V', '--vol-per-unit-cell', action='store', type='float', help='Volume per unit cell (default -1)', default=-1.0)
p.add_option('-r', '--repeat', action='store', help='Supercell repeat (default "1 1 1")', default='1 1 1')
p.add_option('-z', '--z-values', action='store', help='Atomic number values (1 per atom in primitive cell)')

p.add_option('--at-bulk', action='store', help='Input file for primitive bulk unit cell')
p.add_option('--init-args', action='store', help='Potential init_args, e.g. {IP SW}')
p.add_option('--param-file', action='store', help='XML parameter filename')

p.add_option('--relax-lattice', action='store_true', help='Relax lattice of initial structure')
p.add_option('--cij-virial', action='store_true', help='Calculate elastic constants C_ij using virial')
p.add_option('--cij-fit', action='store_true', help='Calculate elastic constants C_ij using fitting of stress vs. strain')
p.add_option('--cij-energy', action='store_true', help='Calculate elastic constants C_ij using fitting of energy vs. strain')
p.add_option('--cij-symmetry', action='store', help='Symmetry name for C_ij fitting (default monoclinic)', default='monoclinic')

p.add_option('--surface-energy', action='store_true', help='Calculate surface energy')
p.add_option('--surface', action='store', help='Generate surface with given Miller indices (e.g. (111)[11b0])')
p.add_option('--surface-vacuum', action='store', help='Amount of vacuum to add to surface cell (default 10 A)', type=float, default=10.0)
p.add_option('--surface-normal', action='store', help='Axis of surface normal: x, y or z (default is y)', default='y')
p.add_option('--relax-surface', action='store_true', help='Calculate relaxed surface energy')
p.add_option('--at-surface', action='store', help='Input file for surface cell')




opt, args = p.parse_args()

if len(args) != 0:
    p.error('No arguments are required but got %r' % args)

if opt.repeat is not None:
    try:
        opt.repeat = [int(x) for x in opt.repeat.split()]
    except ValueError:
        p.error('Cannot parse repeat argument %s' % opt.repeat)

if opt.init_args is None:
    import xml.dom.minidom
    xml = xml.dom.minidom.parse(opt.param_file)
    opt.init_args = 'xml_label=%s' % str(xml.documentElement.tagName)

try:

    if opt.at_bulk is not None:
        at = Atoms(opt.at_bulk)
    else:
        at = structure_from_file(opt.struct, opt.vol_per_atom, opt.vol_per_unit_cell, opt.repeat, opt.z_values)

    pot = Potential(opt.init_args, param_filename=opt.param_file)

    pot.print_()

    at.set_cutoff(pot.cutoff()+1.0)
    at.calc_connect()

    if opt.relax_lattice:
        relaxed = at.copy()
        verbosity_push_decrement()
        pot.minim(relaxed, 'cg', 1e-3, 100, do_pos=True, do_lat=True, use_n_minim=True)
        verbosity_pop()

        # Remove near-zero elements, and exactly restore original symmetries
        at.set_lattice(np.where(abs(relaxed.lattice) > 1e-5, relaxed.lattice, 0), True)

        print 'LATTICE'
        print 'Relaxed lattice / A'
        print at.lattice
        print

        print 'FRACTIONAL POSITIONS'
        for sp, t in zip(at.species.stripstrings(), np.dot(at.g, at.pos)):
           print '%-3s [%8.6f %8.6f %8.6f]' % (sp, t[1], t[2], t[3])

    b, v0 = pot.bulk_modulus(at, minimise_bulk=False)
    print 'BULK MODULUS B=', b, 'GPa'
    print 'EQUILIBRIUM CELL VOLUME V0=', v0, 'A**3'

    if opt.cij_virial:
        c = fzeros((6,6))
        c0 = fzeros((6,6))
        pot.calc_elastic_constants(at, c=c, c0=c0)
        print 'CIJ_VIRIAL_C'
        print c.round(2)*EV_A3_IN_GPA
        print 'C_ij^0 (virial) / GPa ='
        print c0.round(2)*EV_A3_IN_GPA
        print


    if opt.cij_fit:
        c = elastic_constants(pot, at, opt.cij_symmetry, relax=True)
        c0 = elastic_constants(pot, at, opt.cij_symmetry, relax=False)
        print 'CIJ_FIT_C'
        print c.round(2)
        print 'CIJ_FIT_C0'
        print c0.round(2)
        print

    if opt.cij_energy:
        cs = {}

        strain_patterns = {
           'C11':        [1,0,0,0,0,0],
           'C11+C12':    [1,1,0,0,0,0],
           'C11+C13':    [1,0,1,0,0,0],
           'C33':        [0,0,1,0,0,0],
           'C44':        [0,0,0,1,0,0],
           'C66':        [0,0,0,0,0,1]
           }
        
        for label, strain_pattern in strain_patterns.iteritems():
            x = np.linspace(0.995, 1.005, 10)
            energies = []
            for xi in x:
                t = strain_matrix((xi-1.0)*np.array(strain_pattern))
                strained_at = transform(at, t)
                strained_at.set_cutoff(pot.cutoff()+1.0)
                strained_at.calc_connect()
                verbosity_push_decrement()
                pot.minim(strained_at, 'cg', 1e-3, 100,
                          do_pos=True, do_lat=False, use_n_minim=True)
                verbosity_pop()
                pot.calc(strained_at, energy=True)
                energies.append(strained_at.energy)

            a, b, c = np.polyfit(x, energies, 2)
            C = 2*a/v0*EV_A3_IN_GPA/sum(strain_pattern)
            cs[label] = C

        cs['C12'] = cs['C11+C12']-cs['C11']
        cs['C13'] = cs['C11+C13']-cs['C11']
        cs['B'] = 1./3*(cs['C11']+2*cs['C12']) # bulk modulus

        print 'CIJ_ENERGY'
        for label in ['B', 'C11', 'C33', 'C44', 'C66', 'C12', 'C13']:
            print '%s %.1f GPa' % (label, cs[label])

      


    if opt.surface_energy:
        if opt.surface:
            axes = crack_parse_name(opt.surface)
            m = rotation_matrix(at, axes[:,2], axes[:,3])
            bulk = orthorhombic_slab(at, rot=m, verbose=False)

            surface = bulk.copy()
            surface.lattice[2,2] = surface.lattice[2,2] + opt.surface_vacuum
            surface.set_lattice(surface.lattice, False)

        else:
            if not opt.at_surface:
               p.error('If --surface is not given, --at-surface must both be present.')
            bulk = at.copy()
            surface = Atoms(opt.at_surface)


        bulk.set_cutoff(pot.cutoff())
        bulk.calc_connect()
        surface.set_cutoff(pot.cutoff())
        surface.calc_connect()

        pot.calc(bulk, energy=True)
        print 'BULK'
        bulk.write('stdout')

        if opt.relax_surface:
            verbosity_push_decrement()
            pot.minim(surface, 'cg', 1e-6, 100, do_pos=True, do_lat=False)
            verbosity_pop()
            surface.write('relaxed-surface.xyz')

        pot.calc(surface, energy=True)
        print 'SURFACE'
        surface.write('stdout')

        if opt.surface_normal == 'x':
            surface_area = surface.lattice[2,2]*surface.lattice[3,3]
        elif opt.surface_normal == 'y':
            surface_area = surface.lattice[1,1]*surface.lattice[3,3]
        elif opt.surface_normal == 'z':
            surface_area = surface.lattice[1,1]*surface.lattice[2,2]
        else:
            raise RuntimeError('surface-normal should be one of x, y, z')

        gamma = (surface.energy - bulk.energy/bulk.n*surface.n)/(2.0*surface_area)*J_PER_M2
        print 'SURFACE_ENERGY'
        print 'Surface energy: gamma = ', gamma, ' J/m^2'

except RuntimeError, re:
    p.error(str(re))

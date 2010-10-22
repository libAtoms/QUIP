#!/usr/bin/env python

from quippy import *
import optparse

p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-s', '--struct', action='store', help='Structure name (mandatory)')
p.add_option('-o', '--outfile', action='store', help='Output file name (default stdout)', default='stdout')
p.add_option('-v', '--vol-per-atom', action='store', type='float', help='Volume per atom (default -1)', default=-1.0)
p.add_option('-V', '--vol-per-unit-cell', action='store', type='float', help='Volume per unit cell (default -1)', default=-1.0)
p.add_option('-r', '--repeat', action='store', help='Supercell repeat (default "1 1 1")', default='1 1 1')
p.add_option('-z', '--z-values', action='store', help='Atomic number values (1 per atom in primitive cell)')

opt, args = p.parse_args()

if len(args) != 0:
   p.error('No arguments are required')

if opt.repeat is not None:
    try:
        opt.repeat = [int(x) for x in opt.repeat.split()]
    except ValueError:
        p.error('Cannot parse repeat argument %s' % opt.repeat)

try:
    cell = structure_from_file(opt.struct, opt.vol_per_atom, opt.vol_per_unit_cell, opt.repeat, opt.z_values)
except RuntimeError, re:
    p.error(str(re))

cell.write(opt.outfile)

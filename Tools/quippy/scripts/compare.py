#!/usr/bin/env python

from quippy import *
import itertools
import sys
import optparse

TOL = 1e-5

p = optparse.OptionParser(usage='%prog [options] trajectory_A trajectory_B')
p.add_option('-O', '--offset', action='append', type='int',
             help='Frame offset, should be given once for each trajectory')
p.add_option('-P', '--params', action='append',
             help='List of params to compare. Can appear multiple times. Default is empty list.')
p.add_option('-p', '--properties', action='append',
             help='List of properties to compare. Can appear multiple times. Default is positions and forces.')
p.add_option('-m', '--mask', action='store',
             help='Name of property which should be non-zero in region of interest.')
opt, args = p.parse_args()

if len(args) != 2:
    p.error("Number of input trajectories must be 2")

if opt.offset is not None:
    if len(opt.offset) != len(args):
        p.error('Number of offset args should match number of input files')
else:
    opt.offset = [ 0 for i in args ]

if opt.params is None:
    opt.params = []

if opt.properties is None:
    opt.properties = ['pos', 'force']

a0 = Atoms(args[0], frame=opt.offset[0])
if opt.mask is not None:
    mask = getattr(a0, opt.mask) != 0
    atom_range = (mask.nonzero()[0].min(), mask.nonzero()[0].max())
    print 'Restricting atom range to %d..%d' % atom_range
else:
    atom_range = None

field_names = ['Frame']
fmts = ['%12d']
for prop in opt.params + opt.properties:
    value = farray(getattr(a0, prop))
    if value.shape == ():
        field_names.append('%s error' % prop)
        fmts.append('%18.6f')
    else:
        field_names.append('RMS %s error' % prop)
        fmts.append('%18.6f')
        field_names.append('Max %s error' % prop)
        fmts.append('%18.6f')

fmt = ''.join(fmts)
titles = ('%-18s'*len(field_names)) % tuple(field_names)


A = AtomsReader(args[0], range=atom_range, start=opt.offset[0])
B = AtomsReader(args[1], range=atom_range, start=opt.offset[1])

print titles
print '-'*len(titles)
frame = 0
for a, b in itertools.izip(A, B):
    data = []
    data.append(frame)
    frame += 1

    for param in opt.params:
        a_value = farray(getattr(a, param))
        b_value = farray(getattr(b, param))
        if a_value.shape == ():
            data.append(abs(float(a_value - b_value)))
        else:
            data.append(rms_diff(a_value, b_value))
            data.append(abs(a_value-b_value).max())

    if opt.mask is not None:
        mask = getattr(a, opt.mask) != 0
    else:
        mask = Ellipsis
    for prop in opt.properties:
        a_prop = getattr(a, prop)[mask]
        b_prop = getattr(b, prop)[mask]
        data.append(rms_diff(a_prop, b_prop))
        data.append(abs(a_prop - b_prop).max())

    print fmt % tuple(data)

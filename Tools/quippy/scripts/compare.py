from quippy import *
import itertools
import sys
import optparse

TOL = 1e-5

p = optparse.OptionParser(usage='%prog [options] trajectory_A trajectory_B')
p.add_option('-O', '--offset', action='append', type='int',
             help='Frame offset, should be given once for each trajectory')
p.add_option('-p', '--properties', action='append',
             help='List of properties to compare. Can appear multiple times. Default is positions and forces.')
p.add_option('-m', '--mask', action='store',
             help='Name of property which should be non-zero in first frame region of interest.', default='hybrid')
opt, args = p.parse_args()

if len(args) != 2:
    p.error("Number of input trajectories must be 2")

if opt.offset is not None:
    if len(opt.offset) != len(args):
        p.error('Number of offset args should match number of input files')
else:
    opt.offset = [ 0 for i in args ]

if opt.properties is None:
    opt.properties = ['pos', 'force']

field_names = ['Time']
fmts = ['%12.3f']
for prop in opt.properties:
    field_names.append('RMS %s error' % prop)
    fmts.append('%18.6f')
    field_names.append('Max %s error' % prop)
    fmts.append('%18.6f')

fmt = ''.join(fmts)
titles = ('%-18s'*len(field_names)) % tuple(field_names)

a0 = Atoms(args[0], frame=opt.offset[0])
mask = getattr(a0, opt.mask) != 0
atom_range = (mask.nonzero()[0].min(), mask.nonzero()[0].max())
print 'Restricting atom range to %d..%d' % atom_range

A = AtomsReader(args[0], range=atom_range, start=opt.offset[0])
B = AtomsReader(args[1], range=atom_range, start=opt.offset[1])

print titles
print '-'*len(titles)
for a, b in itertools.izip(A, B):
    if abs(a.time - b.time) > TOL:
        p.error("temporal misalignment! %f != %f\n" % (a.time, b.time))

    mask = getattr(a, opt.mask) != 0
    data = [a.time]
    for prop in opt.properties:
        a_prop = getattr(a, prop)[mask]
        b_prop = getattr(b, prop)[mask]
        data.append(rms_diff(a_prop, b_prop))
        data.append(abs(a_prop - b_prop).max())

    print fmt % tuple(data)

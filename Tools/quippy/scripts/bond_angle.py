#!/usr/bin/env python

import sys
import os
import optparse
import itertools
from quippy import AtomsReader

p = optparse.OptionParser(usage='%prog [options] <input file>...')
p.add_option('-b', '--bond', action='append', type='int', nargs=2,
             help='Pair of atoms. Can appear multiple times', default=[])
p.add_option('-a', '--angle', action='append', type='int', nargs=3,
             help='Triplet of atoms. Can appear multiple times', default=[])
p.add_option('-m', '--min-image', action='store_true',
             help='Apply minimum image convention')
p.add_option('-O', '--offset', action='append', type='int',
             help='Frame offset (one per input file)')
opt, args = p.parse_args()

if opt.offset is not None:
    if len(opt.offset) != len(args):
        p.error('Number of -O/--offset options must match number of input files')
else:
    opt.offset = [0]*len(args)

if len(opt.bond) == 0 and len(opt.angle) == 0:
    p.error('One or more -b/--bond pairs or -a/--angle triplets needed')

atoms = set([idx for bond in opt.bond for idx in bond] +
            [idx for angle in opt.angle for idx in angle])
atoms_start = min(atoms)
atoms_stop  = max(atoms)
sources = [AtomsReader(arg, start=offset, range=[atoms_start, atoms_stop]) for arg, offset in zip(args, opt.offset)]

bonds =  [ (i-atoms_start+1, j-atoms_start+1) for (i, j) in opt.bond ]
angles =  [ (i-atoms_start+1, j-atoms_start+1, k-atoms_start+1) for (i, j, k) in opt.angle ]

fmt_string = '%18.6f'*(1 + len(sources)*(len(bonds) + len(angles)))
for configs in itertools.izip(*sources):
    times = []
    bondlengths = []
    angles = []
    for at in configs:
        times.append(at.time)
        for i, j in bonds:
            if opt.min_image:
                bondlength = at.distamce_min_image(i, j)
            else:
                bondlength = (at.pos[:,j] - at.pos[:,i]).norm() 
            bondlengths.append(bondlength)

        for i, j, k in angles:
            if opt.min_image:
                u_ij = at.diff_min_image(i, j)
                u_ik = at.diff_min_image(i, k)
            else:
                u_ij = at.pos[:,j] - at.pos[:,i]
                u_ik = at.pos[:,k] - at.pos[:,i]

            u_ij = u_ij/u_ij.norm()
            u_ik = u_ik/u_ik.norm()
            angle = np.acos(np.dot(u_ij, u_ik))
            angles.append(angle)
                            
    if times.count(times[0]) != len(times):
        p.error('Misalignment: times=%=s' % times)

    print fmt_string % tuple(times[:1] + bondlengths)

#!/usr/bin/env python

import sys
import os
import glob
import optparse

p = optparse.OptionParser(usage='%prog [options] <image pattern>...')
p.add_option('-o', '--output', action='store', help='Output image pattern',
             default='out%05d.jpg')
p.add_option('-m', '--montage-opts', action='store', help='Extra options to pass to montage',
             default='-frame 5 -shadow -geometry +5+5')
p.add_option('-O', '--offset', action='append', type='int', help='Frame offset, should be given once for each pattern')
opt, args = p.parse_args()

globs = []
for pattern in args:
    globs.append(sorted(glob.glob(pattern)))

if opt.offset is not None:
    if len(opt.offset) != len(globs):
        p.error('Number of offset args should match number of input patterns')
    for i in range(len(globs)):
        globs[i] = globs[i][opt.offset[i]:]

#sorted_globs = [ glob for length, glob in sorted([ (len(glob), glob) for glob in globs ]) ]

for i, imgs in enumerate(zip(*globs)):
    infiles = ' '.join(imgs)
    outfile = opt.output % i
    command = 'montage %s -tile 1x%d %s %s' % (opt.montage_opts, len(imgs), infiles, outfile)
    print command
    os.system(command)


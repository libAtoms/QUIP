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

from quippy import *
import sys, optparse, itertools, os
from math import log10, floor, ceil

p = optparse.OptionParser(usage='%prog [options] <input file> [<output file>]')

p.add_option('-r', '--range', action='store', help="""Range of frames to include. Should be either a single frame
number or a slice [start]:[stop][:step]. If -r is omitted,
default is all frames. We use Fortran indexing, so slices include the
"stop" frame. Negative indices count backwards from the end of the file,
with -1 being the last frame.""")
p.add_option('-f', '--format', action='store', help="""Output format, specified by file extension. If format not given,
output format will be the same as input format.""")
p.add_option('-s', '--stem', action='store', help="""Output stem. If not specified, stem is basename of first input file.""")
p.add_option('-n', '--ndigit', action='store', help="""Number of digits in output filenames. If not given, try to count
number of frames in input file(s).""")
p.add_option('-t', '--test', action='store_true', help="""Test mode. Do not write any files, but print filenames.""")

opt, args = p.parse_args()

if len(args) == 0:
   p.error('One or more input files must be specified.')

if opt.range is not None:
   try:
      opt.range = parse_slice(opt.range)
   except:
      p.error('Cannot parse slice "%s" - should be in format [start]:[stop][:step]')
else:
   # Default is all frames
   opt.range = slice(1, None, None)

basename, ext = os.path.splitext(args[0])

if opt.stem is None: opt.stem = basename
if opt.format is None: opt.format = ext[1:]

sources = [AtomsReader(f) for f in args]

if isinstance(opt.range, slice):
    # Try to count number of frames
    nframes = sum([len(s) for s in sources])
    nframes = len(range(*opt.range.indices(nframes)))

    if opt.ndigit is None:
        if nframes is not None:
            opt.ndigit = int(ceil(log10(nframes+1)))
        else:
            opt.ndigit = 5

    out_fmt = '%s%%0%dd.%s' % (opt.stem, opt.ndigit, opt.format)

    # Build a chain of iterators over all input files, skipping frames as appropriate
    atomseq = itertools.islice(itertools.chain.from_iterable(sources),
                               opt.range.start-1, opt.range.stop, opt.range.step)

    for i, at in enumerate(atomseq):
       outfile = out_fmt % i
       if opt.test:
          print outfile
       else:
          at.write(outfile)
else:
    # Single frame
    at = list(itertools.chain.from_iterable(sources))[opt.range-1]
    outfile = '%s.%s' % (opt.stem, opt.format)
    if opt.test:
       print outfile
    else:
       at.write(outfile)

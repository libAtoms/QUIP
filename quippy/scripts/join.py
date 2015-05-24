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
from numpy import *
import sys, optparse, itertools, os
from math import log10, floor

p = optparse.OptionParser(usage='%prog [options] <input file> [<output file>]')

p.add_option('-f', '--format', action='store', help="""Output format, specified by file extension. If format not given,
output format is inferred from OUTPUT if present, otherwise from first input file.""")
p.add_option('-o', '--output', action='store', help="""Output file name, If not specified basename(INPUT)+FORMAT is used.""")

opt, args = p.parse_args()

if len(args) == 0:
   p.error('One or more input files must be specified.')

if opt.output is None: 
    basename, ext = os.path.splitext(args[0])
    if opt.format is None: opt.format = ext[1:]
    opt.output = basename + '.' + opt.format
else:
    basename, ext = os.path.splitext(opt.output)
    if opt.format is None: opt.format = ext[1:]

# Build a chain of iterators over all input files
atomseq = itertools.chain(*[AtomsReader(f) for f in args])

al = AtomsList(atomseq)
al.write(opt.output, format=opt.format)
        

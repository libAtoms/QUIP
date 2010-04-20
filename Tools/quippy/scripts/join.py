#!/usr/bin/env python
from quippy import *
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
atomseq = itertools.chain.from_iterable([AtomsList(f, store=False) for f in args])

al = AtomsList(atomseq)
al.write(opt.output, format=opt.format)
        

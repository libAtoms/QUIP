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
import sys

if len(sys.argv[1:]) == 0:
    print 'Usage: count.py <INPUT FILE>... - counts total number of frames in input files.'
    sys.exit(1)

args = sys.argv[1:]
sources = [AtomsReader(f) for f in args]

nframes = 0
for s in sources:
    # Try to get length from each source. If it doesn't support length, we have
    # to explicitly load all frames from the file.
    try:
        nframes += len(s)
    except ValueError:
        for at in s:
            nframes += 1

print nframes    


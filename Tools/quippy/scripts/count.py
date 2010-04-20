#!/usr/bin/env python

from quippy import *
import sys

if len(sys.argv[1:]) == 0:
    print 'Usage: count.py <INPUT FILE>... - counts total number of frames in input files.'
    sys.exit(1)

args = sys.argv[1:]
sources = [AtomsList(f, store=False) for f in args]

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


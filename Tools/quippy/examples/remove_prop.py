#!/usr/bin/env python

from quippy import *

if len(sys.argv[1:]) != 3:
   print 'Usage: <property name> <property> <infile> <outfile>'
   sys.exit(1)

prop = sys.argv[1]
a = Atoms(sys.argv[2])

if not a.has_property(prop):
   print 'Property "%s" not found' % prop
   sys.exit(1)

a.remove_property(prop)
a.write(sys.argv[3])

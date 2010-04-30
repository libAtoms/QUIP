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

if len(sys.argv[1:]) < 2:
   print 'Usage: param2xml gen.in Minimisation_progress [xml_file]'
   sys.exit(1)

xml_file = sys.stdout
if len(sys.argv[1:]) >  2:
   xml_file = open(sys.argv[3], 'w')


costs, params = sio2.read_minimisation_progress(sys.argv[2])

params = sio2.read_gen_and_param_files(sys.argv[1], params[-1], verbose=False)

xml_file.write(sio2.param_to_xml(params))



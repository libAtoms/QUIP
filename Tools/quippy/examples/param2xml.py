#!/usr/bin/env python

from quippy import *
import sys

if len(sys.argv[1:]) < 2:
   print 'Usage: param2xml gen.in Minimisation_progress [xml_file]'

xml_file = sys.stdout
if len(sys.argv[1:]) >  2:
   xml_file = open(sys.argv[3], 'w')


costs, params = sio2.read_minimisation_progress(sys.argv[2])

params = sio2.read_gen_and_param_files(sys.argv[1], params[-1], verbose=False)

xml_file.write(sio2.param_to_xml(params))



#!/usr/bin/env python

from pyatoms import *
import sys, os

# Determine output format from command name
two_pos = sys.argv[0].find('2')
outformat = sys.argv[0][two_pos+1:two_pos+4]

if not outformat in ['xyz', 'cfg']:
   print "Don't know how to convert to files of type %s" % outformat
   sys.exit(1)


if sys.argv[1] == '-p':
   properties = ['species','pos'] + sys.argv[2].split(':')
   files = sys.argv[3:]
else:
   properties = None
   files = sys.argv[1:]

for filename in files:
   root, ext = os.path.splitext(filename)
   fr = frame_reader(filename)

   if outformat == 'xyz':
      outf = open(root+'.xyz','w')

   for i, a in enumerate(fr):

      sys.stdout.write('Frame %d (%d atoms)               \r' % (i, a.n))
      sys.stdout.flush()

      force_props = filter(lambda name: name.endswith('force'), a.properties.keys())
      for force_prop in force_props:
         a.add_property('norm_'+force_prop, eval('norm(a.'+force_prop+')'))
         if properties is not None and not 'norm_'+force_prop in properties:
            properties.append('norm_'+force_prop)

      if outformat == 'xyz':
         a.write_xyz(outf, properties=properties)
      elif outformat == 'cfg':
         a.write_cfg('%s%05d.cfg' % (root, i), properties=properties)

   print
   
   if outformat == 'xyz':
      outf.close()

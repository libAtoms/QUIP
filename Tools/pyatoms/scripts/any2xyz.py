#!/usr/bin/env python
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HP X
# HP X   pyatoms: atomistic simulations tools
# HP X
# HP X   Copyright James Kermode 2010
# HP X
# HP X   These portions of the source code are released under the GNU General
# HP X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HP X
# HP X   If you would like to license the source code under different terms,
# HP X   please contact James Kermode, james.kermode@gmail.com
# HP X
# HP X   When using this software, please cite the following reference:
# HP X
# HP X   http://www.jrkermode.co.uk/PyAtoms
# HP X
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

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

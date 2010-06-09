#!/usr/bin/env python

from quippy import *
import optparse, sys

p = optparse.OptionParser(usage='%prog [option] [file1] [file2] ...')
p.add_option('-b', '--bond_by_cutoff', action='store_true', default=False, help="""add bond_by_cutoff statement""")
p.add_option('-t', '--type', action='store', help="""property to use for atom type""")
p.add_option('-v', '--value', action='append', help="""property to use for values (may be repeated)""")
p.add_option('-p', '--post_config_command', action='append', help="""commands to add after each config""")
p.add_option('-P', '--post_file_command', action='append', help="""commands to add after each input file""")
p.add_option('-e', '--end_command', action='append', help="""commands to add at the very end""")

opt, args = p.parse_args()

if (len(args) == 0):
  args = ["stdin"]

for file in args:
  if (file == "-"):
    file = "stdin"

  for at in AtomsReader(file):
    print "new_configuration"
    if (hasattr(at,'lattice')):
      print "pbc_a 1 %f %f %f" % ( at.lattice[1,1], at.lattice[2,1], at.lattice[3,1])
      print "pbc_a 2 %f %f %f" % ( at.lattice[1,2], at.lattice[2,2], at.lattice[3,2])
      print "pbc_a 3 %f %f %f" % ( at.lattice[1,3], at.lattice[2,3], at.lattice[3,3])
    if (hasattr(at,'n')):
      for i_at in frange(at.n):
	if (opt.type is not None):
	  type = opt.type
	else:
	  if (hasattr(at,'z')):
	    type = 'z'
	  else:
	    if (hasattr(at,'species')):
	      type = 'species'
	    else:
	      if (hasattr(at,'type')):
		type = 'type'
	      else:
		print "Can't find z, species, or type for atom type"
		os.exit(1)
	print "atom %f %f %f %s" % (at.pos[1,i_at], at.pos[2,i_at], at.pos[3,i_at], getattr(at,type)[i_at]),
	if (opt.value is not None):
	  for iv in range(len(opt.value)):
	    print " value %d %s" % (iv+1, getattr(at,opt.value[iv])[i_at]),
	print ""
    if opt.bond_by_cutoff:
      print "bond_by_cutoff"
    if (opt.post_config_command is not None):
      for cmd in opt.post_config_command:
	print cmd

  if (opt.post_file_command is not None):
    for cmd in opt.post_file_command:
      print cmd

if (opt.end_command is not None):
  for cmd in opt.end_command:
    print cmd

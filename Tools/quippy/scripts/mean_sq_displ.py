#!/usr/bin/env python

from quippy import *
import sys

if (len(sys.argv) != 2):
   sys.stderr.write("Usage: %s infile\n" % sys.argv[0])
   sys.exit(-1)

if (sys.argv[1] == "-h" or sys.argv[1] == "--help"):
   sys.stderr.write("Usage: %s infile\n" % sys.argv[0])
   sys.exit(-1)

ar = AtomsReader(sys.argv[1])
a0 = None
for a_in in ar:
   if (a0 is None):
      # first iteration
      a0 = a_in.copy()
      a = a_in.copy()
      if (not hasattr(a, 'mass')):
	 a.add_property('mass', 0.0, n_cols=1)
	 a.set_atoms(a.Z[:])
      if (not hasattr(a, 'orig_pos')):
	 a.add_property('orig_pos', a0.pos)
	 a.params['orig_CoM'] = a.centre_of_mass(origin=0)
	 orig_CoM = a.params['orig_CoM']
      if (not hasattr(a, 'msd_displ')):
	 a.add_property('msd_displ', 0.0, n_cols=3)
      delta_CoM = fzeros( (3,a.n) )

   # copy from read in struct to persistent struct
   a.pos[:,:] = a_in.pos[:,:]
   a.set_lattice (a_in.lattice, False)

   # undo pbc jumps
   a.undo_pbc_jumps()

   # fix CoM
   CoM = a.centre_of_mass(origin=0)
   delta_CoM[1,:] = CoM[1] - orig_CoM[1]
   delta_CoM[2,:] = CoM[2] - orig_CoM[2]
   delta_CoM[3,:] = CoM[3] - orig_CoM[3]
   a.pos[:,:] = a.pos[:,:] - delta_CoM[:,:]

   # calc msd
   a.msd_displ = a.orig_pos - a.pos
   msd = sum(sum((a.pos-a.orig_pos)**2))/float(a.n)
   a.params['msd'] = msd

   a.write("stdout")

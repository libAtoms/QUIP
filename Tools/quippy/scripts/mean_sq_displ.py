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
for a in ar:
   if (a0 is None):
      # first iteration
      a0 = a.copy()
      if (not hasattr(a, 'mass')):
	 a0.add_property('mass', 0.0, n_cols=1)
	 a0.set_atoms(a0.Z[:])
      if (not hasattr(a0, 'orig_pos')):
	 a0.add_property('orig_pos', a0.pos)
      a0.params['orig_CoM'] = a0.centre_of_mass(origin=0)
      orig_CoM = a0.params['orig_CoM']
      delta_CoM = fzeros( (3,a0.n) )
   else:
      if (hasattr(a, 'prev_pos')):
	 a.prev_pos[:,:] = prev_pos[:,:]
      else:
	 a.add_property('prev_pos', prev_pos)

   if (hasattr(a, 'orig_pos')):
      a.orig_pos[:,:] = a0.orig_pos[:,:]
   else:
      a.add_property('orig_pos', a0.orig_pos)
   a.params['orig_CoM'] = orig_CoM

   if (not hasattr(a, 'msd_displ')):
      a0.add_property('msd_displ', 0.0, n_cols=3)

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

   # save prev_pos
   prev_pos = a.prev_pos.copy()

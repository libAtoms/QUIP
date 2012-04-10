#!/usr/bin/env python

from quippy import *
import sys;

################################################################################

if (len(sys.argv) != 4):
   sys.stderr.write("Usage: %s infile outfile min_step\n" % sys.argv[0])
   sys.exit(1)

ar=AtomsReader(sys.argv[1])
aw=AtomsWriter(sys.argv[2])
msd_min_step=int(sys.argv[3])

iter=0
for at_cur in ar:
   if (iter == 0): # create new atoms structure
      at = at_cur.copy()
      at.add_property("prev_pos", at.pos)
      at.add_property("orig_pos", at.pos)
      at.add_property("msd_displ", 0.0, n_cols=3)
      at.add_property('mass', 0.0, n_cols=1)
      at.set_atoms(at.Z[:])
      at.params['orig_CoM'] = at.centre_of_mass(origin=0)
      at.params['step'] = 0
      delta_CoM = fzeros( (3, at.n) )
   else: # copy current pos into atoms structure
      at.pos[:,:] = at_cur.pos

   at.params['step'] += 1

   # print "at.pos[:,1] pre undo ", at.pos[:,1]
   at.undo_pbc_jumps()
   # print "at.pos[:,1] post undo ", at.pos[:,1]

   at.prev_pos[:,:] = at.pos

   # shift CoM back to where it was
   CoM = at.centre_of_mass(origin=0)
   # set delta_CoM from current CoM, orig_CoM
   delta_CoM[1,:] = CoM[1] - at.params['orig_CoM'][1]
   delta_CoM[2,:] = CoM[2] - at.params['orig_CoM'][2]
   delta_CoM[3,:] = CoM[3] - at.params['orig_CoM'][3]
   # finally move CoM back to where it belongs
   at.pos -= delta_CoM

   # print "orig CoM ", at.params['orig_CoM'], " cur CoM ", CoM
   # print "at.pos[:,1] post CoM shift ", at.pos[:,1]
   # print "final CoM ", at.centre_of_mass(origin=0)

   if (at.params['step'] <= msd_min_step):
      at.orig_pos[:,:] = at.pos[:,:]

   at.msd_displ[:,:] = at.orig_pos-at.pos
   msd = sum(sum((at.pos-at.orig_pos)*(at.pos-at.orig_pos)))/float(at.n)

   # save msd
   at.params['msd'] = msd

   # map_into_cell(at.pos, at.lattice, at.g)

   aw.write(at)

   if ((iter % 100) == 0):
      sys.stderr.write("%1d" % (iter/100))
   elif ((iter % 10) == 0):
      sys.stderr.write(".")

   iter += 1

sys.stderr.write("\n")

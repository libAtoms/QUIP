#!/usr/bin/env python

from quippy import *
import sys;

################################################################################

if (len(sys.argv) != 4 and len(sys.argv) != 5):
   sys.stderr.write("Usage: %s infile outfile min_step [ mask ]\n" % sys.argv[0])
   sys.exit(1)

ar=AtomsReader(sys.argv[1])
aw=AtomsWriter(sys.argv[2])
msd_min_step=int(sys.argv[3])
if (len(sys.argv) == 5):
   mask_str = sys.argv[4]
else:
   mask_str = ""

step=0
for at in ar:
   step += 1

   at.undo_pbc_jumps(persistent=False)
   at.undo_com_motion(persistent=False)
   if (mask_str != ""):
      at.calc_msd(mask=eval(mask_str), reset_msd = (step <= msd_min_step), persistent=False)
   else:
      at.calc_msd(reset_msd = (step <= msd_min_step), persistent=False)

   aw.write(at)

   if ((step % 100) == 0):
      sys.stderr.write("%1d" % (step/100))
   elif ((step % 10) == 0):
      sys.stderr.write(".")

sys.stderr.write("\n")

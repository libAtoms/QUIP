#!/usr/bin/env python

from quippy import *
import sys, optparse;

################################################################################

p = optparse.OptionParser(usage='%prog [ options ]')
p.add_option('-i','--infile', action='store', help="""input file""", default="stdin")
p.add_option('-o','--outfile', action='store', help="""output file""", default="stdout")
p.add_option('-s','--min_step', action='store', help="""frame to start from""", type='int', default=0)
p.add_option('-I','--interval', action='store', help="""interval between calculated configs""", type='int', default=1)
p.add_option('-m','--mask', action='store', help="""mask of atoms to include (must resolve to python logical array)""", default="")
p.add_option('-f','--fake_smooth_pos_mixing', action='store', help="""mixing coefficient for fake smooth pos""", type='float', default=-1.0)

opt, args = p.parse_args()

ar = AtomsReader(opt.infile, step=opt.interval)
aw = AtomsWriter(opt.outfile)

step=0
eff_step=0
for at in ar:
   at.undo_pbc_jumps(persistent=False)
   at.undo_com_motion(persistent=False)
   if (opt.mask != ""):
      at.calc_msd(mask=eval(opt.mask), reset_msd = (step <= opt.min_step), persistent=False)
   else:
      at.calc_msd(reset_msd = (step <= opt.min_step), persistent=False)

   if (opt.fake_smooth_pos_mixing > 0.0):
      at.fake_smooth_pos(opt.fake_smooth_pos_mixing, persistent=False)

   aw.write(at)

   step += opt.interval
   eff_step += 1

   if ((eff_step % 100) == 0):
      sys.stderr.write("%1d" % (eff_step/100))
   elif ((eff_step % 10) == 0):
      sys.stderr.write(".")


if (eff_step >= 10):
   sys.stderr.write("\n")

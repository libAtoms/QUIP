#!/usr/bin/env python

from quippy import *
import sys, optparse;
import numpy as np

p = optparse.OptionParser(usage='%prog [ options ]')
p.add_option('-i','--infile', action='store', help="""input file""", default="stdin")
p.add_option('-o','--outfile', action='store', help="""output file""", default="stdout")
p.add_option('-s','--min_step', action='store', help="""frame to start from""", type='int', default=0)
p.add_option('-I','--interval', action='store', help="""interval between calculated configs""", type='int', default=1)

opt, args = p.parse_args()

if opt.infile == 'stdin' or opt.infile == '-':
   sys.stderr.write("infile can't be stdin or -\n")
   sys.exit(1)

ar=AtomsReader(opt.infile, step=opt.interval)
mean_pos = None
n=0
step=0
step_eff=0
for at in ar:
   if (step >= opt.min_step):
      if mean_pos is None:
	 mean_pos = at.pos.copy()
	 n=1
      else:
	 mean_pos[:,:] += at.pos[:,:]
	 n += 1
   step += opt.interval
   step_eff += 1

   if (step_eff % 100) == 0:
      sys.stderr.write("%d" % (step_eff/100))
   else:
      if (step_eff % 10) == 0:
	 sys.stderr.write(".")

sys.stderr.write("\n")

mean_pos /= float(n)

ar=AtomsReader(opt.infile, step=opt.interval)
aw=AtomsWriter(opt.outfile)
displ = None
n=0
step=0
step_eff=0
for at in ar:
   if (step >= opt.min_step):
      displ = at.pos - mean_pos
      if displ is None:
	 if not hasattr(at, "msd_displ"):
	    at.add_property("msd_displ", displ)
	 if not hasattr(at, "msd_displ_mag"):
	    at.add_property("msd_displ_mag", np.sum(displ*displ, axis=0))
      at.msd_displ[:,:] = displ
      at.msd_displ_mag[:] = np.sum(displ*displ, axis=0)
      aw.write(at)
   step += opt.interval
   step_eff += 1

   if (step_eff % 100) == 0:
      sys.stderr.write("%d" % (step_eff/100))
   elif (step_eff % 10) == 0:
	 sys.stderr.write(".")

sys.stderr.write("\n")

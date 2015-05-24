#!/usr/bin/env python

from quippy import *
from numpy import *
import optparse, sys

p = optparse.OptionParser(usage='%prog file; # splits file by frame, with all output files lattice rescaled isotropically to match the mean volume of the inputs')

opt, args = p.parse_args()

for file in args:
   print file
   v_mean = 0.0
   al = AtomsReader(file)
   n=0
   for at in al:
      v_mean = v_mean + at.cell_volume()
      n += 1
   v_mean = v_mean / real(n)
   nd = int(log(n)/log(10))+1
   fmt = ".%%0%dd" % nd
   al = AtomsReader(file)
   i=0
   for at in al:
      if (hasattr(at,"time")):
	 print "time=%f " % at.time,
      v = at.cell_volume()
      lat = at.lattice
      new_lat = lat*pow(v_mean/v,1.0/3.0)
      at.set_lattice(new_lat, scale_positions=True)
      dot_pos = file.rindex('.')
      new_file=file[:dot_pos]+(fmt % i)+".mean_vol."+file[dot_pos+1:]
      print new_file
      at.write(new_file)
      i += 1

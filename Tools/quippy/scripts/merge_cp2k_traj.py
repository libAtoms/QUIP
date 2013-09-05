#!/usr/bin/env python

import optparse, sys, re

p = optparse.OptionParser(usage='%prog file1 [file2 ... ]')
p.add_option('-c', '--check_time_steps', action='store_true', default=False, help="""check that all configs agree on i and time""")
p.add_option('-f', '--file', action='store', default="stdout", help="""output file""")
opt, args = p.parse_args()

from quippy import *
from itertools import izip

if (len(args) == 0):
  print "Need at least one input file"
  sys.exit(1)

ar = []
for i in range(len(args)):
  ar.append(AtomsReader(args[i]))

at=[]
for i in range(len(args)):
  at.append(None)

co = CInOutput(opt.file,OUTPUT)

config_i = 0
for configs_to_merge in izip(*ar):
  config_i += 1
  sys.stderr.write("%d\n" % config_i)
  for i in range(len(args)):
    if (i == 0):
      at = configs_to_merge[0].copy()

    # figure out property name - should be more clever
    if re.search("-pos[-.]",args[i]) is not None:
      in_name="pos"
      out_name="pos"
      unit_conv=1.0
    elif re.search("-vel[-.]",args[i]) is not None:
      in_name="velo"
      out_name="velo"
      unit_conv=HARTREE*BOHR/HBAR
    elif re.search("-frc[-.]",args[i]) is not None:
      in_name="frc"
      out_name="force"
      unit_conv=HARTREE/BOHR
    elif re.search("-fmlabels[-.]",args[i]) is not None:
      in_name=""
      out_name="label"
      unit_conv=1.0
    else:
      print "Don't know what property comes from file names '%s', aborting" % args[i]
      sys.exit(2)

    # check i and time, if needed
    if (opt.check_time_steps and i > 0):
      if (hasattr(configs_to_merge[i],"i") and hasattr(configs_to_merge[0],"i") and
	  configs_to_merge[i].i != configs_to_merge[0].i):
	print "config %d has mismatching i params %d != %d from files %s %s" % (config_i, configs_to_merge[0].i, configs_to_merge[i].i, 
	  args[0], args[i])
	sys.exit(3)
      if (hasattr(configs_to_merge[i],"time") and hasattr(configs_to_merge[0],"time") and
	  configs_to_merge[i].time != configs_to_merge[0].time):
	print "config %d has mismatching time params %f != %f from files %s %s" % (config_i, configs_to_merge[0].time, configs_to_merge[i].time, 
	  args[0], args[i])
	sys.exit(3)

    # add property
    at.add_property(out_name, unit_conv*getattr(configs_to_merge[i], in_name))
    # add params
    for k in configs_to_merge[i].params:
      at.params[k] = configs_to_merge[i].params[k]

  at.write(co)

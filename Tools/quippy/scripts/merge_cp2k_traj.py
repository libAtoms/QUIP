#!/usr/bin/env python

from quippy import *
import optparse, sys

p = optparse.OptionParser(usage='%prog file1 [file2 ... ]')
p.add_option('-c', '--check_time_steps', action='store_true', default=False, help="""check that all configs agree on i and time""")
p.add_option('-f', '--file', action='store', default="stdout", help="""output file""")
opt, args = p.parse_args()

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

##
configs_to_merge=[]
for i in range(len(ar)):
  configs_to_merge.append(None)
##

config_i = 0
# for configs_to_merge in zip(*ar):
for at0 in ar[0]:
  ##
  configs_to_merge[0] = at0
  for i in range(1,len(ar)):
    configs_to_merge[i] = ar[i].next()
  ##
  config_i += 1
  print config_i
  for i in range(len(args)):
    if (i == 0):
      at = configs_to_merge[0].copy()

    # figure out property name - should be more clever
    if (args[i].count("-pos-") == 1):
      in_name="pos"
      out_name="pos"
      unit_conv=1.0
    elif (args[i].count("-vel-") == 1):
      in_name="velo"
      out_name="velo"
      unit_conv=HARTREE*BOHR/HBAR
    elif (args[i].count("-frc-") == 1):
      in_name="frc"
      out_name="force"
      unit_conv=HARTREE*BOHR
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

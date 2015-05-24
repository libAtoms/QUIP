#!/usr/bin/env python

from quippy import *
import optparse

p=optparse.OptionParser(usage="%prog [option]")
p.add_option('-l','--library', action='store', help="""motif library file""")
p.add_option('-p','--pos', action='store', help="""positions property name for connectivity""")
p.add_option('-f','--file', action='store', help="""file to check, default stdin""")

opt, args = p.parse_args()

if (opt.file is not None):
  file=opt.file
else:
  file = "stdin"

if (file == "-"):
  file = "stdin"

a=Atoms(file)

if (opt.library is not None):
  a.params["Library"] = opt.library

a.calc_connect()

if (opt.pos is not None):
    create_residue_labels_arb_pos(a, do_charmm=True, pos_field_for_connectivity=opt.pos)
else:
  if (hasattr(a,'avgpos')):
    create_residue_labels_arb_pos(a, do_charmm=True)
  else:
    create_residue_labels_arb_pos(a, do_charmm=True, pos_field_for_connectivity='pos')

a.write("stdout")

#!/usr/bin/env python

from quippy import *

a=Atoms("stdin")
a.calc_connect()
if (hasattr(a,'avgpos')):
  create_residue_labels_arb_pos(a, do_CHARMM=True)
else:
  create_residue_labels_arb_pos(a, do_CHARMM=True, pos_field_for_connectivity='pos')

a.write("stdout")

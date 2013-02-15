#!/usr/bin/env python

from quippy import *
import optparse, sys

p = optparse.OptionParser(usage='%prog [option] [file1] [file2] ...')
p.add_option('-n', '--n_images', action='store', type='int', default=1, help="""number of images to add between each pair""")
opt, args = p.parse_args()

ar_i = Atoms(args[0])
ar_ip1 = ar_i.copy()
ar_ip1.undo_pbc_jumps()

for i in frange(len(args)-1):
   ar_in = Atoms(args[i])
   ar_ip1.set_lattice(ar_in.lattice)
   ar_ip1.pos[:,:] = ar_in.pos
   ar_ip1.undo_pbc_jumps()
   dlattice = ar_ip1.lattice - ar_i.lattice
   dpos = ar_ip1.pos - ar_i.pos

   if (i == 1):
      at_out = ar_i.copy()

   for ii in range(opt.n_images+1):
      dii = float(ii)/float(opt.n_images+1)
      at_out.set_lattice(ar_i.lattice+dii*dlattice, False)
      at_out.pos[:,:] = ar_i.pos + dii*dpos
      at_out.write("stdout")

   ar_i = ar_ip1.copy()

ar_i.write("stdout")

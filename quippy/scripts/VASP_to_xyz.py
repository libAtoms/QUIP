#!/usr/bin/python

from quippy import *
from numpy import *
import optparse, sys

p = optparse.OptionParser(usage='%prog [--species S1 S2] POSCAR [OUTCAR]')
p.add_option('-s','--species', action='store', help="""list of species to map to atom types, if needed""")

opt, args = p.parse_args()

if (len(args) == 0):
   print "POSCAR missing"
   sys.exit(1)

if (len(args) > 2):
   print "too many arguments"
   sys.exit(1)

poscar_file = args[0]
outcar = (len(args) == 2)
if (outcar):
   outcar_file = args[1]

p = open(poscar_file, 'r')
comment = p.readline().rstrip()
l = p.readline().strip(); lc_factor=real(l)
l = p.readline().strip(); a1 = real(l.split())
l = p.readline().strip(); a2 = real(l.split())
l = p.readline().strip(); a3 = real(l.split())
l = p.readline().strip(); species = l.split()
try:
   ns = [ int(n) for n in species ]
   no_species = True
except:
   no_species = False

have_species = True
if (no_species):
   for i in range(len(ns)):
      if (opt.species is not None):
	 species_cli = opt.species.split()
	 species[i] = species_cli[i]
      else:
	 have_species = False
	 species[i] = ("%d" % (i+1))
else:
   l = p.readline().strip();
   ns = [ int(n) for n in l.split() ]

dyn_type = p.readline().strip();
coord_type = p.readline().strip();

n=0
for i in range(len(ns)):
   n += ns[i]

lat = fzeros( (3,3) )
lat[:,1] = a1[0:3]
lat[:,2] = a2[0:3]
lat[:,3] = a3[0:3]

at = Atoms(n=n, lattice=lat)

ii = 1
for ti in range(len(ns)):
   for i in range(ns[ti]):
      l = p.readline().strip(); pos = l.split()[0:3]
      at.pos[:,ii] = pos[:]
      at.species[:,ii] = species[ti]
      ii += 1

if (have_species):
   at.set_zs()
else:
   at.Z[:] = [ int("".join(n)) for n in at.species[:] ]

at.write("stdout")

at_cur = at.copy()
lat_cur = at_cur.lattice.copy()
at_i = -1
lat_i = -1
if outcar:
   p = open(outcar_file, 'r')
   for lr in p:
      l=lr.rstrip()
      if (lat_i >= 1 and lat_i <= 3):
	 lat_cur[:,lat_i] = [ float(r) for r in l.replace("-"," -").split()[0:3] ]
	 lat_i += 1
      if (at_i >= 1 and at_i <= at_cur.n):
	 at_cur.pos[:,at_i] = [ float(r) for r in l.replace("-"," -").split()[0:3] ]
	 at_i += 1
      if (l.find("TOTAL-FORCE (eV/Angst)") >= 0):
	 at_i=1
	 p.next()
      if (l.find("direct lattice vectors") >= 0):
	 lat_i=1
      if (at_i == at_cur.n):
	 at_cur.set_lattice(lat_cur, False)
	 for i in frange(at_cur.n):
	    dr = at_cur.diff_min_image(i, at.pos[:,i])
	    at_cur.pos[:,i] = at.pos[:,i] - dr
	 at_cur.write("stdout")

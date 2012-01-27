# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import quippy
from quippy.atoms import Atoms, atoms_reader, AtomsReaders
from quippy.farray import fzeros,frange
import numpy as np
import re, sys

@atoms_reader('OUTCAR')
def VASP_POSCAR_Reader(outcar, species=None):
   """Read a configuration from a VASP OUTCAR file."""
   
   p = open(outcar, 'r')

   re_comment = re.compile("\s*POSCAR:\s*(.+)")
   re_potcar = re.compile("\s*POTCAR:\s*\S+\s+(\S+)")
   re_n_atoms = re.compile("\s*ions per type =\s*((?:\d+\s*)*)")

   energy_i = -1
   at_i = -1
   lat_i = -1
   elements=[]
   n_at = -1
   at_cur = None
   for lr in p:
      l=lr.rstrip()
      if (n_at <= 0):
	 # parse header type things
	 m = re_comment.match(l)
	 if (m is not None):
	    VASP_Comment=m.group(1)
	    # print "got VASP_Comment '%s'" % VASP_Comment
	 m = re_potcar.match(l)
	 if (m is not None):
	    elements.append(m.group(1))
	 m = re_n_atoms.match(l)
	 if (m is not None):
	    # print "got ions per type, groups are:"
	    # print m.groups()
	    lat = fzeros( (3,3) )
	    n_types = [ int(f) for f in m.group(1).split() ]
	    n_at=sum(n_types)
	    at = Atoms(n=n_at, latttice=lat)
	    i_at = 0
	    for type_i in range(len(n_types)):
	       for j in range(n_types[type_i]):
		  i_at += 1
		  # print "set species of atom %d to '%s'" % (i_at, elements[type_i])
		  at.species[i_at] = elements[type_i]
	    at.set_zs()
      else:
	 # parse per-config lattice/pos/force
         if (l.find("direct lattice vectors") >= 0):
	    at_cur = at.copy()
	    lat_cur=fzeros( (3,3) )
            lat_i=1
         elif (lat_i >= 1 and lat_i <= 3):
            lat_cur[:,lat_i] = [ float(r) for r in l.replace("-"," -").split()[0:3] ]
            lat_i += 1
         elif (l.find("TOTAL-FORCE (eV/Angst)") >= 0):
	    if (not hasattr(at_cur, "force")):
	       at_cur.add_property("force", 0.0, n_cols=3)
            at_i=1
            p.next()
         elif (at_i >= 1 and at_i <= at_cur.n):
	    pos_force = [ float(r) for r in l.replace("-"," -").split()[0:6] ]
            at_cur.pos[:,at_i] = pos_force[0:3]
            at_cur.force[:,at_i] = pos_force[3:6]
            at_i += 1
	 elif (l.find("free  energy") >= 0):
	    at_cur.params['Energy'] = float(l.split()[4])
	    energy_i = 1
	    p.next()
	 elif (energy_i == 1):
	    print "energy(sigma->0) line"
	    print l.split()
	    at_cur.params['Energy_sigma_to_zero'] = float(l.split()[6])
	    energy_i += 1
            yield at_cur
         if (at_cur is not None and at_i == at_cur.n):
            at_cur.set_lattice(lat_cur, False)
            for i in frange(at_cur.n):
               dr = at_cur.diff_min_image(i, at.pos[:,i])
               at_cur.pos[:,i] = at.pos[:,i] - dr


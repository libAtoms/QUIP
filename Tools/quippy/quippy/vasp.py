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
from quippy.atoms import Atoms, atoms_reader, AtomsReaders, AtomsWriters
from quippy.farray import fzeros,frange
import numpy as np
import re, sys

@atoms_reader('vasp')
@atoms_reader('POSCAR')
def VASPReader(poscar, outcar=None, species=None):
   """Read a configuration from a VASP POSCAR file.

   Following POSCAR, optionally also read a trajectory from an OUTCAR file."""
   
   p = open(poscar, 'r')
   comment = p.readline().rstrip()
   l = p.readline().strip(); lc_factor=np.real(l)
   l = p.readline().strip(); a1 = np.real(l.split())
   l = p.readline().strip(); a2 = np.real(l.split())
   l = p.readline().strip(); a3 = np.real(l.split())
   l = p.readline().strip(); at_species = l.split()
   try:
      ns = [ int(n) for n in at_species ]
      no_species = True
   except:
      no_species = False

   have_species = True
   if (no_species):
      for i in range(len(ns)):
         if (species is not None):
            species_cli = species.split()
            at_species[i] = species_cli[i-1]
         else:
            have_species = False
            at_species[i] = ("%d" % (i+1))
   else:
      l = p.readline().strip();
      ns = [ int(n) for n in l.split() ]

   l=p.readline().strip()
   if (re.compile("^\s*s", re.IGNORECASE).match(l)):
      dyn_type = l
      coord_type = p.readline().strip();
   else:
      coord_type = l

   n=0
   for i in range(len(ns)):
      n += ns[i]

   lat = fzeros( (3,3) )
   lat[:,1] = a1[0:3]
   lat[:,2] = a2[0:3]
   lat[:,3] = a3[0:3]

   at = Atoms(n=n, lattice=lat)
   at.params['VASP_Comment'] = comment

   coord_direct=re.compile("^\s*d", re.IGNORECASE).match(coord_type);

   ii = 1
   for ti in range(len(ns)):
      for i in range(ns[ti]):
         l = p.readline().strip(); pos = np.array(l.split()[0:3], float); 
	 if (coord_direct):
	    at.pos[:,ii] = np.dot(at.lattice[:,:],pos[:])
	 else:
	    at.pos[:,ii] = pos[:]
         at.species[:,ii] = at_species[ti]
         ii += 1

   if (have_species):
      at.set_zs()
   else:
      at.Z[:] = [ int("".join(n)) for n in at.species[:] ]

   yield at

   at_cur = at.copy()
   lat_cur = at_cur.lattice.copy()
   at_i = -1
   lat_i = -1

   if outcar is not None:
      p = open(outcar, 'r')
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
            yield at_cur

class VaspWriter(object):

   def __init__(self, out, species_list=None):
      self.out=out
      self.opened = False
      if type(self.out) == type(''):
	 if self.out == 'stdout':
	    self.out = sys.stdout
	 else:
	    self.out = open(self.out, 'w')
	    self.opened = True
      if species_list is not None:
	 self.species_list = species_list

   def close(self):
      if self.opened:
	 self.out.close()

   def write(self, at):

      # find numbers of atoms of each type, property name and value to match, and property labels to print above numbers
      atnums = []
      prop_vals = []
      prop_vals_map = {}
      labels = []
      if (hasattr(self, 'species_list')):
	 if (not hasattr(at, 'species')):
	    sys.stderr.write("VASP writer needs species property when species_list is specified")
	    sys.exit(1)
	 property='species'
	 for s in self.species_list:
	    labels.append(s)
	    prop_vals.append(s)
	    prop_vals_map[s] = 1
	    atnums.append((at.species[:].stripstrings() == s).count())
      else:
	 property='Z'
	 atnums_map = {}
	 for i_at in frange(at.n):
	    try:
	       atnums_map["%d" % at.Z[i_at]] += 1
	    except:
	       atnums_map["%d" % at.Z[i_at]] = 1

	 for Z_s in sorted(atnums_map.keys(), key = lambda entry: int(entry)):
	    prop_vals.append(int(Z_s))
	    prop_vals_map[int(Z_s)] = 1
	    labels.append(quippy.ElementName[int(Z_s)])
	    atnums.append(int(atnums_map[Z_s]))

      # check that every atom has a valid type
      for i_at in frange(at.n):
	 try:
	    prop = getattr(at,property)[i_at].stripstrings()
	 except:
	    prop = getattr(at,property)[i_at]

	 if not prop in prop_vals_map:
	    # should probably handle situation when prop isn't a string, but that should never happen
	    sys.stderr.write("Failed to find property %s in prop_vals_map dictionary" % prop)
	    sys.exit(1)

      # Comment
      try:
	 self.out.write(at.params['VASP_Comment']+"\n")
      except:
	 try:
	    self.out.write(at.params['comment']+"\n")
	 except:
	    self.out.write("\n")

      # Lattice
      self.out.write("1.0\n")
      self.out.write("%f %f %f\n" % (at.lattice[1,1], at.lattice[2,1], at.lattice[3,1]))
      self.out.write("%f %f %f\n" % (at.lattice[1,2], at.lattice[2,2], at.lattice[3,2]))
      self.out.write("%f %f %f\n" % (at.lattice[1,3], at.lattice[2,3], at.lattice[3,3]))

      # Numbers of atoms and type labels
      self.out.write(" ".join(labels)+"\n")
      self.out.write(" ".join([("%d" % Z) for Z in atnums])+"\n")

      # Positions
      for i in range(len(prop_vals)):
	 for i_at in frange(at.n):
	    try:
	       match = getattr(at,property)[i_at].stripstrings() == prop_vals[i]
	    except:
	       match = getattr(at,property)[i_at] == prop_vals[i]
	    if (match):
	       self.out.write("%f %f %f   T T T\n" % (at.pos[1,i_at], at.pos[2,i_at], at.pos[3,i_at]))

      # Velocities
      if (hasattr(at, 'velo')):
	 self.out.write("\n")
	 for i in range(len(prop_vals)):
	    for i_at in frange(at.n):
	       try:
		  match = getattr(at,property)[i_at].stripstrings() == prop_vals[i]
	       except:
		  match = getattr(at,property)[i_at] == prop_vals[i]
	       if (match):
		  self.out.write("%f %f %f\n" % (at.velo[1,i_at], at.velo[2,i_at], at.velo[3,i_at]))


AtomsWriters['vasp'] = VaspWriter
AtomsWriters['VASP'] = VaspWriter
AtomsWriters['POSCAR'] = VaspWriter

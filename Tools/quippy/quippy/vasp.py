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
from quippy.atoms import Atoms
from quippy.io import atoms_reader, AtomsReaders, AtomsWriters
from quippy.farray import fzeros,frange
import numpy as np
import re, sys

__all__ = ['VASP_POSCAR_Reader', 'VASP_OUTCAR_Reader', 'VASPWriter']

@atoms_reader('vasp')
@atoms_reader('POSCAR')
@atoms_reader('CONTCAR')
def VASP_POSCAR_Reader(poscar, species=None):
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
   if (len(comment) > 0):
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

class VaspWriter(object):
   """
   Writer for VASP POSCAR format
   """

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

      swapped_a1_a2 = False
      vol = np.dot(at.lattice[:,1],np.cross(at.lattice[:,2],at.lattice[:,3]))
      if (vol < 0.0):
	 t_a1 = at.lattice[:,1].copy()
	 at.lattice[:,1] = at.lattice[:,2]
	 at.lattice[:,2] = t_a1[:]
	 swapped_a1_a2 = True
	 sys.stderr.write("WARNING: swapped a1 and a2 to get positive scalar triple product\n")

      # Comment
      try:
	 self.out.write(at.params['VASP_Comment'])
	 if (swapped_a1_a2):
	    self.out.write(" a1 and a2 swapped relative to input to get positive volume")
	 self.out.write("\n")
      except:
	 try:
	    self.out.write(at.params['comment'])
	 except:
	    self.out.write('')
	 if (swapped_a1_a2):
	    self.out.write(" a1 and a2 swapped relative to input to get positive volume")
	 self.out.write("\n")

      # Lattice
      self.out.write("1.0\n")
      self.out.write("%.12f %.12f %.12f\n" % (at.lattice[1,1], at.lattice[2,1], at.lattice[3,1]))
      self.out.write("%.12f %.12f %.12f\n" % (at.lattice[1,2], at.lattice[2,2], at.lattice[3,2]))
      self.out.write("%.12f %.12f %.12f\n" % (at.lattice[1,3], at.lattice[2,3], at.lattice[3,3]))

      # Numbers of atoms and type labels
      self.out.write(" ".join(labels)+"\n")
      self.out.write(" ".join([("%d" % Z) for Z in atnums])+"\n")
      self.out.write("Selective Dynamics\n")
      if (hasattr(at, 'VASP_Pos_Format')):
	 self.out.write(at.params['VASP_Pos_Format']+"\n")
	 if ((at.params['VASP_Pos_Format'][0] == 'd' or at.params['VASP_Pos_Format'][0] == 'D') and swapped_a1_a2):
	    t_p1 = at.pos[1,:].copy()
	    at.pos[1,:] = at.pos[2,:]
	    at.pos[2,:] = t_p1[:]
      else:
	 self.out.write("Cartesian\n")

      # Positions
      for i in range(len(prop_vals)):
	 for i_at in frange(at.n):
	    try:
	       match = getattr(at,property)[i_at].stripstrings() == prop_vals[i]
	    except:
	       match = getattr(at,property)[i_at] == prop_vals[i]
	    if (match):
	       self.out.write("%.12f %.12f %.12f   T T T\n" % (at.pos[1,i_at], at.pos[2,i_at], at.pos[3,i_at]))

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
		  self.out.write("%.12f %.12f %.12f\n" % (at.velo[1,i_at], at.velo[2,i_at], at.velo[3,i_at]))


AtomsWriters['vasp'] = VaspWriter
AtomsWriters['VASP'] = VaspWriter
AtomsWriters['POSCAR'] = VaspWriter


@atoms_reader('OUTCAR')
def VASP_OUTCAR_Reader(outcar, species=None):
   """Read a configuration from a VASP OUTCAR file."""
   
   if (outcar == 'stdin' or outcar == '-'):
      p = sys.stdin
   else:
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
	    # print "energy(sigma->0) line"
	    # print l.split()
	    at_cur.params['Energy_sigma_to_zero'] = float(l.split()[6])
	    energy_i += 1
            yield at_cur
         if (at_cur is not None and at_i == at_cur.n):
            at_cur.set_lattice(lat_cur, False)
            for i in frange(at_cur.n):
               dr = at_cur.diff_min_image(i, at.pos[:,i])
               at_cur.pos[:,i] = at.pos[:,i] - dr


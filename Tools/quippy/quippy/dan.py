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

from quippy.atoms import Atoms, atoms_reader, AtomsReaders, AtomsWriters
from quippy.farray import frange
from types import StringType
import sys
import numpy as np

class DanWriter(object):

    def __init__(self, out, pos_field='pos', atom_type=None, value=None,
                 vector=None, bond_by_cutoff=None, post_config_command=None,
                 post_file_command=None, end_command=None):
        self.out = out
        self.n_graphs = 0
        self.atom_type = atom_type
        self.pos_field = pos_field
        self.value = value
        self.vector = vector
        self.bond_by_cutoff = bond_by_cutoff
        self.post_config_command = post_config_command
        self.post_file_command = post_file_command
        self.end_command = end_command
        self.opened = False
        if type(self.out) == type(''):
            if self.out == 'stdout':
                self.out = sys.stdout
            else:
                self.out = open(self.out, 'w')
                self.opened = True
        self.first_config = True

    def close(self):
        if self.post_file_command is not None:
	    if (type(self.post_file_command) is StringType):
	        self.post_file_command=[self.post_file_command]
            for cmd in self.post_file_command:
                self.out.write(cmd + "\n")
        if self.end_command is not None:
	    if (type(self.end_command) is StringType):
	        self.end_command=[self.end_command]
            for cmd in self.end_command:
                self.out.write(cmd + "\n")
        if self.n_graphs > 0:
            for ig in range(self.n_graphs):
                self.out.write("graph_range %d %f %f" % (ig, self.graph_min[ig], self.graph_max[ig]) + "\n")
        if self.opened:
            self.out.close()

    def write(self, at):
	if "dan_graph" in at.params:
	    graph_vals = at.params["dan_graph"].split()
	    self.n_graphs = len(graph_vals)

        if self.first_config:
	    if self.n_graphs > 0:
	       self.out.write("n_graphs %d\n" % self.n_graphs)
	       self.graph_min = [1.0e38] * self.n_graphs
	       self.graph_max = [-1.0e38] * self.n_graphs
            self.first_config=False
        self.out.write("new_configuration" + "\n")
        
        self.out.write("pbc_a 1 %f %f %f\n" % ( at.lattice[1,1], at.lattice[2,1], at.lattice[3,1]))
        self.out.write("pbc_a 2 %f %f %f\n" % ( at.lattice[1,2], at.lattice[2,2], at.lattice[3,2]))
        self.out.write("pbc_a 3 %f %f %f\n" % ( at.lattice[1,3], at.lattice[2,3], at.lattice[3,3]))

        for i_at in frange(at.n):
            if self.atom_type is not None:
                atom_type = self.atom_type
            else:
                if (hasattr(at,'z')):
                    atom_type = 'z'
                else:
                    if (hasattr(at,'species')):
                        atom_type = 'species'
                    else:
                        if (hasattr(at,'type')):
                            atom_type = 'type'
                        else:
                            raise ValueError("Can't find z, species, or type for atom type")

	    px=getattr(at,self.pos_field)[1,i_at]
	    py=getattr(at,self.pos_field)[2,i_at]
	    pz=getattr(at,self.pos_field)[3,i_at]
	    try:
	       self.out.write("atom %f %f %f %s" % (px, py, pz, getattr(at,atom_type)[i_at].stripstrings()))
	    except:
	       self.out.write("atom %f %f %f %s" % (px, py, pz, getattr(at,atom_type)[i_at]))

            if self.value is not None:
	        if (type(self.value) is StringType):
		    self.value=[self.value]
                for iv in range(len(self.value)):
                    self.out.write(" value %d %s" % (iv+1, getattr(at,self.value[iv])[i_at]))
            self.out.write("\n")

            if self.vector is not None:
                if hasattr(at,self.vector):
                    self.out.write("vector %f %f %f   %f %f %f\n" % (px, py, pz,
                                                                     getattr(at,self.vector)[1,i_at],
                                                                     getattr(at,self.vector)[2,i_at],
                                                                     getattr(at,self.vector)[3,i_at]))
                        
        if self.n_graphs > 0:
	    for ig in range(len(graph_vals)):
	       val = graph_vals[ig]
	       self.out.write("graph_value %d %f\n" % (ig, float(val)))
	       f_val = float(val)
	       if (f_val < self.graph_min[ig]):
		  self.graph_min[ig] = f_val
	       if (f_val > self.graph_max[ig]):
		  self.graph_max[ig] = f_val
		       
        if self.bond_by_cutoff:
            self.out.write("bond_by_cutoff\n")
            
        if self.post_config_command is not None:
	    if (type(self.post_config_command) is StringType):
	        self.post_config_command=[self.post_config_command]
            for cmd in self.post_config_command:
                self.out.write(cmd + "\n")

AtomsWriters['dan'] = DanWriter

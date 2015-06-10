from pylab import *
from numpy import *

import quippy
from quippy import *

import os

ip = _ip = get_ipython()

def load_atoms(self, arg):
    """Load file into an Atoms object"""    
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0].replace('-','_').replace('.','_')
    self.ex('%s = Atoms("%s"); %s.show()' % (n,f,n))

def load_atoms_list(self, arg):
    """Load file into an AtomsList object"""
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0]
    n = n.replace('-','_').replace('.','_').replace('*','').replace('?','')
    n = n.strip('_')
    self.ex('%s = AtomsList("%s"); %s.show()' % (n,f,n))

def load_atoms_reader(self, arg):
    """Load file into an AtomsReader object"""
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0]
    n = n.replace('-','_').replace('.','_').replace('*','').replace('?','')
    n = n.strip('_')
    self.ex('%s = AtomsReader("%s"); %s.show()' % (n,f,n))

ip.define_magic("at", load_atoms)
ip.define_magic("ar", load_atoms_reader)
ip.define_magic("al", load_atoms_list)

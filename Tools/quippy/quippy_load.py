from pylab import *
from numpy import *

import quippy
from quippy import *

import os
import IPython.ipapi
ip = IPython.ipapi.get()

def load_atoms(self, arg):
    """Load file into an Atoms object"""    
    ip = self.api
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0].replace('-','_').replace('.','_')
    ip.ex('%s = Atoms("%s"); %s.show()' % (n,f,n))

def load_atoms_list(self, arg):
    """Load file into an AtomsList object"""
    ip = self.api
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0]
    n = n.replace('-','_').replace('.','_').replace('*','').replace('?','')
    n = n.strip('_')
    ip.ex('%s = AtomsList("%s"); %s.show()' % (n,f,n))

def load_atoms_reader(self, arg):
    """Load file into an AtomsReader object"""
    ip = self.api
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0]
    n = n.replace('-','_').replace('.','_').replace('*','').replace('?','')
    n = n.strip('_')
    ip.ex('%s = AtomsReader("%s"); %s.show()' % (n,f,n))

ip.expose_magic("at", load_atoms)
ip.expose_magic("ar", load_atoms_reader)
ip.expose_magic("al", load_atoms_list)

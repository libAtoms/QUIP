from pylab import *
from numpy import *

import quippy
from quippy import *
from quippy.plot2d import *

import IPython.ipapi
ip = IPython.ipapi.get()

def load_atoms(self, arg):
    ip = self.api
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0].replace('-','_').replace('.','_')
    ip.ex('%s = Atoms("%s"); %s.show()' % (n,f,n))

def load_atoms_list(self, arg):
    ip = self.api
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0]
    n = n.replace('-','_').replace('.','_').replace('*','').replace('?','')
    n = n.strip('_')
    ip.ex('%s = AtomsList("%s"); %s.loadall(); %s.show(frame=0)' % (n,f,n,n))

ip.expose_magic("at", load_atoms)
ip.expose_magic("al", load_atoms_list)

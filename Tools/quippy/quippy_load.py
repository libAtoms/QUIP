# Set matplotlib in interactive mode with the macosx backend
# THESE MUST BE THE FIRST MATPLOTLIB COMMANDS CALLED!
import matplotlib
matplotlib.use('macosx')
matplotlib.interactive(True)

# Bring all of the numeric and plotting commands to the toplevel namespace
from pylab import *

from numpy import *
import quippy
from quippy import *

print """Welcome to quippy!

help(numpy) -> help on NumPy, Python's basic numerical library.
help(plotting) -> help on plotting commands.
help(quippy) -> help on quippy commands
"""

import IPython.ipapi
ip = IPython.ipapi.get()

def load_atoms(self, arg):
    ip = self.api
    f = arg
    n = os.path.splitext(os.path.basename(arg))[0].replace('-','_').replace('.','_')
    ip.ex('%s = Atoms("%s")' % (n,f))

ip.expose_magic("atoms", load_atoms)

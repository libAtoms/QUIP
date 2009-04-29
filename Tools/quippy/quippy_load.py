# Set matplotlib in interactive mode with the TkAgg backend
# THESE MUST BE THE FIRST MATPLOTLIB COMMANDS CALLED!
import matplotlib
matplotlib.use('TkAgg')
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

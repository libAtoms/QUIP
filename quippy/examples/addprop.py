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

"""Simple example of looping over an input file, adding properties
   and writing to an output file.
"""

from __future__ import print_function, unicode_literals

from quippy import AtomsList
import sys

if len(sys.argv[1:]) != 2:
   print('Usage: addprop.py <INPUT> <OUTPUT>')
   sys.exit(1)

frames = AtomsList(sys.argv[1])

for at in frames[1:]:
  at.add_property('p', 0.0, n_cols=3)
  at.p[:] = frames[0].pos[:]

(frames[1:]).write(sys.argv[2])

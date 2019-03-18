# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2019
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

import sys
import unittest
import os.path
from distutils.util import get_platform
from distutils.sysconfig import get_python_version

if 'QUIP_ROOT' in os.environ:
    quip_root = os.environ['QUIP_ROOT']
else:
    quip_root = os.path.join(os.getcwd(), '..')

try:
    quip_arch = os.environ['QUIP_ARCH']
except KeyError:
    raise RuntimeError('You need to define the architecture using the QUIP_ARCH variable. Check out the arch/ subdirectory.')
print('QUIP_ARCH', quip_arch)

# extend sys.path
sys.path.insert(0, os.path.join(quip_root, 'build', quip_arch))
print(sys.path)

# find tests and run them
# ONLY RUNS ONE NOW, THE ONE THAT IS DONE
suite = unittest.defaultTestLoader.discover(os.getcwd(), pattern='test_pot.py')  # fixme run all tests not just one
unittest.TextTestRunner().run(suite)

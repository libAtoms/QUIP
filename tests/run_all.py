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
import os
import os.path

import quippy
print('Successfully imported quippy')

# Set QUIP_WHEEL_TEST to skip shell script tests that require Fortran binaries
os.environ['QUIP_WHEEL_TEST'] = '1'

# find tests and run them
suite = unittest.defaultTestLoader.discover(os.getcwd())
result = unittest.TextTestRunner(verbosity=2).run(suite)
if result.wasSuccessful():
    sys.exit(0)
else:
    sys.exit(1)


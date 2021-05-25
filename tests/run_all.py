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

platform = '{0}-{1}'.format(get_platform(), get_python_version())
print('platform', platform)

# optionally extend sys.path
quip_test_in_place = True
if 'QUIP_TEST_IN_PLACE' in os.environ:
    quip_test_in_place = int(os.environ['QUIP_TEST_IN_PLACE']

if quip_test_inplace:
    build_dir = os.path.join(quip_root, 'build/{0}/'.format(quip_arch)
    print(f'Adding QUIP build directory {build_dir} to sys.path')
    sys.path.insert(0, build_dir)

import quippy
print('Successfully imported quippy3')

# find tests and run them
suite = unittest.defaultTestLoader.discover(os.getcwd())
result = unittest.TextTestRunner().run(suite)
if result.wasSuccessful():
    sys.exit(0)
else:
    sys.exit(1)


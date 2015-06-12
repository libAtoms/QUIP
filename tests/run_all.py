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

import sys, unittest, os.path, glob
from distutils.util import get_platform
from distutils.sysconfig import get_python_version

if 'QUIP_ROOT' in os.environ:
    quip_root = os.environ['QUIP_ROOT']
else:
    quip_root = os.path.join(os.getcwd(), '..')
quip_arch = os.environ['QUIP_ARCH']
print('QUIP_ARCH', quip_arch)

platform = '{0}-{1}'.format(get_platform(), get_python_version())
print('platform', platform)

# extend sys.path
sys.path.insert(0, os.path.join(quip_root, 'quippy', 'build/{0}/lib.{1}'.format(quip_arch,
                                                                                platform)))
# find all the tests
test_files = glob.glob(os.path.join(os.path.dirname(__file__), 'test*.py'))
test_mods = [ os.path.splitext(os.path.basename(t))[0] for t in test_files ]

for name in test_mods:
   __import__(name)
   globals().update(vars(sys.modules[name]))

unittest.main()
   

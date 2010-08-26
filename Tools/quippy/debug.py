import sys
from custom_commands import *
from numpy.distutils.core import setup

if len(sys.argv[1:]) != 1:
    print 'Usage: debug.py test_case'
    sys.exit(1)

sys.argv.insert(1, 'test')
sys.argv.insert(2, '--verbosity=2')
if not sys.argv[3].startswith('--test='):
    sys.argv[3] = '--test=' + sys.argv[3]

setup(name='quippy', cmdclass = {'test': test},
      options={'build': { 'build_base': 'build.%s' % os.environ['QUIP_ARCH']}})

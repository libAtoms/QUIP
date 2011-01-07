import sys
from custom_commands import *
from numpy.distutils.core import setup
import optparse

p = optparse.OptionParser(usage='%prog [options]')
p.add_option('-e', '--execute', action='store')
p.add_option('-t', '--test', action='store')

opt, arg = p.parse_args()

if (opt.test is None) + (opt.execute is None) != 1:
    p.error('Exactly one of -e and -t must be present.')

sys.argv[1:] = []

if opt.test is not None:
    sys.argv.insert(1, 'test')
    sys.argv.insert(2, '--verbosity=2')
    sys.argv.insert(3,'--test=' + opt.test)

if opt.execute is not None:
    sys.argv.insert(1, 'interact')
    sys.argv.insert(2, '--execute=' + opt.execute)

setup(name='quippy', cmdclass = {'test': test, 'interact' : interact},
      options={'build': { 'build_base': 'build.%s' % os.environ['QUIP_ARCH']}})

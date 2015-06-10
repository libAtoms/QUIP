from __future__ import print_function

import os
import sys
import glob
from unittest import TestLoader, TextTestRunner
from distutils.util import get_platform
from distutils.sysconfig import get_python_version

test_glob = 'test_*.py'

if 'QUIP_ROOT' in os.environ:
    quip_root = os.environ['QUIP_ROOT']
else:
    quip_root = os.path.join(os.getcwd(), '..')
quip_arch = os.environ['QUIP_ARCH']
print('QUIP_ARCH', quip_arch)

platform = '{0}-{1}'.format(get_platform(), get_python_version())
print('platform', platform)

# extend sys.path
sys.path.insert(0, os.path.join(quip_root, 'quippy', 'build.{0}/lib.{1}'.format(quip_arch,
                                                                                platform)))
print('sys.path', sys.path)

test_files = glob.glob(test_glob)
print('test_files', test_files)

test_names = [os.path.splitext(test_file)[0] for test_file in test_files ]
tests = TestLoader().loadTestsFromNames(test_names)
runner = TextTestRunner(verbosity=2)
result = runner.run(tests)

# Exit script with exitcode 1 if any tests failed
if result.failures or result.errors:
    sys.exit(1)


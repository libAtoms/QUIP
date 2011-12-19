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

import unittest, glob, os, sys, quippy
from quippy import QUIP_ROOT, QUIP_ARCH, QUIP_MAKEFILE
from quippytest import *

mpi_n_cores = [1, 2, 4]

class Test_TopLevel(QuippyTestCase):
    pass

def make_test(script, env=None):
    def run_script(self):
        if env is not None: os.environ.update(env)

        predicate = []
        if os.path.exists(script):
            predicate = [L for L in open(script).readlines() if L.startswith('#PREDICATE:')]
        if predicate != []:
            predicate = predicate[0][len('#PREDICATE:'):]
            if not eval(predicate):
                print 'Predicate %s is false, skipping test' % predicate
                return
            
        exitcode = os.system("/bin/bash %s" % script)
        if exitcode != 0:
            self.fail()

    return run_script

scripts =  [script for script in glob.glob(os.path.join(QUIP_ROOT, 'Tests/test_*')) if '.' not in script and not script.endswith('~')]

if 'mpi' in quippy.available_modules:
    for script in scripts:
        for cores in mpi_n_cores:
            setattr(Test_TopLevel, '%s_MPI_%d_cores' % (os.path.basename(script), cores),
                    make_test(script,{'QUIP_ROOT': QUIP_ROOT,
                                      'QUIP_ARCH': QUIP_ARCH,
                                      'MPIRUN': 'mpirun -n %d' % cores}))

    for cores in mpi_n_cores:
        setattr(Test_TopLevel, 'test_python_MPI_%d_cores' % cores,
                make_test('mpirun -n %d python %s/Tools/quippy/tests/test_mpi.py' % (cores, QUIP_ROOT),
                          {'QUIP_ROOT': QUIP_ROOT,
                           'QUIP_ARCH': QUIP_ARCH,
                           'PYTHONPATH': ':'.join(sys.path)}))
        

else:
    for script in scripts:
        setattr(Test_TopLevel, os.path.basename(script), make_test(script, {'QUIP_ROOT': QUIP_ROOT,
                                                                            'QUIP_ARCH': QUIP_ARCH}))



                          

if __name__ == '__main__':
   unittest.main()

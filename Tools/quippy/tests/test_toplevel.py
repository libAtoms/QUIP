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

import unittest, glob, os
from quippytest import *

if 'QUIP_ROOT' in os.environ:

    class Test_TopLevel(QuippyTestCase):
        pass


    def make_test(script):
        def run_script(self):
            exitcode = os.system(script)
            if exitcode != 0:
                self.fail()

        return run_script

    for script in glob.glob(os.path.join(os.environ['QUIP_ROOT'], 'Tests/test_*')):
        if '.' in script or script.endswith('~') : continue # skip data and backup files
        setattr(Test_TopLevel, os.path.basename(script), make_test(script))

                          

if __name__ == '__main__':
   unittest.main()

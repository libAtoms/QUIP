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

# MPI test cases - should pass both in serial and parallel

# They are run in serial as part of the normal test suite, and in
# parallel as 'mpirun -n NPROC python test_mpi.py' from
# test_toplevel.py.

from quippy import *
from quippy.atoms import *
from numpy import *

import unittest, quippy
from quippytest import *

class Test_MPI(QuippyTestCase):
    def setUp(self):
        self.mpi = MPI_context()
        #mainlog.mpi_all_inoutput(True) # for debugging

    def test_atoms_bcast(self):
        a = Atoms("""12
Lattice="3.679445 -0.181265 -0.219095 -0.140937 3.759568 0.175386 -0.223013 -0.011764 9.715448" Properties=species:S:1:pos:R:3:frac_pos:R:3:f:R:3 castep_run_time=40430.63 energy=-9949.25932697 virial="6.38880278817 4.17889707917 7.15229427215 4.17889707917 -0.740653778339 -2.86812348957 7.15229427215 -2.86812348957 -4.60133902559"
O       0.07892115      1.81915384      1.95631150      0.05184500      0.48697900      0.19373900     -0.08661000      1.65235000      0.57685000
O       0.03867388      1.98421766     -1.99187494      0.01775400      0.52796400     -0.21415200     -0.53478000     -0.09838000      1.11345000
O      -0.00808325     -0.07057365      4.35164490      0.02438200     -0.01619200      0.44875200      1.66233000      0.71894000      0.73628000
O      -0.07597036      0.06054846      0.43978436     -0.01735400      0.01540800      0.04459700      0.55106000     -0.33373000     -0.94889000
O       1.88194480      1.90132724      5.13029557      0.56414400      0.53459200      0.53112700     -0.50613000     -0.96190000      2.26224000
O       1.86840060     -0.05023179      2.83158150      0.52666200      0.01298000      0.30309400     -0.37249000      0.16802000     -0.89396000
O       1.83382247     -0.00862026     -2.78189090      0.48244900      0.02010500     -0.27582000      0.48621000      0.09282000     -1.55320000
O       1.90053677      1.80865713     -0.42277565      0.53347300      0.50667500     -0.04063200     -2.57147000     -0.07581000     -1.09444000
Ti      -0.09878726      1.99062810      0.05586401     -0.00682000      0.52914200     -0.00395600      2.77217000     -0.23594000      0.91766000
Ti       0.08835527      0.05867827      2.34945896      0.03940400      0.01826600      0.24238600     -0.43629000     -0.90315000      1.27819000
Ti       1.90982467      1.87755198     -2.42094136      0.52416000      0.52390600     -0.24682200     -0.51711000     -0.31319000     -1.29767000
Ti       1.88176829      0.00352974      4.80843526      0.54323600      0.02871600      0.50665900     -0.44688000      0.28996000     -1.09650000""", format='string')

        b = Atoms()
        if self.mpi.my_proc == 0:
            b = a.copy()
        bcast(self.mpi, b)
        self.assertEqual(a, b)

    def test_parallel_read(self):
        s="""6
EMP2=2.67539000e-03  ELONG=1.46920613e-03  COM=3.16174712e+00 Lattice="20 0 0 0 20 0 0 0 20" Properties=species:S:1:pos:R:3
O	8.433042	7.254127	4.000330
H	9.329070	6.902494	4.019458
H	8.505653	8.201400	4.069496
O	5.939224	8.979049	3.011011
H	5.602166	8.307183	3.670115
H	6.859570	9.215634	3.262673"""
        cio = CInOutput(mpi=self.mpi)
        a = Atoms(s, format='string')
        b = cio.read(str=s)
        self.assertEqual(a,b)
        


if __name__ == '__main__':
   unittest.main()
    
         

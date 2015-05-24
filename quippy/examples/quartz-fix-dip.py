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

from quippy import *
verbosity_set_minimum(VERBOSE)

at = alpha_quartz(**sio2.quartz_params['ASAP_JRK'])

p = Potential('IP TS', """
<TS_params betapol="0.75" cutoff="20.0 20.0 18.0 0.0" cutoff_coulomb="20.0" cutoff_ms="18.0" tolpol="0.0005" yuksmoothlength="10.0" iesr="-1 -1 -1" a_ew="1e-06" n_types="2" gcut="0.0" pred_order="2" maxipol="60" raggio="0.0" tewald="F" yukalpha="0.1"><per_type_data atomic_num="8" pol="14.131863" z="-1.4295594" type="1"></per_type_data><per_pair_data C_pol="0.44302622" atnum_j="8" atnum_i="8" D_ms="0.00030700577" gamma_ms="12.165654" B_pol="1.1221903" R_ms="7.0252019"></per_pair_data><per_type_data atomic_num="14" pol="0.0" z="2.8591188" type="2"></per_type_data><per_pair_data C_pol="-1.5003213" atnum_j="8" atnum_i="14" D_ms="0.0020129372" gamma_ms="11.350477" B_pol="1.973181" R_ms="4.5780828"></per_pair_data><per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="0.33967532" gamma_ms="-0.17694797" B_pol="0.0" R_ms="-0.085202834"></per_pair_data></TS_params>
""")

at.set_cutoff(p.cutoff())
at.calc_connect()


print_title('With polarisation')
p.calc(at, force=True)
f1 = at.force.copy()

print at.dipoles

at.fixdip[:] = True

print_title('Without polarisation')
p.calc(at, force=True)
f2 = at.force.copy()

print abs(f1 - f2).max()

at.pos += random.uniform(size=3*at.n).reshape((3,at.n),order='F')

p.minim(at, 'cg', 1e-3, 100, do_pos=True, do_lat=True)

#from quippy.elasticity import *
#C = elastic_constants(mp, at, 'trigonal_low')



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
from quippy.elasticity import *

def bonds_angles(aq):

   aq.calc_connect()

   lengths = []
   angles = []

   for i in frange(aq.n):
      #if aq.z[i] != 14: continue

      nn = aq.neighbours[i]
      r = [n.distance for n in nn]
      diff = farray([n.diff for n in nn]).T
      cosine = farray([n.cosines for n in nn]).T

      print i, aq.z[i], r
      lengths.extend(r)

      for j in frange(len(nn)):
         for k in frange(len(nn)):
            if j == k: continue
            theta = arccos(dot(cosine[:,j], cosine[:,k]))*180./PI
            angles.append(theta)
            print j, k, theta

   return lengths, angles



def tetrahedron(at, i):
   if at.z[i] != 14: return

   nn = at.neighbours[i]
   if len(nn) != 4: return

   ndiff = farray([n.diff for n in nn]).T

   n1 = ndiff[:,2] - ndiff[:,1]
   n2 = ndiff[:,3] - ndiff[:,1]
   n3 = ndiff[:,4] - ndiff[:,1]

   return n1, n2, n3

aq = Atoms("""9
Lattice="2.420218 -4.191941 -0.000000 2.420218 4.191941 0.000000 -0.000000 0.000000 5.328546" Properties=species:S:1:pos:R:3:Z:I:1
Si              1.12341729     -1.94581532     -1.77618189 14
Si              1.12341701      1.94581566      1.77618193 14
Si             -2.24683410     -0.00000018     -0.00000029 14
O               1.67107820     -0.55763243     -1.19142205 8
O              -0.35261542      1.72601223      2.36094217 8
O              -1.31846300     -1.16837972      0.58475994 8 
O               1.67107818      0.55763217      1.19142214 8
O              -1.31846298      1.16837979     -0.58475995 8
O              -0.35261516     -1.72601222     -2.36094201 8""".split('\n'), format='pupyxyz')


# Beta quartz, from NRL <http://cst-www.nrl.navy.mil/lattice/struk/sio2b.html>
bq = Atoms("""9
Lattice="2.49825000 -4.32709593  0.00000000 2.49825000  4.32709593  0.00000000 0.00000000  0.00000000  5.45700000" Properties=species:S:1:pos:R:3:Z:I:1
Si    1.24912500 -2.16354796  0.00000000 14
Si    1.24912500  2.16354796  3.63800002 14
Si    2.49825000  0.00000000  1.81899998 14
O     1.55591010  0.89830512  2.72850000 8
O    -1.55591010  0.89830512  0.90950002 8
O     0.00000000 -1.79661023  4.54749998 8
O    -1.55591010 -0.89830512  2.72850000 8
O     1.55591010 -0.89830512  0.90950002 8
O     0.00000000  1.79661023  4.54749998 8""".split('\n'), format='pupyxyz')

xml_lda_500k = """
<TS_params betapol="0.75" cutoff="20.0 20.0 18.0 0.0" cutoff_coulomb="20.0" cutoff_ms="18.0" tolpol="0.0005" yuksmoothlength="10.0" iesr="-1 -1 -1" a_ew="1e-06" n_types="2" gcut="0.0" pred_order="2" maxipol="60" raggio="0.0" tewald="F" yukalpha="0.1"><per_type_data atomic_num="8" pol="14.131863" z="-1.4295594" type="1"></per_type_data><per_pair_data C_pol="0.44302622" atnum_j="8" atnum_i="8" D_ms="0.00030700577" gamma_ms="12.165654" B_pol="1.1221903" R_ms="7.0252019"></per_pair_data><per_type_data atomic_num="14" pol="0.0" z="2.8591188" type="2"></per_type_data><per_pair_data C_pol="-1.5003213" atnum_j="8" atnum_i="14" D_ms="0.0020129372" gamma_ms="11.350477" B_pol="1.973181" R_ms="4.5780828"></per_pair_data><per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="0.33967532" gamma_ms="-0.17694797" B_pol="0.0" R_ms="-0.085202834"></per_pair_data></TS_params>
"""

def add_asap_props(at):
   at.add_property('efield', 0.0, n_cols=3)
   at.add_property('dipoles', 0.0, n_cols=3)
   at.add_property('efield_old1', 0.0, n_cols=3)
   at.add_property('efield_old2', 0.0, n_cols=3)
   at.add_property('efield_old3', 0.0, n_cols=3)


if True:
   pot = Potential('IP TS', xml_lda_500k)
   aq.set_cutoff(pot.cutoff())
   aq.calc_connect()

   aq.add_property('efield', 0.0, n_cols=3)
   aq.add_property('dipoles', 0.0, n_cols=3)
   aq.add_property('efield_old1', 0.0, n_cols=3)
   aq.add_property('efield_old2', 0.0, n_cols=3)
   aq.add_property('efield_old3', 0.0, n_cols=3)

   # Ensure we're at the equilibrium geometry
   pot.minim(aq, 'cg', 1e-10, 100, do_pos=True, do_lat=True)

   # Calculate stress, should be very close to zero
   pot.calc(aq, virial=True)
   aq.params['stress'] = -aq.virial*GPA/aq.cell_volume()

   strained_configs = generate_strained_configs(aq, 'trigonal_low')
   stressed_configs = calc_stress(strained_configs, pot, relax=True)
   C, C_err = fit_elastic_constants(stressed_configs, 'trigonal_low', verbose=True, graphics=True)

   bond_length = (1.6067 + 1.6027)/2
   a = bond_length/(sqrt(3.0)/4.0)


   aqc = alpha_quartz_cubic(**sio2.get_quartz_params(aq))

   aqc.calc_connect()

   elastic_fields(aqc, a=a, cij=C, save_reference=True)

   eps = strain_matrix([-1e-3,0,0,0,0,0])
   aqc2 = transform(aqc, eps)

   elastic_fields(aqc2, a=a, cij=C, use_reference=True)

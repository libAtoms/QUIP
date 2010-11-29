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
from quippy.elastic import *

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


def elastic_fields(at, a, cij, save_reference=False, use_reference=False):
   """
   Compute atomistic strain field and linear elastic stress response.

   Stress and strain are stored in compressed Voigt notation::

      at.strain[:,i] = [e_xx,e_yy,e_zz,e_yz,e_xz,e_xy]
      at.stress[:,i] = [sig_xx, sig_yy, sig_zz, sig_yz, sig_xz, sig_xy]

   so that sig = dot(C, strain) in the appropriate reference frame.
   
   Four-fold coordinate atoms within `at` are used to define
   tetrahedra. The deformation of each tetrahedra is determined
   relative to the ideal structure, using `a` as the cubic lattice
   constant (related to bond length by a factor :math:`sqrt{3}/4`).
   This deformation is then split into a strain and a rotation
   using a Polar decomposition.

   If `save_reference` or `use_reference` are True then `at` must have
   a `primitive_index` integer property which is different for each
   atom in the primitive unit cell. `save_reference` causes the
   local strain and rotation for one atom of each primitive type
   to be saved as entries `at.params`. Conversely, `use_reference` uses
   this information and to undo the local strain and rotation.

   The strain is then transformed into the crystal reference frame
   (i.e. x=100, y=010, z=001) to calculate the stress using the `cij`
   matrix of elastic constants. Finally the resulting stress is
   transformed back into the sample frame.

   The stress and strain tensor fields are interpolated to give
   values for atoms which are not four-fold coordinated (for example
   oxygen atoms in silica).
   
   Eigenvalues and eigenvectors of the stress are stored in the
   properties `stress_evals`,`stress_evec1`, `stress_evec2` and
   `stress_evec3`, ordered decreasingly by eigenvalue so that the
   principal eigenvalue and eigenvector are `stress_evals[1,:]` and
   `stress_evec1[:,i]` respectively.
   
   """

   def compute_stress_eig(i):
      D, SigEvecs = linalg.eig(stress_matrix(at.stress[:,i]))

      # Order by descending size of eigenvalues
      sorted_evals, order = zip(*sorted(zip(D, [1,2,3]),reverse=True))
         
      at.stress_eval[:,i]  = sorted_evals
      at.stress_evec1[:,i] = SigEvecs[:,order[0]]
      at.stress_evec2[:,i] = SigEvecs[:,order[1]]
      at.stress_evec3[:,i] = SigEvecs[:,order[2]]


   # We want nearest neighbour connectivity only
   save_cutoff, save_use_uniform_cutoff = at.cutoff, at.use_uniform_cutoff
   at.set_cutoff_factor(1.2)
   at.calc_connect()

   if (save_reference or use_reference) and not at.has_property('primitive_index'):
      raise ValueError('Property "primitive_index" missing from Atoms object')

   # Add various properties to store results
   at.add_property('strain', 0.0, n_cols=6)
   at.add_property('stress', 0.0, n_cols=6)
   at.add_property('stress_eval', 0.0, n_cols=3)
   at.add_property('stress_evec1', 0.0, n_cols=3)
   at.add_property('stress_evec2', 0.0, n_cols=3)
   at.add_property('stress_evec3', 0.0, n_cols=3)

   rotXYZ = fzeros((3,3))
   rotXYZ[1,2] = 1.0
   rotXYZ[2,3] = 1.0
   rotXYZ[3,1] = 1.0
   E = fidentity(3)
   

   # Consider first atoms with four neighbours
   for i in frange(at.n):
      neighb = at.neighbours[i]

      if len(neighb) == 4:
         #print '\nProcessing atom', i

         # Consider neighbours in order of their index within the primitive cell
         (j1,i1), (j2,i2), (j3,i3), (j4,i4) = sorted((at.primitive_index[n.j],i) for i,n in fenumerate(neighb))

         # Find cubic axes from neighbours
         n1 = neighb[i2].diff - neighb[i1].diff
         n2 = neighb[i3].diff - neighb[i1].diff
         n3 = neighb[i4].diff - neighb[i1].diff
         #print 'n1', n1.norm(), n1
         #print 'n2', n2.norm(), n2
         #print 'n3', n3.norm(), n3

         E[:,1] = (n1 + n2 - n3)/a
         E[:,2] = (n2 + n3 - n1)/a
         E[:,3] = (n3 + n1 - n2)/a

         #print 'E', E

         # Kill near zero elements
         E[abs(E) < 1e-6] = 0.0

         if (E < 0.0).all(): E = -E

         # Find polar decomposition: E = S*R where S is symmetric, 
         # and R is a rotation
         #
         #  EEt = E*E', EEt = VDV' D diagonal, S = V D^1/2 V', R = S^-1*E
         EEt = dot(E, E.T)

         D, V = linalg.eig(EEt)

         S = farray(dot(dot(V, diag(sqrt(D))), V.T))
         R = farray(dot(dot(dot(V, diag(D**-0.5)), V.T), E))

         #print 'S:', S
         #print 'R:', R

         if save_reference:
            key = 'strain_inv_%d' % at.primitive_index[i]
            if key not in at.params:
               at.params[key] = linalg.inv(S)

         if use_reference:
            S = dot(S, at.params['strain_inv_%d' % at.primitive_index[i]])
            
         #print 'S after apply S0_inv:', S

         # Strain in rotated coordinate system
         at.strain[:,i] = strain_vector(S)

         # Test for permutations - check which way x points
         RtE = dot(R.T, E)
         if RtE[2,1] > RtE[1,1] and RtE[2,1] > RtE[3,1]:
            R = dot(rotXYZ, R)
         elif RtE[3,1] > RtE[1,1] and RtE[3,1] > RtE[2,1]:
            R = dot(rotXYZ.T, R)

         if save_reference:
            key = 'rotation_inv_%d' % at.primitive_index[i]
            if key not in at.params:
               at.params[key] = R.T
               
         if use_reference:
            R = dot(R, at.params['rotation_inv_%d' % at.primitive_index[i]])

         print 'R after apply R0t:', R            

         # Rotate to crystal coordinate system to apply Cij matrix
         RtSR = dot(dot(R.T, S), R)
         print 'RtSR:', RtSR

         sig = stress_matrix(dot(cij, strain_vector(RtSR)))

         print 'sig', sig

         # Rotate back to local coordinate system
         RsigRt = dot(dot(R, sig), R.T)

         # Symmetrise stress tensor
         RsigRt = (RsigRt + RsigRt.T)/2.0
         at.stress[:,i] = stress_vector(RsigRt)

         compute_stress_eig(i)

   # For atoms without 4 neighbours, interpolate stress and strain fields
   for i in frange(at.n):
      neighb = at.neighbours[i]
      if len(neighb) == 4: continue

      at.strain[:,i] = at.strain[:,[n.j for n in neighb]].mean(axis=1)
      at.stress[:,i] = at.stress[:,[n.j for n in neighb]].mean(axis=1)

      compute_stress_eig(i)

   # Restore original neighbour cutoff
   at.cutoff, at.use_uniform_cutoff = save_cutoff, save_use_uniform_cutoff

   

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
<ASAP_params betapol="0.75" cutoff="20.0 20.0 18.0 0.0" cutoff_coulomb="20.0" cutoff_ms="18.0" tolpol="0.0005" yuksmoothlength="10.0" iesr="-1 -1 -1" a_ew="1e-06" n_types="2" gcut="0.0" pred_order="2" maxipol="60" raggio="0.0" tewald="F" yukalpha="0.1"><per_type_data atomic_num="8" pol="14.131863" z="-1.4295594" type="1"></per_type_data><per_pair_data C_pol="0.44302622" atnum_j="8" atnum_i="8" D_ms="0.00030700577" gamma_ms="12.165654" B_pol="1.1221903" R_ms="7.0252019"></per_pair_data><per_type_data atomic_num="14" pol="0.0" z="2.8591188" type="2"></per_type_data><per_pair_data C_pol="-1.5003213" atnum_j="8" atnum_i="14" D_ms="0.0020129372" gamma_ms="11.350477" B_pol="1.973181" R_ms="4.5780828"></per_pair_data><per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="0.33967532" gamma_ms="-0.17694797" B_pol="0.0" R_ms="-0.085202834"></per_pair_data></ASAP_params>
"""

def add_asap_props(at):
   at.add_property('efield', 0.0, n_cols=3)
   at.add_property('dipoles', 0.0, n_cols=3)
   at.add_property('efield_old1', 0.0, n_cols=3)
   at.add_property('efield_old2', 0.0, n_cols=3)
   at.add_property('efield_old3', 0.0, n_cols=3)


if True:
   pot = Potential('IP ASAP2', xml_lda_500k)
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

   aqc.set_cutoff_factor(1.2)
   aqc.calc_connect()

   elastic_fields(aqc, a=a, cij=C, save_reference=True)

   eps = strain_matrix([-1e-3,0,0,0,0,0])
   aqc2 = transform(aqc, eps)

   elastic_fields(aqc2, a=a, cij=C, use_reference=True)

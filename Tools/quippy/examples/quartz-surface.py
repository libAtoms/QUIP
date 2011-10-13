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
from quippy.surface import *
from quippy.structure_tools import *

def add_asap_props(at):
   at.add_property('efield', 0.0, n_cols=3)
   at.add_property('dipoles', 0.0, n_cols=3)
   at.add_property('efield_old1', 0.0, n_cols=3)
   at.add_property('efield_old2', 0.0, n_cols=3)
   at.add_property('efield_old3', 0.0, n_cols=3)



# Screened ASAP quartz parameters
quartz_params = {'a': 4.8403809707320216,
                 'c': 5.3285240037002248,
                 'u': 0.46417561617105912,
                 'x': 0.41174271054205958,
                 'y': 0.27872745399831672,
                 'z': 0.10973603276909905}

# CASTEP LDA quartz parameters
#quartz_params = sio2.quartz_params['CASTEP_LDA']

aq = alpha_quartz(**quartz_params)

nx = 1
ny = 3
nz = 1
d  = 15.0
h = 15.0

quartz_bulk = {}
quartz_surface = {}

# (0001)[-110] and (0001)[-100]
# (10-10)[001] and (10-10)[010]

indices = [ ((0,0,0,1), (-1,1,0), None),
            ((0,0,0,1), (-1,0,0), None),
            ((1,0,-1,0), (0,0,1), [0.0, 0.07, 0.0]),
            ((1,0,-1,0), (0,1,0), [0.0, 0.07, 0.0])]

for (y, z, shift) in indices:
   quartz_bulk[(y,z)] = orthorhombic_slab(aq, rot=rotation_matrix(aq, y, z), shift=shift)

   quartz_surface[(y,z)]  = supercell(quartz_bulk[(y,z)], nx, ny, nz)
   quartz_surface[(y,z)].lattice[2,2] += d
   quartz_surface[(y,z)].set_lattice(quartz_surface[(y,z)].lattice, False)

   quartz_surface[(y,z)].params['axes'] = rotation_matrix(aq, y, z).T

# (10-11)[010]

rot = rotation_matrix(aq, y=[1,0,-1,1], z=[0,1,0])
quartz_surface[((1,0,-1,1),(0,1,0))]  = orthorhombic_slab(aq, rot=rot, periodicity=[0.0, h, 0.0],
                                                          shift=[0.0, 0.0, 0.0], vacuum=[0.0, d, 0.0], verbose=False)
quartz_surface[((1,0,-1,1),(0,1,0))].params['axes'] = rot.T

# (10-11)[21-2]

rot = rotation_matrix(aq, y=[1,0,-1,1], z=[2,1,-2])
quartz_surface[((1,0,-1,1),(2,1,-2))]  = orthorhombic_slab(aq, rot=rot, periodicity=[0.0, h, 0.0],
                                                          shift=[0.0, 0.0, 0.0], vacuum=[0.0, d, 0.0], verbose=False)
quartz_surface[((1,0,-1,1),(2,1,-2))].params['axes'] = rot.T


# (10-1-1)[010]

rot = rotation_matrix(aq, y=[1,0,-1,-1], z=[0,1,0])
quartz_surface[((1,0,-1,-1),(0,1,0))] = orthorhombic_slab(aq, rot=rot, periodicity=[0.0, h + 2.0, 0.0], vacuum=[0.0, d, 0.0],
                                                          shift=[0.0, 0.0, 0.0], verbose=False)
quartz_surface[((1,0,-1,-1),(0,1,0))].params['axes'] = rot.T

# (10-1-1)[212]

rot = rotation_matrix(aq, y=[1,0,-1,-1], z=[2,1,2])
quartz_surface[((1,0,-1,-1),(2,1,2))] = orthorhombic_slab(aq, rot=rot, periodicity=[0.0, h + 2.0, 0.0], vacuum=[0.0, d, 0.0],
                                                          shift=[0.0, 0.0, 0.0], verbose=False)
quartz_surface[((1,0,-1,-1),(2,1,2))].params['axes'] = rot.T


xml_ewald="""<?xml version="1.0" encoding="iso-8859-1"?>
<TS_params betapol="0.75" cutoff="18.0 0.0 0.0 0.0" tolpol="0.0005" yuksmoothlength="0.0" iesr="-1 -1 -1" a_ew="1e-06" n_types="2" gcut="0.0" pred_order="2" maxipol="60" raggio="0.0" tewald="T" yukalpha="0.0"><per_type_data atomic_num="8" pol="8.89378" z="-1.38257" type="1"></per_type_data><per_pair_data C_pol="0.0" atnum_j="8" atnum_i="8" D_ms="0.00024748" gamma_ms="12.07092" B_pol="0.0" R_ms="7.17005"></per_pair_data><per_type_data atomic_num="14" pol="0.0" z="2.76514" type="2"></per_type_data><per_pair_data C_pol="-1.50435" atnum_j="8" atnum_i="14" D_ms="0.0019033" gamma_ms="11.1523" B_pol="2.02989" R_ms="4.6371"></per_pair_data><per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="-0.0020846" gamma_ms="10.45517" B_pol="0.0" R_ms="5.75038"></per_pair_data></TS_params>"""

xml_lda_500k = """
<TS_params betapol="0.75" cutoff="20.0 20.0 18.0 0.0" cutoff_coulomb="20.0" cutoff_ms="18.0" tolpol="0.0005" yuksmoothlength="10.0" iesr="-1 -1 -1" a_ew="1e-06" n_types="2" gcut="0.0" pred_order="2" maxipol="60" raggio="0.0" tewald="F" yukalpha="0.1"><per_type_data atomic_num="8" pol="14.131863" z="-1.4295594" type="1"></per_type_data><per_pair_data C_pol="0.44302622" atnum_j="8" atnum_i="8" D_ms="0.00030700577" gamma_ms="12.165654" B_pol="1.1221903" R_ms="7.0252019"></per_pair_data><per_type_data atomic_num="14" pol="0.0" z="2.8591188" type="2"></per_type_data><per_pair_data C_pol="-1.5003213" atnum_j="8" atnum_i="14" D_ms="0.0020129372" gamma_ms="11.350477" B_pol="1.973181" R_ms="4.5780828"></per_pair_data><per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="0.33967532" gamma_ms="-0.17694797" B_pol="0.0" R_ms="-0.085202834"></per_pair_data></TS_params>
"""


pot = Potential('IP TS', param_str=xml_lda_500k)
#pot = Potential('IP ASAP', xml_ewald)

gamma = {}
gamma_relaxed = {}

add_asap_props(aq)
aq.set_cutoff(pot.cutoff()+2.0)


do_relax = False
do_md = False

movie = CInOutput('relax.xyz', OUTPUT)


C = FortranArray([[  8.31662343e+01,   1.30178406e+01,   1.44884078e+01,
                -1.67795141e+01,   2.78491006e-05,   0.00000000e+00],
              [  1.30178406e+01,   8.31662343e+01,   1.44884078e+01,
                 1.67795141e+01,  -2.78491006e-05,   0.00000000e+00],
              [  1.44884078e+01,   1.44884078e+01,   1.18895940e+02,
                 0.00000000e+00,   0.00000000e+00,   0.00000000e+00],
              [ -1.67795141e+01,   1.67795141e+01,   0.00000000e+00,
                 5.91792377e+01,   0.00000000e+00,  -2.78491006e-05],
              [  2.78491006e-05,  -2.78491006e-05,   0.00000000e+00,
                 0.00000000e+00,   5.91792377e+01,  -1.67795141e+01],
              [  0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
                -2.78491006e-05,  -1.67795141e+01,   3.50741969e+01]]).T


for k in quartz_surface.keys():
   at = quartz_surface[k]
   at.params['YoungsModulus'] = youngs_modulus(C, at.axes[:,2])
   at.params['PoissonRatio_yx'] = poisson_ratio(C, at.axes[:,2], at.axes[:,1])
   at.params['PoissonRatio_yz'] = poisson_ratio(C, at.axes[:,2], at.axes[:,3])
   if k in quartz_bulk:
      quartz_bulk[k].params.update(at.params)
   

for index in sorted(quartz_surface.keys()):
   print index

   surface = quartz_surface[index].copy()
   add_asap_props(surface)

   surface.set_cutoff(pot.cutoff()+2.0)   
   gamma[index] = surface_energy(pot, aq, surface)*J_PER_M2
   print 'gamma[%s] = %f' % (index, gamma[index])
   surface.write('quartz_%d%d%d%d_%d%d%d.xyz' % tuple(list(index[0]) + list(index[1])))

   if do_relax:
      if do_md:
         ds = DynamicalSystem(surface)
         traj = AtomsList(ds.run(pot, dt=0.5, n_steps=1000, connect_interval=10))
         traj.loadall()
         surface = traj[-1].copy()

      pot.minim(surface, 'cg', 1e-2, 1000, do_lat=False, do_pos=True, do_print=True, print_cinoutput=movie)
      gamma_relaxed[index] = surface_energy(pot, aq, surface)*J_PER_M2
      print 'gamma_relaxed[%s] = %f' % (index, gamma_relaxed[index])


movie.close()

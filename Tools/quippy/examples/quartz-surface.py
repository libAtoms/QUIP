from quippy import *

def surface_energy(pot, bulk, surface, dir=2):
   bulk.calc_connect()
   pot.calc(bulk, calc_energy=True, calc_force=True)

   surface.calc_connect()
   pot.calc(surface, calc_energy=True, calc_force=True)

   area = surface.cell_volume()/surface.lattice[dir,dir]

   return (surface.energy - bulk.energy)/(2.0*area)
   

def crack_axes_from_miller_indices(surf, line):

   lattice = fzeros((3,3))
   lattice[:,1] = (0.5, -0.5*sqrt(3.0), 0.0)
   lattice[:,2] = (0.5,  0.5*sqrt(3.0), 0.0)
   lattice[:,3] = (0.0,    0.0,         1.0)
   
   if len(surf) == 4:
      # remove redundant rhombehedral index
      h, k, i, l = surf
      surf = dot(lattice, (h, k, l))

   if len(line) == 4:
      # remove redundant rhombehedral index
      h, k, i, l = line
      line = dot(lattice, (h, k, l))

   axes = fzeros((3,3))

   axes[:,2] = surf
   axes[:,3] = line

   # Right handed coordinate system, x = y .cross. z
   axes[:,1] = cross(surf, line)

   return axes
      


axes = farray([[1.0, 0.0, 0.0],
               [0.0, 0.0, 1.0],
               [0.0, 1.0, 0.0]])


quartz_params = {'a': 4.8403809707320216,
                 'c': 5.3285240037002248,
                 'u': 0.46417561617105912,
                 'x': 0.41174271054205958,
                 'y': 0.27872745399831672,
                 'z': 0.10973603276909905}

J_PER_M2 = ELEM_CHARGE*1.0e20

xml_ewald="""<?xml version="1.0" encoding="iso-8859-1"?>
<ASAP_params betapol="0.75" cutoff="18.0 0.0 0.0 0.0" tolpol="0.0005" yuksmoothlength="0.0" iesr="-1 -1 -1" a_ew="1e-06" n_types="2" gcut="0.0" pred_order="2" maxipol="60" raggio="0.0" tewald="T" yukalpha="0.0"><per_type_data atomic_num="8" pol="8.89378" z="-1.38257" type="1"></per_type_data><per_pair_data C_pol="0.0" atnum_j="8" atnum_i="8" D_ms="0.00024748" gamma_ms="12.07092" B_pol="0.0" R_ms="7.17005"></per_pair_data><per_type_data atomic_num="14" pol="0.0" z="2.76514" type="2"></per_type_data><per_pair_data C_pol="-1.50435" atnum_j="8" atnum_i="14" D_ms="0.0019033" gamma_ms="11.1523" B_pol="2.02989" R_ms="4.6371"></per_pair_data><per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="-0.0020846" gamma_ms="10.45517" B_pol="0.0" R_ms="5.75038"></per_pair_data></ASAP_params>"""

xml_lda_500k = """
<ASAP_params betapol="0.75" cutoff="20.0 20.0 18.0 0.0" cutoff_coulomb="20.0" cutoff_ms="18.0" tolpol="0.0005" yuksmoothlength="10.0" iesr="-1 -1 -1" a_ew="1e-06" n_types="2" gcut="0.0" pred_order="2" maxipol="60" raggio="0.0" tewald="F" yukalpha="0.1"><per_type_data atomic_num="8" pol="14.131863" z="-1.4295594" type="1"></per_type_data><per_pair_data C_pol="0.44302622" atnum_j="8" atnum_i="8" D_ms="0.00030700577" gamma_ms="12.165654" B_pol="1.1221903" R_ms="7.0252019"></per_pair_data><per_type_data atomic_num="14" pol="0.0" z="2.8591188" type="2"></per_type_data><per_pair_data C_pol="-1.5003213" atnum_j="8" atnum_i="14" D_ms="0.0020129372" gamma_ms="11.350477" B_pol="1.973181" R_ms="4.5780828"></per_pair_data><per_pair_data C_pol="0.0" atnum_j="14" atnum_i="14" D_ms="0.33967532" gamma_ms="-0.17694797" B_pol="0.0" R_ms="-0.085202834"></per_pair_data></ASAP_params>
"""

nx = 1  # repeats in x
ny = 3
nz = 1  # repeats in z
d  = 15.0 # Angstrom

pot = Potential('IP ASAP2', xml_lda_500k)
#pot = Potential('IP ASAP', xml_ewald)

gamma = []
gamma_relaxed = []

movie = CInOutput('relax.xyz', OUTPUT)

verbosity_set_minimum(NORMAL)

mp = MetaPotential('Simple', pot)

bulk = slab_nx_ny_nz(axes, nx=nx, ny=ny, nz=nz, lat_type="alpha_quartz", **quartz_params)

bulk.add_property('efield', 0.0, n_cols=3)
bulk.add_property('dipoles', 0.0, n_cols=3)
bulk.add_property('efield_old1', 0.0, n_cols=3)
bulk.add_property('efield_old2', 0.0, n_cols=3)
bulk.add_property('efield_old3', 0.0, n_cols=3)

bulk.set_cutoff(pot.cutoff()+2.0)
bulk.calc_connect()

surface = bulk.copy()
surface.calc_connect()

surface.lattice[2,2] += d
surface.set_lattice(surface.lattice)

gamma.append(surface_energy(mp, bulk, surface)*J_PER_M2)

ds = DynamicalSystem(surface)

traj = AtomsList(ds.run(pot, dt=0.5, n_steps=1000, connect_interval=10))
traj.loadall()

surface = traj[-1].copy()

mp.minim(surface, 'cg', 1e-2, 1000, do_lat=False, do_pos=True, do_print=True, print_cinoutput=movie)

gamma_relaxed.append(surface_energy(mp, bulk, surface)*J_PER_M2)

movie.close()


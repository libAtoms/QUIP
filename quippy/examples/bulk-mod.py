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

"""Fit Birch equation of state to a set of energies at different cell volumes.

   scipy is required for the least-squares non-linear optimisation routine.

   Input: set of Atoms configuations with energy and pressure in params dict

   castep_geom_loader() reads last frame from a set of CASTEP geometry optimisations.
   These should be carried out with variable ion positions and cell parameters
   at a series of fixed pressures.
   
   sample_input() routine uses SW potential to calcualte energy and virial.
   This is not quite correct since the the stress is not constrained to be hydrostatic,
   but works okay for small deviations from the equilibrium volume
   """

from pylab import *
from quippy import *
from scipy.optimize import leastsq

do_higher_order = False

def castep_geom_loader():
   a_ref = sio2.alpha_quartz()
   prange = ('-1', '-0.4', '-0.2', '0', '0.2', '0.4', '1')
   for pressure in prange:
      a = AtomsList('alpha_quartz.%s/alpha_quartz.%s.geom' % (pressure, pressure), atoms_ref=a_ref)[-1]
      a.params['pressure'] = float(pressure)
      yield a

def sample_input():
   at0 = diamond(5.43, 14)
   eps = 1e-3

   xml="""
   <SW_params n_types="2" label="PRB_31_plus_H">
   <comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984), extended for other elements </comment>
   <per_type_data type="1" atomic_num="1" />
   <per_type_data type="2" atomic_num="14" />
   <per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
         p="0" q="0" a="1.0" sigma="1.0" eps="0.0" />
   <per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827"
         p="4" q="0" a="1.25" sigma="2.537884" eps="2.1672" />
   <per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
         p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

   <!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
   <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />

   <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   <per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
         lambda="21.0" gamma="1.20" eps="2.1675" />
   </SW_params>
   """

   pot = Potential('IP SW', param_str=xml)
  
   for i in frange(-5,5):
      at = at0.copy()
      at.set_lattice(at.lattice*(1.0 + eps*i), True)
      at.set_cutoff(pot.cutoff()+1.0)
      at.calc_connect()
      pot.minim(at, 'cg', 1e-7, 100, do_pos=True, do_lat=False)
      pot.calc(at, energy=True, virial=True)
      at.pressure = at.virial.trace()*at.cell_volume()/(3.0*EV_A3_IN_GPA)
      print at.virial
      yield at
      

def birch_pressure(vo,bo,bop,v):
   return 3.0*bo/2.0*((vo/v)**(7.0/3.0) - (vo/v)**(5.0/3.0))*(1.0 + 0.75*(bop-4.0)*((vo/v)**(2.0/3.0)-1))

def birch_energy(vo,eo,bo,bop,v):
   t = (vo/v)**.6666666666666667 - 1.0
   e = eo + 1.125*bo*vo*t*t*(1.0+0.5*(bop-4.0)*t)
   return e

def birch_energy_2(vo,eo,bo,bop,b4,b5,v):
   t = (vo/v)**.6666666666666667 - 1.0
   e2 = eo + 1.125*bo*vo*t*t*(1.0 + .5*t*(bop-4.0 + t*(b4 + t*b5)))
   return e2

def errfunc(p, volume, energy):
   return birch_energy(*p, v=volume) - energy

def errfunc_2(p,volume, energy):
   vo, eo, bo, bop, b4, b5 = p
   return birch_energy_2(vo, eo, bo, bop, b4, b5, volume) - energy


##al = AtomsList(castep_geom_loader())
al = AtomsList(sample_input())

for at in al:
   print at.pressure, at.cell_volume(), at.energy

cellv = farray([at.cell_volume() for at in al])
energy = farray([at.energy for at in al])
pressure = farray([at.pressure for at in al])

# Fit vo, eo, bo, bop
(vo,eo,bo,bop), success = leastsq(errfunc, (cellv[energy.argmin()],energy.min(), 30.0, 4.0),
                                  args=(cellv, energy))

assert success > 0

print 'Volume vo =', vo, 'A^3'
print 'Energy eo =', eo, 'eV'
print 'Bulk modulus bo = ', bo, 'eV/A^3 =', bo*EV_A3_IN_GPA, 'GPa'
print 'dP/dV (P=0) bop = ', bop

if do_higher_order:
   # Higher order fit, varying vo, eo, bo, bop, b4, b5
   p2 = (cellv[energy.argmin()],energy.min(), 30.0, 4.0, 200.0, 800.0)
   (vo2, eo2, bo2, bop2, b4, b5),cov_x, infodict,mesg, success = leastsq(errfunc_2, p2[:],
                                                            args=(cellv, energy), full_output=1)

   assert success > 0
   
   print
   print 'Higher order fit:'
   print 'Volume vo2 =', vo2, 'A^3'
   print 'Energy eo2 =', eo2, 'eV'
   print 'Bulk modulus bo2 = ', bo2, 'eV/A^3 =', bo2*EV_A3_IN_GPA, 'GPa'
   print 'dP/dV (P=0) bop2 = ', bop2
   print 'b4 =', b4
   print 'b5 =', b5
   

volgrid = linspace(cellv.min(), cellv.max(), 100)
figure(1); clf()
plot(cellv, energy, 'o', label='Calculation')
plot(volgrid, birch_energy(vo,eo,bo,bop,v=volgrid), label='Fit to third-order Birch polynomial')
xlabel(r'Cell Volume / $\AA^3$')
ylabel(r'Energy / eV')
title(r'Bulk Modulus')

if do_higher_order:
   plot(volgrid, birch_energy_2(vo2,eo2,bo2,bop2,b4,b5,v=volgrid), label='Birch Fit 2')

legend(loc='upper left')

figure(2); clf()
plot(pressure, cellv, 'o', label='Calculation')
xlabel('Pressure / GPa')
ylabel(r'Cell Volume / $\AA^3$')
title(r'Equation of State')
plot(birch_pressure(vo,bo,bop,volgrid)*EV_A3_IN_GPA, volgrid, label='Fit to third-order Birch polynomial')
legend()

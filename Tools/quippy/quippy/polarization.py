"""This module contains utility functions for polarizable potentials."""

from quippy.units import BOHR, HARTREE, PI
from quippy.farray import fzeros, frange, fidentity
import numpy as np

def born_effective_charge(pot, at0, dx=1e-5, args_str=None):
   """
   Calculate Born effective charges for all atoms in at0

   Potential must be polarizable, i.e. compute dipole moments.
   """

   born_tensor = fzeros((3,3,at0.n))

   restart = True
   for i in frange(at0.n):
      for j in frange(3):

         at = at0.copy()
         at.pos[j,i] -= dx
         at.calc_connect()

         pot.calc(at, force=True, restart=restart, args_str=args_str)
         restart = True

         dip1 = fzeros(3)
         for k in frange(at.n):
            dip1 += at.dipoles[k] + at.charge[k]*at.pos[:,k]

         at = at0.copy()
         at.pos[j,i] += dx
         at.calc_connect()

         pot.calc(at, force=True, restart=restart, args_str=args_str)
         
         dip2 = fzeros(3)
         for k in frange(at.n):
            dip2 += at.dipoles[k] + at.charge[k]*at.pos[:,k]

         born_tensor[:,j,i] = (dip2 - dip1)/(dx*2.0)

   return born_tensor
         

def epsilon_infty(pot, at, deltafield=0.001, zerotol=1e-5):
   """
   Calculate dielectric tensor.

   Potential must be polarisable and allow an external electic field to be applied.
   """

   dielectric_tensor = fzeros((3,3))
   celldip = fzeros((3,3))

   # Calculation with no applied field
   pot.calc(at, force=True, restart=True, applied_efield=False)
   celldip0 = at.dipoles.sum(axis=2)/BOHR

   at.add_property('ext_efield', 0., n_cols=3)

   # Now we apply field along each of x,y,z in turn, and
   # calculate overall cell dipole moment
   for i in frange(3):

      at.ext_efield[:] = 0.0
      at.ext_efield[i,:] += deltafield/(BOHR/HARTREE)

      pot.calc(at, force=True, applied_efield=True)
      celldip[:,i] = at.dipoles.sum(axis=2)/BOHR
      dielectric_tensor[:,i] = celldip[:,i] - celldip0
      
   dielectric_tensor = 4.0*PI*dielectric_tensor/deltafield/(at.cell_volume()/BOHR**3) + fidentity(3)
   at.ext_efield[:] = 0.0

   dielectric_tensor[dielectric_tensor < zerotol] = 0.0

   return dielectric_tensor
      

def calc_effective_charge_vectors(a, born_tensor):

   effective_charge = fzeros(3)
   for i in frange(a.n):
      p_norm = a.phonon[i]/np.sqrt(np.dot(a.phonon[i],a.phonon[i]))
      disp = (p_norm/np.sqrt(ElementMass[str(a.species[i])]*MASSCONVERT))/BOHR
      print disp
      effective_charge = effective_charge + dot(born_tensor[:,:,i], disp)

   return effective_charge

def screened_effective_charge(born, eps):
   """
   Compute screened effective charge tensor from Born and dielectric tensors
   """

   screened = fzeros((3,3,3,3))
   for i in frange(3):
      for j in frange(3):
         for k in frange(3):
            for l in frange(3):
               if eps[k,l] == 0.0: continue
               screened[i,j,k,l] = born[i,j]/np.sqrt(eps[k,l])

   return screened

def dipole_moment(at):
    return np.dot(at.pos.T, at.charge)


def silica_dipole_moment(at, mask=None, q_si=1.0):
    if mask is None:
        mask = [True for i in at.indices]
    return q_si*(at.pos[:,(at.z == 14) & mask].sum(axis=2) -
                 0.5*at.pos[:,(at.z == 8) & mask].sum(axis=2) +
                 0.25*at.pos[:,(at.z == 1) & mask].sum(axis=2))

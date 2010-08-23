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

"""Simple interface between :mod:`quippy` and :mod:`ase` <https://wiki.fysik.dtu.dk/ase>"""

def atoms_to_ase(at):
   """Return a copy of `at` as an :class:`ase.Atoms` instance.
   `at` should be a :class:`quippy.Atoms` instance."""
   
   from numpy import array
   from quippy import Atoms as quippy_Atoms, ElementMass
   from ase import Atoms as ase_Atoms

   if not isinstance(at, quippy_Atoms):
      raise ValueError('at (%r) is not an instance of "quippy.Atoms" class' % at)

   mass = None
   if at.has_property('mass'): mass = at.mass

   charge = None
   if at.has_property('charge'): charge = at.charge

   ase_at = ase_Atoms(symbols=list(at.species.stripstrings()),
                      positions=array(at.pos.T),
                      masses=mass,
                      magmoms=None,
                      charges=charge,
                      cell=array(at.lattice.T),
                      pbc=True)

   if at.has_property('velo'):
      ase_at.set_velocities(at.velo.T)

   return ase_at

                
def atoms_from_ase(at):
   """Return a copy of `at` as a :class:`quippy.Atoms` instance.
   `at` should be an :class:`ase.Atoms` instance, or something compatible.
   """
   
   from quippy import Atoms as quippy_Atoms, ElementMass
   from ase import Atoms as ase_Atoms

   qat = quippy_Atoms(n=at.get_number_of_atoms(),
                      lattice=at.get_cell().T)

   qat.pos[:] = at.get_positions().T
   qat.set_atoms(at.get_atomic_numbers())

   if at.has('masses'):
      qat.add_property('mass', at.get_masses())

   if at.has('momenta'):
      qat.add_property('velo', at.get_velocities().T)

   if at.has('charges'):
      qat.add_property('charge', at.get_charges())

   return qat

from quippy import Atoms
class ASECompatibleAtoms(Atoms):
   """Subclass of quippy.Atoms which is compatible with the ASE Atoms class."""

   def __init__(self, *args, **kwargs):
      try:
         from ase import Atoms as ase_Atoms
         Atoms.__init__(atoms_from_ase(ase_Atoms(*args, **kwargs)))
      except ImportError:
         raise ImportError('ASE not installed.')
   

class CalculatorMixin(object):
   def update(self, atoms):
      qat = atoms_from_ase(atoms)
      if not hasattr(self, 'quippy_at') or self.quippy_at != qat:
         from numpy import array, zeros
         from quippy import GPA

         self.quippy_at = qat.copy()
         self.quippy_at.set_cutoff(self.cutoff())
         self.quippy_at.calc_connect()

         energy = array(0.0)
         force = zeros((3,self.quippy_at.n), order='F')
         virial = zeros((3,3), order='F')

         self.calc(self.quippy_at, force=force, energy=energy, virial=virial)

         self.energy = float(energy)
         self.forces = array(force.T)
         self.stress = -array(virial)*GPA/qat.cell_volume()
      
   def get_potential_energy(self, atoms):
      self.update(atoms)
      return self.energy

   def get_forces(self, atoms):
      self.update(atoms)
      return self.forces.copy()

   def get_numeric_forces(self, atoms):
      from numpy import zeros
      
      self.update(atoms)
      df = zeros((3,self.quippy_at.n), order='F')
      self.calc(self.quippy_at, force=df, force_fd=True)
      return df.T

   def get_stress(self, atoms):
      self.update(atoms)
      return self.stress

from quippy import Potential
class PotentialCalculator(Potential, CalculatorMixin):
   pass


# If the ASE is installed, register AtomsReader from ase.Atoms instances
# and add a 'to_ase()' method to quippy Atoms class
try:
   from ase import Atoms as ase_Atoms
   from quippy import Atoms, Potential, AtomsReaders

   def atoms_from_ase_gen(ase_at):
      yield atoms_from_ase(ase_at)

   AtomsReaders[ase_Atoms] = atoms_from_ase_gen
   AtomsReaders['ase'] = atoms_from_ase_gen
   Atoms.to_ase = atoms_to_ase
   del atoms_from_ase_gen
   
except ImportError:
   pass


from ase import Atoms as aseAtoms
from quippy import (Atoms as quippyAtoms,
                    PROPERTY_INT, PROPERTY_REAL, PROPERTY_STR, PROPERTY_LOGICAL)
import numpy as np

class MixedAtoms(aseAtoms, quippyAtoms):

   def __new__(cls, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None,
                 constraint=None,
                 calculator=None):
      return aseAtoms.__new__(cls, symbols, positions, numbers, tags, momenta, masses,
                        magmoms, charges, scaled_positions, cell, pbc, constraint, calculator)

   def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None,
                 constraint=None,
                 calculator=None):

      if cell is None:
         lattice = np.eye(3)
      else:
         lattice = cell.T
      quippyAtoms.__init__(self, n=len(symbols), lattice=lattice)
      aseAtoms.__init__(self, symbols, positions, numbers, tags, momenta, masses,
                        magmoms, charges, scaled_positions, cell, pbc,
                        constraint, calculator)
      self.set_atoms(self.z) # initialise species from z

   def new_array(self, name, a, dtype=None, shape=None):
      """Add new array.

      If *shape* is not *None*, the shape of *a* will be checked.
      Overridden to store the data in the Table `self.data`.
      """

      type_map = {'i': PROPERTY_INT,
                  'f': PROPERTY_REAL,
                  'S': PROPERTY_STR,
                  'b': PROPERTY_LOGICAL}

      name_map = {'positions': 'pos',
                  'masses': 'mass',
                  'numbers': 'z'}
            
      a = np.array(a, dtype=dtype)

      if name not in name_map and name in self.arrays:
         raise RuntimeError

      for b in self.arrays.values():
         if len(a) != len(b):
            raise ValueError('Array has wrong length: %d != %d.' %
                             (len(a), len(b)))
         break

      if shape is not None and a.shape[1:] != shape:
         raise ValueError('Array has wrong shape %s != %s.' %
                          (a.shape, (a.shape[0:1] + shape)))

      quippy_name = name_map.get(name, name)

      self.add_property(quippy_name, a.T, property_type=type_map[a.dtype.kind])

   def get_cell(self):
      return self.lattice.view(np.ndarray).T.copy()

   def set_cell(self, cell, scale_positions=False, fix=None):
      quippyAtoms.set_lattice(self, cell.T, scale_positions)
      aseAtoms.set_cell(self, cell, scale_positions, fix)

   def _getcell(self):
      return self.lattice.view(np.ndarray).T

   cell = property(_getcell, set_cell)

   def write(self, filename, format=None, properties=None, *args, **kwargs):
      try:
         quippyAtoms.write(self, filename, format=format, properties=properties, *args, **kwargs)
      except ValueError:
         aseAtoms.write(self, filename, format=format)


   def _update_hook(self):
      quippyAtoms._update_hook(self)
      self.arrays = {}     
      for name in self._props:
         name_map = {'pos': 'positions',
                     'mass': 'masses',
                     'z': 'numbers'}
         ase_name = name_map.get(name, name)
         self.arrays[ase_name] = getattr(self, name).view(np.ndarray).T

   



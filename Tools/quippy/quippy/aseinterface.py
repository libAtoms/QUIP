"""Simple interface between :mod:`quippy` and :mod:`ase` <ttps://wiki.fysik.dtu.dk/ase>"""

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

   
class QuippyCalculator(object):

   def __init__(self, pot):
      """Construct a QuippyCalculator object from `pot`, which
      should either be an instance of :class:`quippy.Potential` or
      :class:`quippy.MetaPotential`"""
      
      self.pot = pot
      self.quippy_at = None

   def update(self, atoms):
      qat = atoms_from_ase(atoms)
      if self.quippy_at is None or self.quippy_at != qat:
         self.calculate(qat)

   def calculate(self, qat):
      from numpy import array, zeros
      from quippy import GPA
      
      self.quippy_at = qat.copy()
      self.quippy_at.set_cutoff(self.pot.cutoff())
      self.quippy_at.calc_connect()
      
      energy = array(0.0)
      force = zeros((3,self.quippy_at.n), order='F')
      virial = zeros((3,3), order='F')

      self.pot.calc(self.quippy_at, f=force, e=energy, virial=virial)

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
      self.pot.calc(self.quippy_at, df=df)
      return df.T

   def get_stress(self, atoms):
      self.update(atoms)
      return self.stress

import castep, os

class CastepCalculator(object):
   def __init__(self, cell=None, param=None, castep_exec=None, stem=None):
      """Construct a CastepCalculator object from `cell` and
      `param` templates."""
      
      if isinstance(cell, str):
         self.cell = castep.CastepCell(cell)
      else:
         self.cell = cell

      if isinstance(param, str):
         self.param = castep.CastepParam(param)
      else:
         self.param = param

      if stem is None:
         self.stem = 'castep_calculator'
      else:
         self.stem = stem

      self.castep_exec = castep_exec
      self.quippy_at = None
      self.n = 0


   def update(self, atoms):
      qat = atoms_from_ase(atoms)
      if self.quippy_at is None or self.quippy_at != qat:
         self.calculate(qat)

   def calculate(self, qat):
      from numpy import array, zeros
      from quippy import GPA
      
      self.quippy_at = qat.copy()
      self.cell.update_from_atoms(self.quippy_at)

      stem = '%s_%05d' % (self.stem, self.n)
      castep.run_castep(self.cell, self.param, stem, self.castep_exec, test_mode=True)
      self.n += 1

      try:
         print 'Reading from file %s' % (stem+'.castep')
         result = Atoms(stem+'.castep', atoms_ref=self.quippy_at)
         print self.quippy_at.pos
         print result.pos
         self.energy = float(result.energy)
         self.forces = array(result.force.T)
         if hasattr(result, 'virial'):
            self.stress = -array(result.virial)*GPA/result.cell_volume()
         else:
            self.stress = zeros((3,3))
      except IOError:
         self.energy = 0.0
         self.forces = zeros((qat.n,3))
         self.stress = zeros((3,3))

      
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
      self.pot.calc(self.quippy_at, df=df)
      return df.T

   def get_stress(self, atoms):
      self.update(atoms)
      return self.stress
   


# If the ASE is installed, register AtomsReader from ase.Atoms instances
# and add a 'to_ase()' method to quippy Atoms and Potential classes
try:
   from ase import Atoms as ase_Atoms
   from quippy import Atoms, Potential, AtomsReaders

   def atoms_from_ase_gen(ase_at):
      yield atoms_from_ase(ase_at)

   AtomsReaders[ase_Atoms] = atoms_from_ase_gen
   AtomsReaders['ase'] = atoms_from_ase_gen
   Atoms.to_ase = atoms_to_ase
   del atoms_from_ase_gen

   def potential_to_ase(self):
      return QuippyCalculator(self)

   Potential.to_ase = potential_to_ase
   del potential_to_ase
   
except ImportError:
   pass

   
      

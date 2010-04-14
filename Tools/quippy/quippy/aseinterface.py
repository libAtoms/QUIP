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

         self.calc(self.quippy_at, f=force, e=energy, virial=virial)

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
      self.calc(self.quippy_at, df=df)
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

   
      

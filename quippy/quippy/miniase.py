# Copyright 2008, 2009 CAMd
# (see accompanying license files for details).

# ASE is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 2.1 of the License, or
# (at your option) any later version.
#
# ASE is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with ASE.  If not, see <http://www.gnu.org/licenses/>.

## subset of ASE containing just Atom and Atoms classes and necessary data

from math import cos, sin

import numpy as np

#from ase.atom import Atom

## from ase/data/__init__.py

chemical_symbols = ['X',  'H',  'He', 'Li', 'Be',
                    'B',  'C',  'N',  'O',  'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P',  'S',  'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V',  'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I',  'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U',  'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

atomic_numbers = {}
for Z, symbol in enumerate(chemical_symbols):
    atomic_numbers[symbol] = Z

atomic_names = [
    '', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Unnilquadium', 'Unnilpentium', 'Unnilhexium']

atomic_masses = np.array([
   0.00000, # X
   1.00794, # H
   4.00260, # He
   6.94100, # Li
   9.01218, # Be
  10.81100, # B
  12.01100, # C
  14.00670, # N
  15.99940, # O
  18.99840, # F
  20.17970, # Ne
  22.98977, # Na
  24.30500, # Mg
  26.98154, # Al
  28.08550, # Si
  30.97376, # P
  32.06600, # S
  35.45270, # Cl
  39.94800, # Ar
  39.09830, # K
  40.07800, # Ca
  44.95590, # Sc
  47.88000, # Ti
  50.94150, # V
  51.99600, # Cr
  54.93800, # Mn
  55.84700, # Fe
  58.93320, # Co
  58.69340, # Ni
  63.54600, # Cu
  65.39000, # Zn
  69.72300, # Ga
  72.61000, # Ge
  74.92160, # As
  78.96000, # Se
  79.90400, # Br
  83.80000, # Kr
  85.46780, # Rb
  87.62000, # Sr
  88.90590, # Y
  91.22400, # Zr
  92.90640, # Nb
  95.94000, # Mo
    np.nan, # Tc
 101.07000, # Ru
 102.90550, # Rh
 106.42000, # Pd
 107.86800, # Ag
 112.41000, # Cd
 114.82000, # In
 118.71000, # Sn
 121.75700, # Sb
 127.60000, # Te
 126.90450, # I
 131.29000, # Xe
 132.90540, # Cs
 137.33000, # Ba
 138.90550, # La
 140.12000, # Ce
 140.90770, # Pr
 144.24000, # Nd
    np.nan, # Pm
 150.36000, # Sm
 151.96500, # Eu
 157.25000, # Gd
 158.92530, # Tb
 162.50000, # Dy
 164.93030, # Ho
 167.26000, # Er
 168.93420, # Tm
 173.04000, # Yb
 174.96700, # Lu
 178.49000, # Hf
 180.94790, # Ta
 183.85000, # W
 186.20700, # Re
 190.20000, # Os
 192.22000, # Ir
 195.08000, # Pt
 196.96650, # Au
 200.59000, # Hg
 204.38300, # Tl
 207.20000, # Pb
 208.98040, # Bi
    np.nan, # Po
    np.nan, # At
    np.nan, # Rn
    np.nan, # Fr
 226.02540, # Ra
    np.nan, # Ac
 232.03810, # Th
 231.03590, # Pa
 238.02900, # U
 237.04820, # Np
    np.nan, # Pu
    np.nan, # Am
    np.nan, # Cm
    np.nan, # Bk
    np.nan, # Cf
    np.nan, # Es
    np.nan, # Fm
    np.nan, # Md
    np.nan, # No
    np.nan])# Lw

## contents of ase/atom.py

"""This module defines the Atom object."""

import warnings
import numpy as np

#from ase.data import atomic_numbers, chemical_symbols, atomic_masses

#         singular,    plural,     default value
names = {'position': ('positions', np.zeros(3)),
         'number':   ('numbers',   0),
         'tag':      ('tags',      0),
         'momentum': ('momenta',   np.zeros(3)),
         'mass':     ('masses',    None),
         'magmom':   ('magmoms',   0.0),
         'charge':   ('charges',   0.0)
         }


def atomproperty(name, doc):
    """Helper function to easily create Atom attribute property."""

    def getter(self):
        return self.get(name)

    def setter(self, value):
        self.set(name, value)

    def deleter(self):
        self.delete(name)

    return property(getter, setter, deleter, doc)


def xyzproperty(index):
    """Helper function to easily create Atom XYZ-property."""

    def getter(self):
        return self.position[index]

    def setter(self, value):
        self.position[index] = value

    return property(getter, setter, doc='XYZ'[index] + '-coordinate')


class Atom(object):
    """Class for representing a single atom.

    Parameters:

    symbol: str or int
        Can be a chemical symbol (str) or an atomic number (int).
    position: sequence of 3 floats
        Atomi position.
    tag: int
        Special purpose tag.
    momentum: sequence of 3 floats
        Momentum for atom.
    mass: float
        Atomic mass in atomic units.
    magmom: float or 3 floats
        Magnetic moment.
    charge: float
        Atomic charge.
    """
    __slots__ = ['data', 'atoms', 'index']

    def __init__(self, symbol='X', position=(0, 0, 0),
                 tag=None, momentum=None, mass=None,
                 magmom=None, charge=None,
                 atoms=None, index=None):

        self.data = d = {}

        if atoms is None:
            # This atom is not part of any Atoms object:
            if isinstance(symbol, str):
                d['number'] = atomic_numbers[symbol]
            else:
                d['number'] = symbol
            d['position'] = np.array(position, float)
            d['tag'] = tag
            if momentum is not None:
                momentum = np.array(momentum, float)
            d['momentum'] = momentum
            d['mass'] = mass
            if magmom is not None:
                magmom = np.array(magmom, float)
            d['magmom'] = magmom
            d['charge'] = charge

        self.index = index
        self.atoms = atoms

    def __repr__(self):
        s = "Atom('%s', %s" % (self.symbol, list(self.position))
        for name in ['tag', 'momentum', 'mass', 'magmom', 'charge']:
            value = self.get_raw(name)
            if value is not None:
                if isinstance(value, np.ndarray):
                    value = value.tolist()
                s += ', %s=%s' % (name, value)
        if self.atoms is None:
            s += ')'
        else:
            s += ', index=%d)' % self.index
        return s

    def cut_reference_to_atoms(self):
        """Cut reference to atoms object."""
        for name in names:
            self.data[name] = self.get_raw(name)
        self.index = None
        self.atoms = None

    def get_raw(self, name):
        """Get attribute, return None if not explicitely set."""
        if name == 'symbol':
            return chemical_symbols[self.get_raw('number')]

        if self.atoms is None:
            return self.data[name]

        plural = names[name][0]
        if plural in self.atoms.arrays:
            return self.atoms.arrays[plural][self.index]
        else:
            return None

    def get(self, name):
        """Get attribute, return default if not explicitely set."""
        value = self.get_raw(name)
        if value is None:
            if name == 'mass':
                value = atomic_masses[self.number]
            else:
                value = names[name][1]
        return value

    def set(self, name, value):
        """Set attribute."""
        if name == 'symbol':
            name = 'number'
            value = atomic_numbers[value]

        if self.atoms is None:
            assert name in names
            self.data[name] = value
        else:
            plural, default = names[name]
            if plural in self.atoms.arrays:
                array = self.atoms.arrays[plural]
                if name == 'magmom' and array.ndim == 2:
                    assert len(value) == 3
                array[self.index] = value
            else:
                if name == 'magmom' and np.asarray(value).ndim == 1:
                    array = np.zeros((len(self.atoms), 3))
                elif name == 'mass':
                    array = self.atoms.get_masses()
                else:
                    default = np.asarray(default)
                    array = np.zeros((len(self.atoms),) + default.shape,
                                     default.dtype)
                array[self.index] = value
                self.atoms.new_array(plural, array)

    def delete(self, name):
        """Delete attribute."""
        assert self.atoms is None
        assert name not in ['number', 'symbol', 'position']
        self.data[name] = None

    symbol = atomproperty('symbol', 'Chemical symbol')
    number = atomproperty('number', 'Atomic number')
    position = atomproperty('position', 'XYZ-coordinates')
    tag = atomproperty('tag', 'Integer tag')
    momentum = atomproperty('momentum', 'XYZ-momentum')
    mass = atomproperty('mass', 'Atomic mass')
    magmom = atomproperty('magmom', 'Initial magnetic moment')
    charge = atomproperty('charge', 'Atomic charge')
    x = xyzproperty(0)
    y = xyzproperty(1)
    z = xyzproperty(2)

    def _get(self, name):
        """Helper function for deprecated get methods."""
        warnings.warn('Use atom.%s' % name, stacklevel=3)
        return getattr(self, name)

    def _set(self, name, value):
        """Helper function for deprecated set methods."""
        warnings.warn('Use atom.%s = ...' % name, stacklevel=3)
        setattr(self, name, value)

    def get_symbol(self): return self._get('symbol')
    def get_atomic_number(self): return self._get('number')
    def get_position(self): return self._get('position')
    def get_tag(self): return self._get('tag')
    def get_momentum(self): return self._get('momentum')
    def get_mass(self): return self._get('mass')
    def get_initial_magnetic_moment(self): return self._get('magmom')
    def get_charge(self): return self._get('charge')

    def set_symbol(self, value): self._set('symbol', value)
    def set_atomic_number(self, value): self._set('number', value)
    def set_position(self, value): self._set('position', value)
    def set_tag(self, value): self._set('tag', value)
    def set_momentum(self, value): self._set('momentum', value)
    def set_mass(self, value): self._set('mass', value)
    def set_initial_magnetic_moment(self, value): self._set('magmom', value)
    def set_charge(self, value): self._set('charge', value)


## pasted from ase.units

from math import pi, sqrt

# Constants from Konrad Hinsen's PhysicalQuantities module (1986 CODATA):
_c = 299792458.              # speed of light, m/s
_mu0 = 4.e-7 * pi            # permeability of vacuum
_eps0 = 1 / _mu0 / _c**2     # permittivity of vacuum
_Grav = 6.67259e-11          # gravitational constant
_hplanck = 6.6260755e-34     # Planck constant, J s
_hbar = _hplanck / (2 * pi)  # Planck constant / 2pi, J s
_e = 1.60217733e-19          # elementary charge
_me = 9.1093897e-31          # electron mass
_mp = 1.6726231e-27          # proton mass
_Nav = 6.0221367e23          # Avogadro number
_k = 1.380658e-23            # Boltzmann constant, J/K
_amu = 1.6605402e-27         # atomic mass unit, kg

Ang = Angstrom = 1.0
nm = 10.0
Bohr = 4e10 * pi * _eps0 * _hbar**2 / _me / _e**2  # Bohr radius

eV = 1.0
Hartree = _me * _e**3 / 16 / pi**2 / _eps0**2 / _hbar**2
kJ = 1000.0 / _e
kcal = 4.184 * kJ
mol = _Nav
Rydberg = 0.5 * Hartree
Ry = Rydberg
Ha = Hartree

second = 1e10 * sqrt(_e / _amu)
fs = 1e-15 * second

kB = _k / _e                 # Boltzmann constant, eV/K

Pascal = (1 / _e) / 1e30  # J/m^3
GPa = 1e9 * Pascal

Debye = 1e11 *_e * _c
alpha = _e**2 / (4 * pi * _eps0) / _hbar / _c # fine structure constant

# Derived atomic units that have no assigned name:
_aut = _hbar / (alpha**2 * _me * _c**2)      # atomic unit of time, s
_auv =  _e**2 / _hbar / (4 * pi * _eps0)     # atomic unit of velocity, m/s
_auf = alpha**3 * _me**2 * _c**3 / _hbar     # atomic unit of force, N
_aup = alpha**5 * _me**4 * _c**5 / _hbar**3  # atomic unit of pressure, Pa

AUT = second * _aut

del pi, sqrt

#from ase.data import atomic_numbers, chemical_symbols, atomic_masses
#import ase.units as units

"""Definition of the Atoms class.

This module defines the central object in the ASE package: the Atoms
object.
"""

class Atoms(object):
    """Atoms object.

    The Atoms object can represent an isolated molecule, or a
    periodically repeated structure.  It has a unit cell and
    there may be periodic boundary conditions along any of the three
    unit cell axes.

    Information about the atoms (atomic numbers and position) is
    stored in ndarrays.  Optionally, there can be information about
    tags, momenta, masses, magnetic moments and charges.

    In order to calculate energies, forces and stresses, a calculator
    object has to attached to the atoms object.

    Parameters:

    symbols: str (formula) or list of str
        Can be a string formula, a list of symbols or a list of
        Atom objects.  Examples: 'H2O', 'COPt12', ['H', 'H', 'O'],
        [Atom('Ne', (x, y, z)), ...].
    positions: list of xyz-positions
        Atomic positions.  Anything that can be converted to an
        ndarray of shape (n, 3) will do: [(x1,y1,z1), (x2,y2,z2),
        ...].
    scaled_positions: list of scaled-positions
        Like positions, but given in units of the unit cell.
        Can not be set at the same time as positions.
    numbers: list of int
        Atomic numbers (use only one of symbols/numbers).
    tags: list of int
        Special purpose tags.
    momenta: list of xyz-momenta
        Momenta for all atoms.
    masses: list of float
        Atomic masses in atomic units.
    magmoms: list of float or list of xyz-values
        Magnetic moments.  Can be either a single value for each atom
        for collinear calculations or three numbers for each atom for
        non-collinear calculations.
    charges: list of float
        Atomic charges.
    cell: 3x3 matrix
        Unit cell vectors.  Can also be given as just three
        numbers for orthorhombic cells.  Default value: [1, 1, 1].
    pbc: one or three bool
        Periodic boundary conditions flags.  Examples: True,
        False, 0, 1, (1, 1, 0), (True, False, False).  Default
        value: False.
    constraint: constraint object(s)
        Used for applying one or more constraints during structure
        optimization.
    calculator: calculator object
        Used to attach a calculator for calculating energies and atomic
        forces.
    info: dict of key-value pairs
        Dictionary of key-value pairs with additional information
        about the system.  The following keys may be used by ase:

          - spacegroup: Spacegroup instance
          - unit_cell: 'conventional' | 'primitive' | int | 3 ints
          - adsorbate_info:

        Items in the info attribute survives copy and slicing and can
        be store to and retrieved from trajectory files given that the
        key is a string, the value is picklable and, if the value is a
        user-defined object, its base class is importable.  One should
        not make any assumptions about the existence of keys.

    Examples:

    These three are equivalent:

    >>> d = 1.104  # N2 bondlength
    >>> a = Atoms('N2', [(0, 0, 0), (0, 0, d)])
    >>> a = Atoms(numbers=[7, 7], positions=[(0, 0, 0), (0, 0, d)])
    >>> a = Atoms([Atom('N', (0, 0, 0)), Atom('N', (0, 0, d)])

    FCC gold:

    >>> a = 4.05  # Gold lattice constant
    >>> b = a / 2
    >>> fcc = Atoms('Au',
    ...             cell=[(0, b, b), (b, 0, b), (b, b, 0)],
    ...             pbc=True)

    Hydrogen wire:

    >>> d = 0.9  # H-H distance
    >>> L = 7.0
    >>> h = Atoms('H', positions=[(0, L / 2, L / 2)],
    ...           cell=(d, L, L),
    ...           pbc=(1, 0, 0))
    """

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None,
                 constraint=None,
                 calculator=None,
                 info=None):

        atoms = None

        if hasattr(symbols, 'GetUnitCell'):
            from ase.old import OldASEListOfAtomsWrapper
            atoms = OldASEListOfAtomsWrapper(symbols)
            symbols = None
        elif hasattr(symbols, 'get_positions'):
            atoms = symbols
            symbols = None
        elif (isinstance(symbols, (list, tuple)) and
              len(symbols) > 0 and isinstance(symbols[0], Atom)):
            # Get data from a list or tuple of Atom objects:
            data = [[atom.get_raw(name) for atom in symbols]
                    for name in
                    ['position', 'number', 'tag', 'momentum',
                     'mass', 'magmom', 'charge']]
            atoms = self.__class__(None, *data)
            symbols = None

        if atoms is not None:
            # Get data from another Atoms object:
            if scaled_positions is not None:
                raise NotImplementedError
            if symbols is None and numbers is None:
                numbers = atoms.get_atomic_numbers()
            if positions is None:
                positions = atoms.get_positions()
            if tags is None and atoms.has('tags'):
                tags = atoms.get_tags()
            if momenta is None and atoms.has('momenta'):
                momenta = atoms.get_momenta()
            if magmoms is None and atoms.has('magmoms'):
                magmoms = atoms.get_initial_magnetic_moments()
            if masses is None and atoms.has('masses'):
                masses = atoms.get_masses()
            if charges is None and atoms.has('charges'):
                charges = atoms.get_charges()
            if cell is None:
                cell = atoms.get_cell()
            if pbc is None:
                pbc = atoms.get_pbc()
            if constraint is None:
                constraint = [c.copy() for c in atoms.constraints]
            if calculator is None:
                calculator = atoms.get_calculator()

        self.arrays = {}

        if symbols is None:
            if numbers is None:
                if positions is not None:
                    natoms = len(positions)
                elif scaled_positions is not None:
                    natoms = len(scaled_positions)
                else:
                    natoms = 0
                numbers = np.zeros(natoms, int)
            self.new_array('numbers', numbers, int)
        else:
            if numbers is not None:
                raise ValueError(
                    'Use only one of "symbols" and "numbers".')
            else:
                self.new_array('numbers', symbols2numbers(symbols), int)

        if cell is None:
            cell = np.eye(3)
        self.set_cell(cell)

        if positions is None:
            if scaled_positions is None:
                positions = np.zeros((len(self.arrays['numbers']), 3))
            else:
                positions = np.dot(scaled_positions, self._cell)
        else:
            if scaled_positions is not None:
                raise RuntimeError('Both scaled and cartesian positions set!')
        self.new_array('positions', positions, float, (3,))

        self.set_constraint(constraint)
        self.set_tags(default(tags, 0))
        self.set_momenta(default(momenta, (0.0, 0.0, 0.0)))
        self.set_masses(default(masses, None))
        self.set_initial_magnetic_moments(default(magmoms, 0.0))
        self.set_charges(default(charges, 0.0))
        if pbc is None:
            pbc = False
        self.set_pbc(pbc)

        if info is None:
            self.info = {}
        else:
            self.info = dict(info)

        self.adsorbate_info = {}

        self.set_calculator(calculator)

    def set_calculator(self, calc=None):
        """Attach calculator object."""
        if hasattr(calc, '_SetListOfAtoms'):
            from ase.old import OldASECalculatorWrapper
            calc = OldASECalculatorWrapper(calc, self)
        if hasattr(calc, 'set_atoms'):
            calc.set_atoms(self)
        self._calc = calc

    def get_calculator(self):
        """Get currently attached calculator object."""
        return self._calc

    def _del_calculator(self):
        self._calc = None

    calc = property(get_calculator, set_calculator, _del_calculator,
                    doc='Calculator object.')

    def set_constraint(self, constraint=None):
        """Apply one or more constrains.

        The *constraint* argument must be one constraint object or a
        list of constraint objects."""
        if constraint is None:
            self._constraints = []
        else:
            if isinstance(constraint, (list, tuple)):
                self._constraints = constraint
            else:
                self._constraints = [constraint]

    def _get_constraints(self):
        return self._constraints

    def _del_constraints(self):
        self._constraints = []

    constraints = property(_get_constraints, set_constraint, _del_constraints,
                           'Constraints of the atoms.')

    def set_cell(self, cell, scale_atoms=False, fix=None):
        """Set unit cell vectors.

        Parameters:

        cell :
            Unit cell.  A 3x3 matrix (the three unit cell vectors) or
            just three numbers for an orthorhombic cell.
        scale_atoms : bool
            Fix atomic positions or move atoms with the unit cell?
            Default behavior is to *not* move the atoms (scale_atoms=False).

        Examples:

        Two equivalent ways to define an orthorhombic cell:

        >>> a.set_cell([a, b, c])
        >>> a.set_cell([(a, 0, 0), (0, b, 0), (0, 0, c)])

        FCC unit cell:

        >>> a.set_cell([(0, b, b), (b, 0, b), (b, b, 0)])
        """

        if fix is not None:
            raise TypeError('Please use scale_atoms=%s' % (not fix))

        cell = np.array(cell, float)
        if cell.shape == (3,):
            cell = np.diag(cell)
        elif cell.shape != (3, 3):
            raise ValueError('Cell must be length 3 sequence or '
                             '3x3 matrix!')
        if scale_atoms:
            M = np.linalg.solve(self._cell, cell)
            self.arrays['positions'][:] = np.dot(self.arrays['positions'], M)
        self._cell = cell

    def get_cell(self):
        """Get the three unit cell vectors as a 3x3 ndarray."""
        return self._cell.copy()

    def get_reciprocal_cell(self):
        """Get the three reciprocal lattice vectors as a 3x3 ndarray.

        Note that the commonly used factor of 2 pi for Fourier
        transforms is not included here."""

        rec_unit_cell = np.linalg.inv(self.get_cell()).transpose()
        return rec_unit_cell

    def set_pbc(self, pbc):
        """Set periodic boundary condition flags."""
        if isinstance(pbc, int):
            pbc = (pbc,) * 3
        self._pbc = np.array(pbc, bool)

    def get_pbc(self):
        """Get periodic boundary condition flags."""
        return self._pbc.copy()

    def new_array(self, name, a, dtype=None, shape=None):
        """Add new array.

        If *shape* is not *None*, the shape of *a* will be checked."""

        if dtype is not None:
            a = np.array(a, dtype)
        else:
            a = a.copy()

        if name in self.arrays:
            raise RuntimeError

        for b in self.arrays.values():
            if len(a) != len(b):
                raise ValueError('Array has wrong length: %d != %d.' %
                                 (len(a), len(b)))
            break

        if shape is not None and a.shape[1:] != shape:
            raise ValueError('Array has wrong shape %s != %s.' %
                             (a.shape, (a.shape[0:1] + shape)))

        self.arrays[name] = a

    def get_array(self, name, copy=True):
        """Get an array.

        Returns a copy unless the optional argument copy is false.
        """
        if copy:
            return self.arrays[name].copy()
        else:
            return self.arrays[name]

    def set_array(self, name, a, dtype=None, shape=None):
        """Update array.

        If *shape* is not *None*, the shape of *a* will be checked.
        If *a* is *None*, then the array is deleted."""

        b = self.arrays.get(name)
        if b is None:
            if a is not None:
                self.new_array(name, a, dtype, shape)
        else:
            if a is None:
                del self.arrays[name]
            else:
                a = np.asarray(a)
                if a.shape != b.shape:
                    raise ValueError('Array has wrong shape %s != %s.' %
                                     (a.shape, b.shape))
                b[:] = a

    def has(self, name):
        """Check for existence of array.

        name must be one of: 'tags', 'momenta', 'masses', 'magmoms',
        'charges'."""
        return name in self.arrays

    def set_atomic_numbers(self, numbers):
        """Set atomic numbers."""
        self.set_array('numbers', numbers, int, ())

    def get_atomic_numbers(self):
        """Get integer array of atomic numbers."""
        return self.arrays['numbers'].copy()

    def set_chemical_symbols(self, symbols):
        """Set chemical symbols."""
        self.set_array('numbers', symbols2numbers(symbols), int, ())

    def get_chemical_symbols(self, reduce=False):
        """Get list of chemical symbol strings.

        If reduce is True, a single string is returned, where repeated
        elements have been contracted to a single symbol and a number.
        E.g. instead of ['C', 'O', 'O', 'H'], the string 'CO2H' is returned.
        """
        if not reduce:
            # XXX
            return [chemical_symbols[Z] for Z in self.arrays['numbers']]
        else:
            num = self.get_atomic_numbers()
            N = len(num)
            dis = np.concatenate(([0], np.arange(1, N)[num[1:] != num[:-1]]))
            repeat = np.append(dis[1:], N) - dis
            symbols = ''.join([chemical_symbols[num[d]] + str(r) * (r != 1)
                               for r, d in zip(repeat, dis)])
            return symbols

    def set_tags(self, tags):
        """Set tags for all atoms."""
        self.set_array('tags', tags, int, ())

    def get_tags(self):
        """Get integer array of tags."""
        if 'tags' in self.arrays:
            return self.arrays['tags'].copy()
        else:
            return np.zeros(len(self), int)

    def set_momenta(self, momenta):
        """Set momenta."""
        if len(self.constraints) > 0 and momenta is not None:
            momenta = np.array(momenta)  # modify a copy
            for constraint in self.constraints:
                constraint.adjust_forces(self.arrays['positions'], momenta)
        self.set_array('momenta', momenta, float, (3,))

    def set_velocities(self, velocities):
        """Set the momenta by specifying the velocities."""
        self.set_momenta(self.get_masses()[:, np.newaxis] * velocities)

    def get_momenta(self):
        """Get array of momenta."""
        if 'momenta' in self.arrays:
            return self.arrays['momenta'].copy()
        else:
            return np.zeros((len(self), 3))

    def set_masses(self, masses='defaults'):
        """Set atomic masses.

        The array masses should contain a list of masses.  In case
        the masses argument is not given or for those elements of the
        masses list that are None, standard values are set."""

        if masses == 'defaults':
            masses = atomic_masses[self.arrays['numbers']]
        elif isinstance(masses, (list, tuple)):
            newmasses = []
            for m, Z in zip(masses, self.arrays['numbers']):
                if m is None:
                    newmasses.append(atomic_masses[Z])
                else:
                    newmasses.append(m)
            masses = newmasses
        self.set_array('masses', masses, float, ())

    def get_masses(self):
        """Get array of masses."""
        if 'masses' in self.arrays:
            return self.arrays['masses'].copy()
        else:
            return atomic_masses[self.arrays['numbers']]

    def set_initial_magnetic_moments(self, magmoms=None):
        """Set the initial magnetic moments.

        Use either one or three numbers for every atom (collinear
        or non-collinear spins)."""

        if magmoms is None:
            self.set_array('magmoms', None)
        else:
            magmoms = np.asarray(magmoms)
            self.set_array('magmoms', magmoms, float, magmoms.shape[1:])

    def get_initial_magnetic_moments(self):
        """Get array of initial magnetic moments."""
        if 'magmoms' in self.arrays:
            return self.arrays['magmoms'].copy()
        else:
            return np.zeros(len(self))

    def get_magnetic_moments(self):
        """Get calculated local magnetic moments."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        if self._calc.get_spin_polarized():
            return self._calc.get_magnetic_moments(self)
        else:
            return np.zeros(len(self))

    def get_magnetic_moment(self):
        """Get calculated total magnetic moment."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        if self._calc.get_spin_polarized():
            return self._calc.get_magnetic_moment(self)
        else:
            return 0.0

    def set_charges(self, charges):
        """Set charges."""
        self.set_array('charges', charges, float, ())

    def get_charges(self):
        """Get array of charges."""
        if 'charges' in self.arrays:
            return self.arrays['charges'].copy()
        else:
            return np.zeros(len(self))

    def set_positions(self, newpositions):
        """Set positions."""
        positions = self.arrays['positions']
        if self.constraints:
            newpositions = np.asarray(newpositions, float)
            for constraint in self.constraints:
                constraint.adjust_positions(positions, newpositions)

        self.set_array('positions', newpositions, shape=(3,))

    def get_positions(self):
        """Get array of positions."""
        return self.arrays['positions'].copy()

    def get_calculation_done(self):
        """Let the calculator calculate its thing,
           using the current input.
           """
        if self.calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        self.calc.initialize(self)
        self.calc.calculate(self)

    def get_potential_energy(self):
        """Calculate potential energy."""
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_potential_energy(self)

    def get_potential_energies(self):
        """Calculate the potential energies of all the atoms.

        Only available with calculators supporting per-atom energies
        (e.g. classical potentials).
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_potential_energies(self)

    def get_kinetic_energy(self):
        """Get the kinetic energy."""
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return 0.0
        return 0.5 * np.vdot(momenta, self.get_velocities())

    def get_velocities(self):
        """Get array of velocities."""
        momenta = self.arrays.get('momenta')
        if momenta is None:
            return None
        m = self.arrays.get('masses')
        if m is None:
            m = atomic_masses[self.arrays['numbers']]
        return momenta / m.reshape(-1, 1)

    def get_total_energy(self):
        """Get the total energy - potential plus kinetic energy."""
        return self.get_potential_energy() + self.get_kinetic_energy()

    def get_forces(self, apply_constraint=True):
        """Calculate atomic forces.

        Ask the attached calculator to calculate the forces and apply
        constraints.  Use *apply_constraint=False* to get the raw
        forces."""

        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        forces = self._calc.get_forces(self)
        if apply_constraint:
            for constraint in self.constraints:
                constraint.adjust_forces(self.arrays['positions'], forces)
        return forces

    def get_stress(self):
        """Calculate stress tensor.

        Returns an array of the six independent components of the
        symmetric stress tensor, in the traditional order
        (s_xx, s_yy, s_zz, s_yz, s_xz, s_xy).
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        stress = self._calc.get_stress(self)
        shape = getattr(stress, 'shape', None)
        if shape == (3, 3):
            return np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                             stress[1, 2], stress[0, 2], stress[0, 1]])
        else:
            # Hopefully a 6-vector, but don't check in case some weird
            # calculator does something else.
            return stress

    def get_stresses(self):
        """Calculate the stress-tensor of all the atoms.

        Only available with calculators supporting per-atom energies and
        stresses (e.g. classical potentials).  Even for such calculators
        there is a certain arbitrariness in defining per-atom stresses.
        """
        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        return self._calc.get_stresses(self)

    def get_dipole_moment(self):
        """Calculate the electric dipole moment for the atoms object.

        Only available for calculators which has a get_dipole_moment()
        method."""

        if self._calc is None:
            raise RuntimeError('Atoms object has no calculator.')
        try:
            dipole = self._calc.get_dipole_moment(self)
        except AttributeError:
            raise AttributeError(
                'Calculator object has no get_dipole_moment method.')
        return dipole

    def copy(self):
        """Return a copy."""
        import copy
        atoms = self.__class__(cell=self._cell, pbc=self._pbc, info=self.info)

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a.copy()
        atoms.constraints = copy.deepcopy(self.constraints)
        atoms.adsorbate_info = copy.deepcopy(self.adsorbate_info)
        return atoms

    def __len__(self):
        return len(self.arrays['positions'])

    def get_number_of_atoms(self):
        """Returns the number of atoms.

        Equivalent to len(atoms) in the standard ASE Atoms class.
        """
        return len(self)

    def __repr__(self):
        num = self.get_atomic_numbers()
        N = len(num)
        if N == 0:
            symbols = ''
        elif N <= 60:
            symbols = self.get_chemical_symbols(reduce=True)
        else:
            symbols = ''.join([chemical_symbols[Z] for Z in num[:15]]) + '...'
        s = "%s(symbols='%s', " % (self.__class__.__name__, symbols)
        for name in self.arrays:
            if name == 'numbers':
                continue
            s += '%s=..., ' % name
        if (self._cell - np.diag(self._cell.diagonal())).any():
            s += 'cell=%s, ' % self._cell.tolist()
        else:
            s += 'cell=%s, ' % self._cell.diagonal().tolist()
        s += 'pbc=%s, ' % self._pbc.tolist()
        if len(self.constraints) == 1:
            s += 'constraint=%s, ' % repr(self.constraints[0])
        if len(self.constraints) > 1:
            s += 'constraint=%s, ' % repr(self.constraints)
        if self._calc is not None:
            s += 'calculator=%s(...), ' % self._calc.__class__.__name__
        return s[:-2] + ')'

    def __add__(self, other):
        atoms = self.copy()
        atoms += other
        return atoms

    def extend(self, other):
        """Extend atoms object by appending atoms from *other*."""
        if isinstance(other, Atom):
            other = self.__class__([other])

        n1 = len(self)
        n2 = len(other)

        for name, a1 in self.arrays.items():
            a = np.zeros((n1 + n2,) + a1.shape[1:], a1.dtype)
            a[:n1] = a1
            a2 = other.arrays.get(name)
            if a2 is not None:
                a[n1:] = a2
            self.arrays[name] = a

        for name, a2 in other.arrays.items():
            if name in self.arrays:
                continue
            a = np.empty((n1 + n2,) + a2.shape[1:], a2.dtype)
            a[n1:] = a2
            if name == 'masses':
                a[:n1] = self.get_masses()
            else:
                a[:n1] = 0

            self.set_array(name, a)

        return self

    __iadd__ = extend

    def append(self, atom):
        """Append atom to end."""
        self.extend(self.__class__([atom]))

    def __getitem__(self, i):
        """Return a subset of the atoms.

        i -- scalar integer, list of integers, or slice object
        describing which atoms to return.

        If i is a scalar, return an Atom object. If i is a list or a
        slice, return an Atoms object with the same cell, pbc, and
        other associated info as the original Atoms object. The
        indices of the constraints will be shuffled so that they match
        the indexing in the subset returned.

        """
        if isinstance(i, int):
            natoms = len(self)
            if i < -natoms or i >= natoms:
                raise IndexError('Index out of range.')

            return Atom(atoms=self, index=i)

        import copy
        #from ase.constraints import FixConstraint

        atoms = self.__class__(cell=self._cell, pbc=self._pbc, info=self.info)
        # TODO: Do we need to shuffle indices in adsorbate_info too?
        atoms.adsorbate_info = self.adsorbate_info

        atoms.arrays = {}
        for name, a in self.arrays.items():
            atoms.arrays[name] = a[i].copy()

        # Constraints need to be deepcopied, since we need to shuffle
        # the indices
        atoms.constraints = copy.deepcopy(self.constraints)
        condel = []
        for con in atoms.constraints:
            pass
            #if isinstance(con, FixConstraint):
            #    try:
            #        con.index_shuffle(i)
            #    except IndexError:
            #        condel.append(con)
        for con in condel:
            atoms.constraints.remove(con)
        return atoms

    def __delitem__(self, i):
        #from ase.constraints import FixAtoms
        #check_constraint = np.array([isinstance(c, FixAtoms)
        #                             for c in self._constraints])
        if len(self._constraints) > 0:  #and not check_constraint.all():
            raise RuntimeError('Remove constraint using set_constraint() ' +
                               'before deleting atoms.')
        mask = np.ones(len(self), bool)
        mask[i] = False
        for name, a in self.arrays.items():
            self.arrays[name] = a[mask]
        if len(self._constraints) > 0:
            for n in range(len(self._constraints)):
                self._constraints[n].delete_atom(range(len(mask))[i])

    def pop(self, i=-1):
        """Remove and return atom at index *i* (default last)."""
        atom = self[i]
        atom.cut_reference_to_atoms()
        del self[i]
        return atom

    def __imul__(self, m):
        """In-place repeat of atoms."""
        if isinstance(m, int):
            m = (m, m, m)

        M = np.product(m)
        n = len(self)

        for name, a in self.arrays.items():
            self.arrays[name] = np.tile(a, (M,) + (1,) * (len(a.shape) - 1))

        positions = self.arrays['positions']
        i0 = 0
        for m0 in range(m[0]):
            for m1 in range(m[1]):
                for m2 in range(m[2]):
                    i1 = i0 + n
                    positions[i0:i1] += np.dot((m0, m1, m2), self._cell)
                    i0 = i1

        if self.constraints is not None:
            self.constraints = [c.repeat(m, n) for c in self.constraints]

        self._cell = np.array([m[c] * self._cell[c] for c in range(3)])

        return self

    def repeat(self, rep):
        """Create new repeated atoms object.

        The *rep* argument should be a sequence of three positive
        integers like *(2,3,1)* or a single integer (*r*) equivalent
        to *(r,r,r)*."""

        atoms = self.copy()
        atoms *= rep
        return atoms

    __mul__ = repeat

    def translate(self, displacement):
        """Translate atomic positions.

        The displacement argument can be a float an xyz vector or an
        nx3 array (where n is the number of atoms)."""

        self.arrays['positions'] += np.array(displacement)

    def center(self, vacuum=None, axis=None):
        """Center atoms in unit cell.

        Centers the atoms in the unit cell, so there is the same
        amount of vacuum on all sides.

        Parameters:

        vacuum (default: None): If specified adjust the amount of
        vacuum when centering.  If vacuum=10.0 there will thus be 10
        Angstrom of vacuum on each side.

        axis (default: None): If specified, only act on the specified
        axis.  Default: Act on all axes.
        """
        # Find the orientations of the faces of the unit cell
        c = self.get_cell()
        dirs = np.zeros_like(c)
        for i in range(3):
            dirs[i] = np.cross(c[i - 1], c[i - 2])
            dirs[i] /= np.sqrt(np.dot(dirs[i], dirs[i]))  # normalize
            if np.dot(dirs[i], c[i]) < 0.0:
                dirs[i] *= -1

        # Now, decide how much each basis vector should be made longer
        if axis is None:
            axes = (0, 1, 2)
        else:
            axes = (axis,)
        p = self.arrays['positions']
        longer = np.zeros(3)
        shift = np.zeros(3)
        for i in axes:
            p0 = np.dot(p, dirs[i]).min()
            p1 = np.dot(p, dirs[i]).max()
            height = np.dot(c[i], dirs[i])
            if vacuum is not None:
                lng = (p1 - p0 + 2 * vacuum) - height
            else:
                lng = 0.0  # Do not change unit cell size!
            top = lng + height - p1
            shf = 0.5 * (top - p0)
            cosphi = np.dot(c[i], dirs[i]) / np.sqrt(np.dot(c[i], c[i]))
            longer[i] = lng / cosphi
            shift[i] = shf / cosphi

        # Now, do it!
        translation = np.zeros(3)
        for i in axes:
            nowlen = np.sqrt(np.dot(c[i], c[i]))
            self._cell[i] *= 1 + longer[i] / nowlen
            translation += shift[i] * c[i] / nowlen
        self.arrays['positions'] += translation

    def get_center_of_mass(self, scaled=False):
        """Get the center of mass.

        If scaled=True the center of mass in scaled coordinates
        is returned."""
        m = self.arrays.get('masses')
        if m is None:
            m = atomic_masses[self.arrays['numbers']]
        com = np.dot(m, self.arrays['positions']) / m.sum()
        if scaled:
            return np.linalg.solve(self._cell.T, com)
        else:
            return com

    def get_moments_of_inertia(self, vectors=False):
        '''Get the moments of inertia

        The three principal moments of inertia are computed from the
        eigenvalues of the inertial tensor. periodic boundary
        conditions are ignored. Units of the moments of inertia are
        amu*angstrom**2.
        '''

        com = self.get_center_of_mass()
        positions = self.get_positions()
        positions -= com  # translate center of mass to origin
        masses = self.get_masses()

        #initialize elements of the inertial tensor
        I11 = I22 = I33 = I12 = I13 = I23 = 0.0
        for i in range(len(self)):
            x, y, z = positions[i]
            m = masses[i]

            I11 += m * (y**2 + z**2)
            I22 += m * (x**2 + z**2)
            I33 += m * (x**2 + y**2)
            I12 += -m * x * y
            I13 += -m * x * z
            I23 += -m * y * z

        I = np.array([[I11, I12, I13],
                      [I12, I22, I23],
                      [I13, I23, I33]])

        evals, evecs = np.linalg.eig(I)
        if vectors:
            return evals, evecs.transpose()
        else:
            return evals

    def rotate(self, v, a=None, center=(0, 0, 0), rotate_cell=False):
        """Rotate atoms.

        Rotate the angle *a* around the vector *v*.  If *a* is not
        given, the length of *v* is used as the angle.  If *a* is a
        vector, then *v* is rotated into *a*.  The point at *center*
        is fixed.  Use *center='COM'* to fix the center of mass.
        Vectors can also be strings: 'x', '-x', 'y', ... .

        Examples:

        Rotate 90 degrees around the z-axis, so that the x-axis is
        rotated into the y-axis:

        >>> a = pi / 2
        >>> atoms.rotate('z', a)
        >>> atoms.rotate((0, 0, 1), a)
        >>> atoms.rotate('-z', -a)
        >>> atoms.rotate((0, 0, a))
        >>> atoms.rotate('x', 'y')
        """

        norm = np.linalg.norm
        v = string2vector(v)
        if a is None:
            a = norm(v)
        if isinstance(a, (float, int)):
            v /= norm(v)
            c = cos(a)
            s = sin(a)
        else:
            v2 = string2vector(a)
            v /= norm(v)
            v2 /= norm(v2)
            c = np.dot(v, v2)
            v = np.cross(v, v2)
            s = norm(v)
            # In case *v* and *a* are parallel, np.cross(v, v2) vanish
            # and can't be used as a rotation axis. However, in this
            # case any rotation axis perpendicular to v2 will do.
            eps = 1e-7
            if s < eps:
                v = np.cross((0, 0, 1), v2)
            if norm(v) < eps:
                v = np.cross((1, 0, 0), v2)
            assert norm(v) >= eps
            if s > 0:
                v /= s

        if isinstance(center, str) and center.lower() == 'com':
            center = self.get_center_of_mass()

        p = self.arrays['positions'] - center
        self.arrays['positions'][:] = (c * p -
                                       np.cross(p, s * v) +
                                       np.outer(np.dot(p, v), (1.0 - c) * v) +
                                       center)
        if rotate_cell:
            rotcell = self.get_cell()
            rotcell[:] = (c * rotcell -
                          np.cross(rotcell, s * v) +
                          np.outer(np.dot(rotcell, v), (1.0 - c) * v))
            self.set_cell(rotcell)

    def rotate_euler(self, center=(0, 0, 0), phi=0.0, theta=0.0, psi=0.0):
        """Rotate atoms via Euler angles.

        See e.g http://mathworld.wolfram.com/EulerAngles.html for explanation.

        Parameters:

        center :
            The point to rotate about. A sequence of length 3 with the
            coordinates, or 'COM' to select the center of mass.
        phi :
            The 1st rotation angle around the z axis.
        theta :
            Rotation around the x axis.
        psi :
            2nd rotation around the z axis.

        """
        if isinstance(center, str) and center.lower() == 'com':
            center = self.get_center_of_mass()
        else:
            center = np.array(center)
        # First move the molecule to the origin In contrast to MATLAB,
        # numpy broadcasts the smaller array to the larger row-wise,
        # so there is no need to play with the Kronecker product.
        rcoords = self.positions - center
        # First Euler rotation about z in matrix form
        D = np.array(((cos(phi), sin(phi), 0.),
                      (-sin(phi), cos(phi), 0.),
                      (0., 0., 1.)))
        # Second Euler rotation about x:
        C = np.array(((1., 0., 0.),
                      (0., cos(theta), sin(theta)),
                      (0., -sin(theta), cos(theta))))
        # Third Euler rotation, 2nd rotation about z:
        B = np.array(((cos(psi), sin(psi), 0.),
                      (-sin(psi), cos(psi), 0.),
                      (0., 0., 1.)))
        # Total Euler rotation
        A = np.dot(B, np.dot(C, D))
        # Do the rotation
        rcoords = np.dot(A, np.transpose(rcoords))
        # Move back to the rotation point
        self.positions = np.transpose(rcoords) + center

    def get_dihedral(self, list):
        """Calculate dihedral angle.

        Calculate dihedral angle between the vectors list[0]->list[1]
        and list[2]->list[3], where list contains the atomic indexes
        in question.
        """

        # vector 0->1, 1->2, 2->3 and their normalized cross products:
        a = self.positions[list[1]] - self.positions[list[0]]
        b = self.positions[list[2]] - self.positions[list[1]]
        c = self.positions[list[3]] - self.positions[list[2]]
        bxa = np.cross(b, a)
        bxa /= np.sqrt(np.vdot(bxa, bxa))
        cxb = np.cross(c, b)
        cxb /= np.sqrt(np.vdot(cxb, cxb))
        angle = np.vdot(bxa, cxb)
        # check for numerical trouble due to finite precision:
        if angle < -1:
            angle = -1
        if angle > 1:
            angle = 1
        angle = np.arccos(angle)
        if np.vdot(bxa, c) > 0:
            angle = 2 * np.pi - angle
        return angle

    def set_dihedral(self, list, angle, mask=None):
        """
        set the dihedral angle between vectors list[0]->list[1] and
        list[2]->list[3] by changing the atom indexed by list[3]
        if mask is not None, all the atoms described in mask
        (read: the entire subgroup) are moved

        example: the following defines a very crude
        ethane-like molecule and twists one half of it by 30 degrees.

        >>> atoms = Atoms('HHCCHH', [[-1, 1, 0], [-1, -1, 0], [0, 0, 0],
                                     [1, 0, 0], [2, 1, 0], [2, -1, 0]])
        >>> atoms.set_dihedral([1,2,3,4],7*pi/6,mask=[0,0,0,1,1,1])
        """
        # if not provided, set mask to the last atom in the
        # dihedral description
        if mask is None:
            mask = np.zeros(len(self))
            mask[list[3]] = 1
        # compute necessary in dihedral change, from current value
        current = self.get_dihedral(list)
        diff = angle - current
        # do rotation of subgroup by copying it to temporary atoms object
        # and then rotating that
        axis = self.positions[list[2]] - self.positions[list[1]]
        center = self.positions[list[2]]
        # recursive object definition might not be the most elegant thing,
        # more generally useful might be a rotation function with a mask?
        group = self.__class__()
        for i in range(len(self)):
            if mask[i]:
                group += self[i]
        group.translate(-center)
        group.rotate(axis, diff)
        group.translate(center)
        # set positions in original atoms object
        j = 0
        for i in range(len(self)):
            if mask[i]:
                self.positions[i] = group[j].get_position()
                j += 1

    def rotate_dihedral(self, list, angle, mask=None):
        """Rotate dihedral angle.

        Complementing the two routines above: rotate a group by a
        predefined dihedral angle, starting from its current
        configuration
        """
        start = self.get_dihedral(list)
        self.set_dihedral(list, angle + start, mask)

    def rattle(self, stdev=0.001, seed=42):
        """Randomly displace atoms.

        This method adds random displacements to the atomic positions,
        taking a possible constraint into account.  The random numbers are
        drawn from a normal distribution of standard deviation stdev.

        For a parallel calculation, it is important to use the same
        seed on all processors!  """

        rs = np.random.RandomState(seed)
        positions = self.arrays['positions']
        self.set_positions(positions +
                           rs.normal(scale=stdev, size=positions.shape))

    def get_distance(self, a0, a1, mic=False):
        """Return distance between two atoms.

        Use mic=True to use the Minimum Image Convention.
        """

        R = self.arrays['positions']
        D = R[a1] - R[a0]
        if mic:
            Dr = np.linalg.solve(self._cell.T, D)
            D = np.dot(Dr - np.round(Dr) * self._pbc, self._cell)
        return np.linalg.norm(D)

    def set_distance(self, a0, a1, distance, fix=0.5):
        """Set the distance between two atoms.

        Set the distance between atoms *a0* and *a1* to *distance*.
        By default, the center of the two atoms will be fixed.  Use
        *fix=0* to fix the first atom, *fix=1* to fix the second
        atom and *fix=0.5* (default) to fix the center of the bond."""

        R = self.arrays['positions']
        D = R[a1] - R[a0]
        x = 1.0 - distance / np.linalg.norm(D)
        R[a0] += (x * fix) * D
        R[a1] -= (x * (1.0 - fix)) * D

    def get_scaled_positions(self):
        """Get positions relative to unit cell.

        Atoms outside the unit cell will be wrapped into the cell in
        those directions with periodic boundary conditions so that the
        scaled coordinates are between zero and one."""

        scaled = np.linalg.solve(self._cell.T, self.arrays['positions'].T).T
        for i in range(3):
            if self._pbc[i]:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test
                scaled[:, i] %= 1.0
                scaled[:, i] %= 1.0
        return scaled

    def set_scaled_positions(self, scaled):
        """Set positions relative to unit cell."""
        self.arrays['positions'][:] = np.dot(scaled, self._cell)

    def get_temperature(self):
        """Get the temperature. in Kelvin"""
        ekin = self.get_kinetic_energy() / len(self)
        return ekin / (1.5 * kB)

    def get_isotropic_pressure(self, stress):
        """Get the current calculated pressure, assume isotropic medium.
            in Bar
        """
        if type(stress) == type(1.0) or type(stress) == type(1):
            return -stress * 1e-5 / Pascal
        elif stress.shape == (3, 3):
            return (-(stress[0, 0] + stress[1, 1] + stress[2, 2]) / 3.0) * \
                    1e-5 / Pascal
        elif stress.shape == (6,):
            return (-(stress[0] + stress[1] + stress[2]) / 3.0) * \
                   1e-5 / Pascal
        else:
            raise ValueError('The external stress has the wrong shape.')

    def __eq__(self, other):
        """Check for identity of two atoms objects.

        Identity means: same positions, atomic numbers, unit cell and
        periodic boundary conditions."""
        try:
            a = self.arrays
            b = other.arrays
            return (len(self) == len(other) and
                    (a['positions'] == b['positions']).all() and
                    (a['numbers'] == b['numbers']).all() and
                    (self._cell == other.cell).all() and
                    (self._pbc == other.pbc).all())
        except AttributeError:
            return NotImplemented

    def __ne__(self, other):
        eq = self.__eq__(other)
        if eq is NotImplemented:
            return eq
        else:
            return not eq

    __hash__ = None

    def get_volume(self):
        """Get volume of unit cell."""
        return abs(np.linalg.det(self._cell))

    def _get_positions(self):
        """Return reference to positions-array for in-place manipulations."""
        return self.arrays['positions']

    def _set_positions(self, pos):
        """Set positions directly, bypassing constraints."""
        self.arrays['positions'][:] = pos

    positions = property(_get_positions, _set_positions,
                         doc='Attribute for direct ' +
                         'manipulation of the positions.')

    def _get_atomic_numbers(self):
        """Return reference to atomic numbers for in-place
        manipulations."""
        return self.arrays['numbers']

    numbers = property(_get_atomic_numbers, set_atomic_numbers,
                       doc='Attribute for direct ' +
                       'manipulation of the atomic numbers.')

    def _get_cell(self):
        """Return reference to unit cell for in-place manipulations."""
        return self._cell

    cell = property(_get_cell, set_cell, doc='Attribute for direct ' +
                       'manipulation of the unit cell.')

    def _get_pbc(self):
        """Return reference to pbc-flags for in-place manipulations."""
        return self._pbc

    pbc = property(_get_pbc, set_pbc,
                   doc='Attribute for direct manipulation ' +
                   'of the periodic boundary condition flags.')

    def get_name(self):
        """Return a name extracted from the elements."""
        elements = {}
        for a in self:
            try:
                elements[a.symbol] += 1
            except:
                elements[a.symbol] = 1
        name = ''
        for element in elements:
            name += element
            if elements[element] > 1:
                name += str(elements[element])
        return name

    def write(self, filename, format=None):
        """Write yourself to a file."""

        raise NotImplementedError
        
        #from ase.io import write
        #write(filename, self, format)

    def edit(self):
        """Modify atoms interactively through ag viewer.

        Conflicts leading to undesirable behaviour might arise
        when matplotlib has been pre-imported with certain
        incompatible backends and while trying to use the
        plot feature inside the interactive ag. To circumvent,
        please set matplotlib.use('gtk') before calling this
        method.
        """
        raise NotImplementedError
    
        #from ase.gui.images import Images
        #from ase.gui.gui import GUI
        images = Images([self])
        gui = GUI(images)
        gui.run()
        # use atoms returned from gui:
        # (1) delete all currently available atoms
        self.set_constraint()
        for z in range(len(self)):
            self.pop()
        edited_atoms = gui.images.get_atoms(0)
        # (2) extract atoms from edit session
        self.extend(edited_atoms)
        self.set_constraint(edited_atoms._get_constraints())
        self.set_cell(edited_atoms.get_cell())
        self.set_initial_magnetic_moments(edited_atoms.get_magnetic_moments())
        self.set_tags(edited_atoms.get_tags())
        return


def string2symbols(s):
    """Convert string to list of chemical symbols."""
    n = len(s)

    if n == 0:
        return []

    c = s[0]

    if c.isdigit():
        i = 1
        while i < n and s[i].isdigit():
            i += 1
        return int(s[:i]) * string2symbols(s[i:])

    if c == '(':
        p = 0
        for i, c in enumerate(s):
            if c == '(':
                p += 1
            elif c == ')':
                p -= 1
                if p == 0:
                    break
        j = i + 1
        while j < n and s[j].isdigit():
            j += 1
        if j > i + 1:
            m = int(s[i + 1:j])
        else:
            m = 1
        return m * string2symbols(s[1:i]) + string2symbols(s[j:])

    if c.isupper():
        i = 1
        if 1 < n and s[1].islower():
            i += 1
        j = i
        while j < n and s[j].isdigit():
            j += 1
        if j > i:
            m = int(s[i:j])
        else:
            m = 1
        return m * [s[:i]] + string2symbols(s[j:])
    else:
        raise ValueError


def symbols2numbers(symbols):
    if isinstance(symbols, str):
        symbols = string2symbols(symbols)
    numbers = []
    for s in symbols:
        if isinstance(s, str):
            numbers.append(atomic_numbers[s])
        else:
            numbers.append(s)
    return numbers


def string2vector(v):
    if isinstance(v, str):
        if v[0] == '-':
            return -string2vector(v[1:])
        w = np.zeros(3)
        w['xyz'.index(v)] = 1.0
        return w
    return np.array(v, float)


def default(data, dflt):
    """Helper function for setting default values."""
    if data is None:
        return None
    elif isinstance(data, (list, tuple)):
        newdata = []
        allnone = True
        for x in data:
            if x is None:
                newdata.append(dflt)
            else:
                newdata.append(x)
                allnone = False
        if allnone:
            return None
        return newdata
    else:
        return data

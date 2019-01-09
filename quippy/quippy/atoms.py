"""
Initial try.
Atoms object and related stuff by Tamas Stenczel on py3 branch
"""

import ase
import numpy as np

import _quippy
import quippy


class Atoms(ase.Atoms):

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None, celldisp=None,
                 constraint=None,
                 calculator=None,
                 info=None):
        ase.Atoms.__init__(self, symbols=symbols, positions=positions, numbers=numbers, tags=tags, momenta=momenta,
                           masses=masses, magmoms=magmoms, charges=charges, scaled_positions=scaled_positions,
                           cell=cell, pbc=pbc, celldisp=celldisp, constraint=constraint, calculator=calculator,
                           info=info)

        self._fortran_atoms = quippy.atoms_types_module.Atoms(len(self), self.get_lattice())

        pass

    def get_lattice(self):
        """ Fortran compatibility

            fortran Atoms uses the transpose of the ase.Atoms.cell as lattice"""
        return self.get_cell().transpose()

    def set_lattice(self, lattice):
        """ Fortran compatibility

        fortran Atoms uses the transpose of the ase.Atoms.cell as lattice"""
        lattice = np.array(lattice)
        self.set_cell(lattice.transpose())

    # ase.Atoms functions which need to be rewritten here

    def __add__(self, other):
        pass

    def extend(self, other):
        """Extend atoms object by appending atoms from *other*."""
        pass

    def append(self, atom):
        """Append atom to end."""
        # self.extend(self.__class__([atom]))
        pass

    def __getitem__(self, i):
        pass

    def __delitem__(self, i):
        pass

    def pop(self, i=-1):
        """Remove and return atom at index *i* (default last)."""
        pass

    def __imul__(self, m):
        """In-place repeat of atoms."""
        pass

    def repeat(self, rep):
        """Create new repeated atoms object.

        The *rep* argument should be a sequence of three positive
        integers like *(2,3,1)* or a single integer (*r*) equivalent
        to *(r,r,r)*."""
        pass

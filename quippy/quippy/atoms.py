"""
Initial try.
Atoms object and related stuff by Tamas Stenczel on py3 branch
"""

import ase
import numpy as np

import _quippy
import quippy


class Atoms(ase.Atoms):

    """
    There is a internal fortran atoms object, which shall not be manipulated by the user.
    The main holder of all information should be the oject itself, the fortran object is only manipulated when there is
    work to do with it and then the output of a calculation is used to update the main object.
    This is especially nice because the length of the fortran object is defined, cannot change on the way. (Can it?)


    Lattice: thought about using decorators and a property, but having the same thing in two places is nay-nay
        - only functionality it would add is self.lattice being an accessible thing

    """

    def __init__(self, symbols=None,
                 positions=None, numbers=None,
                 tags=None, momenta=None, masses=None,
                 magmoms=None, charges=None,
                 scaled_positions=None,
                 cell=None, pbc=None, celldisp=None,
                 constraint=None,
                 calculator=None,
                 info=None):

        # do the ase object initialisation
        ase.Atoms.__init__(self, symbols=symbols, positions=positions, numbers=numbers, tags=tags, momenta=momenta,
                           masses=masses, magmoms=magmoms, charges=charges, scaled_positions=scaled_positions,
                           cell=cell, pbc=pbc, celldisp=celldisp, constraint=constraint, calculator=calculator,
                           info=info)
        # constuct the fortran object to have it, but don't populate it yet, only when needed
        self._fortran_atoms = quippy.atoms_types_module.Atoms(len(self), self.get_lattice())

        pass

    def get_lattice(self):
        """Fortran compatibility
        fortran Atoms uses the transpose of the ase.Atoms.cell as lattice"""
        return self.get_cell().transpose()

    def set_lattice(self, lattice):
        """Fortran compatibility
        fortran Atoms uses the transpose of the ase.Atoms.cell as lattice"""
        lattice = np.array(lattice)
        self.set_cell(lattice.transpose())

    def _update_fortran_object(self):
        """Updates the internal fortran atoms object with the current """
        self._fortran_atoms = quippy.atoms_types_module.Atoms(len(self), self.get_lattice())
        # pass on all data to it
        pass

    def _update_ase_object(self):
        """Update the main object with the data of the fortran object."""
        pass

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

"""
Conversions between ase and fortran atoms objects
"""


import ase
import numpy as np

import _quippy
import quippy


__all__ = ['ase_to_fortran']


def ase_to_fortran(ase_atoms: ase.Atoms, fortran_atoms=None):
    """
    Converter to put the info from an ase atoms object into a fortran atoms object.
    Copies everything to make sure there is not linking back.

    Checks if the

    :param ase_atoms:
    :param fortran_atoms:
    :return:
    """

    if fortran_atoms is not None:
        if isinstance(fortran_atoms, quippy.atoms_types_module.Atoms):
            # check if the length matches, otherwise make a new one in place of that
            if len(ase_atoms) != fortran_atoms.n:
                # need to regenerate the fortran atoms object
                fortran_atoms = quippy.atoms_types_module.Atoms(len(ase_atoms), ase_atoms.get_cell().T.copy())
        else:
            # TODO: decide if we want an error here or just move on and not care, as if it was None. \
            # error is better te make awareness, but makes useage harder; decide with James!!!
            raise TypeError('fortran_atoms argument is not of valid type, cannot work with it')

    else:
        # need to regenerate the fortran atoms object
        fortran_atoms = quippy.atoms_types_module.Atoms(len(ase_atoms), ase_atoms.get_cell().transpose())

    fortran_atoms.pos[:] = ase_atoms.get_positions().T.copy()
    fortran_atoms.is_periodic[:] = ase_atoms.get_pbc()   # fixme this is not making sure it is a fortran compatible arr.
    fortran_atoms.z[:] = ase_atoms.numbers

    # go through all properties
    return fortran_atoms



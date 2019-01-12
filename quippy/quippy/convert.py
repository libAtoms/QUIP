# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2019
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


"""
Conversions between ase and fortran atoms objects
"""


import ase
import numpy as np

import _quippy
import quippy


__all__ = ['ase_to_quip']


def ase_to_quip(ase_atoms: ase.Atoms, quip_atoms=None):
    """
    Converter to put the info from an ase atoms object into a quip atoms object.
    Copies everything to make sure there is not linking back.

    Checks if the

    :param ase_atoms:
    :param quip_atoms:
    :return:
    """

    if quip_atoms is not None:
        if isinstance(quip_atoms, quippy.atoms_types_module.Atoms):
            # check if the length matches, otherwise make a new one in place of that
            if len(ase_atoms) != quip_atoms.n:
                # need to regenerate the quip atoms object
                quip_atoms = quippy.atoms_types_module.Atoms(len(ase_atoms), ase_atoms.get_cell().T.copy())
        else:
            # TODO: decide if we want an error here or just move on and not care, as if it was None. \
            # error is better te make awareness, but makes useage harder; decide with James!!!
            raise TypeError('quip_atoms argument is not of valid type, cannot work with it')

    else:
        # need to regenerate the quip atoms object
        quip_atoms = quippy.atoms_types_module.Atoms(len(ase_atoms), ase_atoms.get_cell().transpose())

    quip_atoms.pos[:] = ase_atoms.get_positions().T.copy()
    quip_atoms.is_periodic[:] = ase_atoms.get_pbc()   # fixme this is not making sure it is a quip compatible arr.
    quip_atoms.z[:] = ase_atoms.numbers

    # go through all properties
    return quip_atoms



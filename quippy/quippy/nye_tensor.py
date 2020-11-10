import numpy as np

from .nye_tensor_module import calc_nye_tensor
from .convert import ase_to_quip

__all__ = ['nye_tensor']

def nye_tensor(atoms, bulk_reference, cutoff):
    """
    Compute Nye tensor for atomic structure

    Args:
        atoms (ase.atoms.Atoms): atomsation structure
        bulk_reference ([ase.atoms.Atoms): bulk unit cell
        cutoff (float): cutoff distance for neighbour list

    Returns:
        alpha: array of shape (3, 3, len(atoms))
    """
    ref_lat_quippy = ase_to_quip(bulk_reference)
    atoms_quippy = ase_to_quip(atoms)
    
    ref_lat_quippy.set_cutoff(cutoff)
    atoms_quippy.set_cutoff(cutoff)
    
    alpha = np.zeros((3, 3, len(atoms)), order='F')
    calc_nye_tensor(atoms_quippy,
                    ref_lat_quippy,
                    alpha)
    return alpha
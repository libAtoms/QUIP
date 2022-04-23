# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Portions of this code were written by
# HQ X     Tamas K. Stenczel, James Kermode
# HQ X
# HQ X   Copyright 2019
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   https://warwick.ac.uk/fac/sci/eng/staff/jrk
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


import quippy
from ase import Atoms
import numpy as np
from ase.io.extxyz import key_val_dict_to_str, key_val_str_to_dict
from quippy.convert import get_dict_arrays

__all__ = ['Descriptor']


def convert_atoms_types_iterable_method(method):
    """
    Decorator to transparently convert ASEAtoms objects into quippy Atoms, and
    to transparently iterate over a list of Atoms objects...

    Taken from py2 version
    """

    def wrapper(self, at, *args, **kw):
        if isinstance(at, quippy.atoms_types_module.Atoms):
            return method(self, at, *args, **kw)
        elif isinstance(at, Atoms):
            _quip_at = quippy.convert.ase_to_quip(at,add_arrays=True)
            return method(self, _quip_at, *args, **kw)
        else:
            return [wrapper(self, atelement, *args, **kw) for atelement in at]

    return wrapper


class Descriptor:

    def __init__(self, args_str=None, **init_kwargs):
        """
        Initialises Descriptor object and calculate number of dimensions and
        permutations.


        properties:
        - cutoff

        calculateable:
        - sizes: `n_desc, n_cross, n_index = desc.sizes(_quip_atoms)`

        """

        if args_str is None:
            args_str = key_val_dict_to_str(init_kwargs)
        else:
            args_str += ' ' + key_val_dict_to_str(init_kwargs)

        # intialise the wrapped object and hide it from the user
        self._quip_descriptor = quippy.descriptors_module.descriptor(args_str)

        # kept for compatibility with older version
        # super convoluted though :D should just rethink it at some point
        self.n_dim = self.dimensions()
        self.n_perm = self.get_n_perm()

    def __len__(self):
        return self.dimensions()

    def dimensions(self):
        return self._quip_descriptor.dimensions()

    def get_n_perm(self):
        return self._quip_descriptor.n_permutations()

    def permutations(self):
        # this is equivalent to the py2 behaviour, giving an array of the permutations
        permutation_arr = np.zeros((self.n_dim, self.n_perm), dtype='int32', order='F')
        self._quip_descriptor.permutations(permutation_arr)
        return np.copy(permutation_arr.T, order='C')

    def cutoff(self):
        # TODO: decide if adding @property is a good idea
        # could be like get_cutoff()
        return self._quip_descriptor.cutoff()

    @convert_atoms_types_iterable_method
    def sizes(self, at, args_str=None, cutoff=None):
        """
        Replicating the QUIP method, is used in the rest of the methods
        """

        # calc connectivity on the atoms object with the internal one
        self._calc_connect(at, cutoff)

        mask = None
        if args_str is not None:
            args_str_dict = key_val_str_to_dict(args_str)
            if "atom_mask_name" in args_str_dict:
                try:
                    mask = get_dict_arrays(at.properties)[args_str_dict["atom_mask_name"]].astype(bool)
                except Exception:
                    raise KeyError
        n_descriptors, n_cross = self._quip_descriptor.sizes(at,mask=mask)

        return n_descriptors, n_cross

    @convert_atoms_types_iterable_method
    def count(self, at, args_str=None):
        """
        Returns how many descriptors of this type are found in the Atoms
        object.
        """
        # fixme: is the decorator needed now?
        return self.sizes(at,args_str=args_str)[0]

    @convert_atoms_types_iterable_method
    def _calc_connect(self, at, cutoff=None):
        """
        Internal method for calculating connectivity on a quip_atoms object

        Ideally called only on quip_atoms object, but put in decorator to make sure
        :param at:
        :return:
        """

        if cutoff is not None:
            at.set_cutoff(cutoff)
        else:
            # setting to +1 is arbitrary here, the point is to set to something a bit higher than the descriptor's
            if at.cutoff < self.cutoff() + 1:
                at.set_cutoff(self.cutoff() + 1)

        # TODO: add logic to skip this if it has been calculated already for speedup
        at.calc_connect()

    @convert_atoms_types_iterable_method
    def calc_descriptor(self, at, args_str=None, cutoff=None, **calc_kwargs):
        """
        Calculates all descriptors of this type in the Atoms object, and
        returns the array of descriptor values. Does not compute gradients; use
        calc(at, grad=True, ...) for that.

        """
        try:
            return self.calc(at, False, args_str, cutoff, **calc_kwargs)['data']
        except KeyError:
            return []

    @convert_atoms_types_iterable_method
    def calc(self, at, grad=False, args_str=None, cutoff=None, **calc_kwargs):
        """
        Calculates all descriptors of this type in the Atoms object, and
        gradients if grad=True. Results can be accessed dictionary- or
        attribute-style; 'descriptor' contains descriptor values,
        'descriptor_index_0based' contains the 0-based indices of the central
        atom(s) in each descriptor, 'grad' contains gradients,
        'grad_index_0based' contains indices to gradients (descriptor, atom).
        Cutoffs and gradients of cutoffs are also returned.

        """

        # arg string and calc_args
        if args_str is None:
            args_str = key_val_dict_to_str(calc_kwargs)
        else:
            # new, for compatibility: merged if both given
            args_str += ' ' + key_val_dict_to_str(calc_kwargs)

        # calc connectivity on the atoms object with the internal one
        self._calc_connect(at, cutoff)

        # descriptor calculation
        descriptor_out_raw = self._quip_descriptor.calc(at, do_descriptor=True, do_grad_descriptor=grad,
                                                        args_str=args_str)

        # unpack to a list of dicts
        count = self.count(at, args_str=args_str)
        descriptor_out = dict()
        for i in range(count):
            # unpack to dict with the specific converter function
            mono_dict = quippy.convert.descriptor_data_mono_to_dict(descriptor_out_raw.x[i])

            # add to the result
            for key, val in mono_dict.items():
                if key in descriptor_out.keys():
                    descriptor_out[key].append(val)
                else:
                    descriptor_out[key] = [val]

        # make numpy arrays out of them
        for key, val in descriptor_out.items():
            # merge the arrays according to shape
            if key in ['data', 'ci']:
                descriptor_out[key] = np.concatenate(val, axis=0)
            elif key in ['covariance_cutoff', 'has_data']:
                descriptor_out[key] = np.array(val)
            elif key in ["pos", "grad_covariance_cutoff"]:
                # corresponds to the gradients
                descriptor_out[key] = np.concatenate([x.T for x in val])
            elif key == "grad_data":
                descriptor_out[key] = np.transpose(np.concatenate(val, axis=2), axes=(2, 1, 0))

        if "ii" in descriptor_out.keys():
            grad_index_0based = []
            for idx, ii_perdesc in enumerate(descriptor_out["ii"]):
                for ii_item in ii_perdesc:
                    grad_index_0based.append([descriptor_out["ci"][idx], ii_item])
            # same as in py2, makes iteration of gradient easier
            descriptor_out["grad_index_0based"] = np.array(grad_index_0based) - 1

        if count > 0:
            descriptor_out['data'] = descriptor_out['data'].reshape((count, -1))
        else:
            descriptor_out['data'] = np.array([[]])

        # This is a dictionary now and hence needs to be indexed as one, unlike the old version
        return descriptor_out

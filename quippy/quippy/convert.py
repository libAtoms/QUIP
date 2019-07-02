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


__all__ = ['ase_to_quip', 'descriptor_data_mono_to_dict']


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


def descriptor_data_mono_to_dict(desc_data_mono):
    """
    Returns a dictionary out of the descriptor_data_mono object with all info it contained.
    :param desc_data_mono:
    :return:
    """

    # if not isinstance(desc_data_mono, quippy.descriptors_module.descriptor_data_mono):
    #     raise TypeError('Not descriptor_data_mono given')
    #
    # # TODO: only take the ones actually needed, this is good for debuggin now though
    # out_data_dict = {'data': desc_data_mono.data,
    #                  'grad_data': desc_data_mono.grad_data,
    #                  'ci': desc_data_mono.ci,
    #                  'ii': desc_data_mono.ii,
    #                  'pos': desc_data_mono.pos,
    #                  'has_data': desc_data_mono.has_data,
    #                  'has_grad_data': desc_data_mono.has_grad_data,
    #                  'covariance_cutoff': desc_data_mono.covariance_cutoff,
    #                  'grad_covariance_cutoff': desc_data_mono.grad_covariance_cutoff}
    #
    # return out_data_dict

    if not isinstance(desc_data_mono, quippy.descriptors_module.descriptor_data_mono):
        raise TypeError('Not descriptor_data_mono given')

    # TODO: only take the ones actually needed, this is good for debuggin now though
    out_data_dict = dict()
    # 'has_grad_data':
    # 'grad_data': desc_data_mono.grad_data,

    try:
        out_data_dict['has_grad_data'] = desc_data_mono.has_grad_data
        #
    except:
        # could be None too
        # out_data_dict['has_grad_data'] = None
        # out_data_dict['grad_data'] = None
        # pass
        pass
        #print('failed has_grad_data')

    try:
        out_data_dict['grad_data'] = desc_data_mono.grad_data
    except:
        pass
        #print('failed grad_data')

    try:
        out_data_dict['ii'] = desc_data_mono.ii
        # 'ii': desc_data_mono.ii,
    except:
        pass
        #print('failed ii')

    try:
        out_data_dict['pos'] = desc_data_mono.pos
        # 'pos': desc_data_mono.pos,
    except:
        pass
        #print('failed pos')

    try:
        out_data_dict['grad_covariance_cutoff'] = desc_data_mono.grad_covariance_cutoff
        # 'grad_covariance_cutoff': desc_data_mono.grad_covariance_cutoff,
    except:
        pass
        #print('failed grad_covariance_cutoff')

    try:
        out_data_dict['covariance_cutoff'] = desc_data_mono.covariance_cutoff
        # 'covariance_cutoff': desc_data_mono.covariance_cutoff,
    except:
        pass
        #print('failed covariance_cutoff')

    try:
        out_data_dict['data'] = desc_data_mono.data
    except:
        # 'data': desc_data_mono.data,
        pass
        #print('failed data')

    try:
        out_data_dict['has_data'] = desc_data_mono.has_data
    except:
        # 'data': desc_data_mono.has_data,
        pass
        #print('failed has_data')

    try:
        out_data_dict['ci'] = desc_data_mono.ci
        # 'ii': desc_data_mono.ci,
    except:
        pass
        #print('failed ci')

    return out_data_dict

# for a try at a nicer implementation:
# def exception(function):
#     """
#     A decorator that wraps the passed in function and logs
#     exceptions should one occur
#     """
#
#     @wraps(function)
#     def wrapper(*args, **kwargs):
#         try:
#             return function(*args, **kwargs)
#         except:
#            pass
#
#     return wrapper
#
#

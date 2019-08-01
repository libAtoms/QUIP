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
import inspect
from copy import deepcopy as cp

import numpy as np
import ase

import f90wrap.runtime

import _quippy
import quippy


__all__ = ['ase_to_quip', 'descriptor_data_mono_to_dict', 'velocities_ase_to_quip', 'velocities_quip_to_ase', 'set_doc']


# conversion between ase and quip mass, taken from Fortran source
MASSCONVERT = 103.6426957074462


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
            # raise an error for the wrong object given
            raise TypeError('quip_atoms argument is not of valid type, cannot work with it')

    else:
        # need to regenerate the quip atoms object
        quip_atoms = quippy.atoms_types_module.Atoms(len(ase_atoms), ase_atoms.get_cell().transpose())

    quip_atoms.pos[:] = ase_atoms.get_positions().T.copy()
    quip_atoms.is_periodic[:] = ase_atoms.get_pbc()   # fixme this is not making sure it is a quip compatible arr.
    quip_atoms.z[:] = ase_atoms.numbers

    if ase_atoms.has('momenta'):
        # if ase atoms has momenta then add velocities to the quip object
        # workaround for the interfaces not behaving properly in the wrapped code, see f90wrap issue #86
        _quippy.f90wrap_atoms_add_property_real_2da(this=quip_atoms._handle, name='velo',
                                                    value=velocities_ase_to_quip(ase_atoms.get_velocities()))

    # go through all properties
    return quip_atoms


def velocities_ase_to_quip(velocities):
    """
    Convert the ASE velocities to QUIP velocities

    :param velocities: velocities obtained from ase, with Atoms.get_velocities()
    :return:
    """

    return (velocities / np.sqrt(MASSCONVERT)).T


def velocities_quip_to_ase(velocities):
    """
    Convert the QUIP velocities to ASE velocities

    :param velocities: velocities obtained from quip, with quip_atom.velo[:]
    :return:
    """

    return (velocities * np.sqrt(MASSCONVERT)).T


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
        # print('failed has_grad_data')

    try:
        out_data_dict['grad_data'] = desc_data_mono.grad_data
    except:
        pass
        # print('failed grad_data')

    try:
        out_data_dict['ii'] = desc_data_mono.ii
        # 'ii': desc_data_mono.ii,
    except:
        pass
        # print('failed ii')

    try:
        out_data_dict['pos'] = desc_data_mono.pos
        # 'pos': desc_data_mono.pos,
    except:
        pass
        # print('failed pos')

    try:
        out_data_dict['grad_covariance_cutoff'] = desc_data_mono.grad_covariance_cutoff
        # 'grad_covariance_cutoff': desc_data_mono.grad_covariance_cutoff,
    except:
        pass
        # print('failed grad_covariance_cutoff')

    try:
        out_data_dict['covariance_cutoff'] = desc_data_mono.covariance_cutoff
        # 'covariance_cutoff': desc_data_mono.covariance_cutoff,
    except:
        pass
        # print('failed covariance_cutoff')

    try:
        out_data_dict['data'] = desc_data_mono.data
    except:
        # 'data': desc_data_mono.data,
        pass
        # print('failed data')

    try:
        out_data_dict['has_data'] = desc_data_mono.has_data
    except:
        # 'data': desc_data_mono.has_data,
        pass
        # print('failed has_data')

    try:
        out_data_dict['ci'] = desc_data_mono.ci
        # 'ii': desc_data_mono.ci,
    except:
        pass
        # print('failed ci')

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



def get_dict_arrays(fdict):
    """Takes the arrays from a quippy dictionary. Copies.

    Probably fails if there are non-array elements in the dictionary"""

    if not isinstance(fdict, quippy.dictionary_module.Dictionary):
        raise TypeError('fdict argument is not a quippy.dictionary_module.Dictionary')

    arrays = {}
    for i in range(1, fdict.n + 1):
        key = fdict.get_key(i)
        key = key.strip().decode('ascii')
        # fixme: fails for non_array elements. Make universal: compatible with array or scalar content in dictionary
        try:    # this is an unsufficient temporary fix
            value = f90wrap.runtime.get_array(f90wrap.runtime.sizeof_fortran_t,
                                              fdict._handle, _quippy.f90wrap_dictionary__array__, key)
            arrays[key] = value.copy()
        except ValueError:
            value = fdict.get_value(key)
            try:
                # normally it is an tuple, because the error arf from fortran is converted to output
                arrays[key] = cp(value[0])
            except TypeError:
                arrays[key] = cp(value)

    return arrays


def set_doc(doc, extra):
    def wrap(method):
        method.__doc__ = update_doc_string(doc, extra)
        return method
    return wrap

def update_doc_string(doc, extra, sections=None, signature=None):
    """
    Insert `extra` in the docstring `doc`, before the first matching section

    Searches for each section heading in the list `sections` in turn.
    If sections is not given, the default is `['Parameters', 'See also']`.
    If not sections are found, extra text is appended to end of docstring.
    """

    if sections is None:
       sections = ['Parameters', 'See also']

    try:
       doc = inspect.cleandoc(doc)
       extra = inspect.cleandoc(extra)
    except AttributeError:
       pass

    extra = '\n' + extra + '\n'

    lines = doc.split('\n')

    if signature is not None:
        lines[0] = signature

    for section in sections:
       indices = [i for i, line in enumerate(lines) if line == section ]
       if len(indices) == 1:
          break
    else:
        indices = [len(lines)-1] # insert at end

    index, = indices
    doc = '\n'.join([line.rstrip() for line in lines[:index] + extra.split('\n') + lines[index:]])
    doc = doc.replace('\n\n\n', '\n\n')

    return doc

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


"""
Conversions between ase and fortran atoms objects
"""
import inspect
from copy import deepcopy as cp

import quippy._quippy as _quippy
import ase
import f90wrap.runtime
import numpy as np
import quippy

__all__ = ['ase_to_quip', 'descriptor_data_mono_to_dict', 'velocities_ase_to_quip', 'velocities_quip_to_ase', 'set_doc']

# conversion between ase and quip mass, taken from Fortran source
MASSCONVERT = 103.6426957074462


def ase_to_quip(ase_atoms: ase.Atoms, quip_atoms=None, add_arrays=None, add_info=None):
    """
    Converter to put the info from an ase atoms object into a quip atoms object.
    Copies everything to make sure there is not linking back.

    Notes on add_arrays and add_info:
        - overwriting a parameter is not possible yet
        - only float arrays can be added, integers are converted to floats by fortran, fails for strings
        - keys can only be strings, as the fortran dictionary will not accept anything else,\
        integer keys are converted to strings
        - possible types:
            None - only defaults: pos, Z, cell, pbc, momenta (if exists)
            str  - single key
            list - all the list elements are added
            True - all of the arrays

    Copying from ase.Atoms.info supports:
    =============== ===== ===== ==== =======================================================================
    type            0D    1D    2D   Note
    =============== ===== ===== ==== =======================================================================
    real            no*   yes   yes  scalars kept as arrays of one dimension with one element
    int             yes   yes   yes
    bool            no*   yes   NO   zero dim. kept as array of one dimension with one element
    =============== ===== ===== ==== =======================================================================

    Copying from ase.Atoms.arrays supports:
    =============== ===== ==== =======================================================================
    type            1D    2D   Note
    =============== ===== ==== =======================================================================
    real            yes   yes  scalars kept as arrays of one dimension with one element
    int             yes   yes
    bool            yes   NO   zero dim. kept as array of one dimension with one element
    =============== ===== ==== =======================================================================

    Logical arrays are kept as ones and zeros, and normally retrieved as int32.

    :param ase_atoms:
    :param quip_atoms:
    :param add_arrays: keys to take from ase.Atoms.arrays
    :param add_info:  keys to take from ase.Atoms.info
    :return:
    """

    lattice = ase_atoms.get_cell().T.copy()
    if quip_atoms is not None:
        if isinstance(quip_atoms, quippy.atoms_types_module.Atoms):
            # check if the length matches, otherwise make a new one in place of that
            if len(ase_atoms) != quip_atoms.n:
                # need to regenerate the quip atoms object
                quip_atoms = quippy.atoms_types_module.Atoms(len(ase_atoms), lattice)
            else:
                # but the cell needs to be set anyways
                quip_atoms.set_lattice(lattice, scale_positions=False)
        else:
            # raise an error for the wrong object given
            raise TypeError('quip_atoms argument is not of valid type, cannot work with it')

    else:
        # need to regenerate the quip atoms object
        quip_atoms = quippy.atoms_types_module.Atoms(len(ase_atoms), lattice)

    quip_atoms.pos[:] = ase_atoms.get_positions().T.copy()
    quip_atoms.is_periodic[:] = ase_atoms.get_pbc()
    quip_atoms.z[:] = ase_atoms.numbers
    quip_atoms.set_atoms(quip_atoms.z)  # set species and mass

    if ase_atoms.has('momenta'):
        # if ase atoms has momenta then add velocities to the quip object
        # workaround for the interfaces not behaving properly in the wrapped code, see f90wrap issue #86
        _quippy.f90wrap_atoms_add_property_real_2da(this=quip_atoms._handle, name='velo',
                                                    value=velocities_ase_to_quip(ase_atoms.get_velocities()))

    def key_spec_to_list(keyspec, default, exclude=()):
        if keyspec is True:
            # taking all the array keys that are not handled elsewhere
            keyspec = set(default.keys())
            [keyspec.discard(used_key) for used_key in exclude]
            keyspec = list(keyspec)
        elif isinstance(keyspec, str):
            # if only one is given as a string
            keyspec = [keyspec]
        elif isinstance(keyspec, list) or isinstance(keyspec, np.ndarray):
            keyspec = list(keyspec)
        else:
            # fixme: decide what to do here, now it is just not adding anything
            keyspec = []

        return keyspec

    # go through all properties for issue#170
    if add_arrays is not None:
        add_arrays = key_spec_to_list(add_arrays, ase_atoms.arrays, exclude=['numbers', 'positions', 'momenta'])
        for info_name in add_arrays:
            try:
                value = np.array(ase_atoms.arrays[info_name])
            except KeyError:
                # fixme: give some warning here if needed
                continue
            add_property_array(quip_atoms, info_name, value)

    if add_info is not None:
        add_info = key_spec_to_list(add_info, ase_atoms.info, exclude=[])
        for info_name in add_info:
            try:
                value = np.array(ase_atoms.info[info_name])
            except KeyError:
                # fixme: give some warning here if needed
                continue
            add_param_value(quip_atoms, info_name, value)

    return quip_atoms


def add_param_value(quip_atoms, name, value):
    """
    Adds property to a quip atoms object in params dictionary, from ase.Atoms.info

    Supports:
    =============== ===== ===== ==== =======================================================================
    type            0D    1D    2D   Note
    =============== ===== ===== ==== =======================================================================
    real            no*   yes   yes  scalars kept as arrays of one dimension with one element
    int             yes   yes   yes
    bool            no*   yes   NO   zero dim. kept as array of one dimension with one element
    =============== ===== ===== ==== =======================================================================

    Note: logical arrays are kept as ones and zeros, and normally retrieved as int32.

    :param quip_atoms:
    :param name:
    :param value:
    :return:
    """
    # to make sure it is a numpy array, so we can use the dtype and shape of it
    if not isinstance(value, np.ndarray):
        value = np.array(value)

    # decide type
    arr_dtype_kind = value.dtype.kind
    dim = len(value.shape)
    if arr_dtype_kind == 'b':
        # only 1D works, 0->1
        fortran_type_name = 'l'
        if dim >= 2:
            raise TypeError('2d logical array is not supported')
        elif dim == 0:
            # change to a one dimension instead
            value = np.atleast_1d(value)
            dim = len(value.shape)  # update dimension!!!
    elif arr_dtype_kind in ['u', 'i']:
        # 0,1,2 are all working
        fortran_type_name = 'i'
    elif arr_dtype_kind == 'f':
        # only 1,2 D works, 0->1
        fortran_type_name = 'r'
        if dim == 0:
            # change to a one dimension instead
            value = np.atleast_1d(value)
            dim = len(value.shape)  # update dimension!!!
    else:
        # so it is one of:
        # c complex floating - point
        # m timedelta
        # M datetime
        # O object
        # V void
        # strings not supported yet, f90wrap_atoms_add_property_str needs some
        raise TypeError('given dtype ({}) is not supported'.format(arr_dtype_kind))

    # decide dim
    if dim == 0:
        add_property_method = getattr(_quippy, 'f90wrap_dictionary_set_value_{}'.format(fortran_type_name))
    elif dim == 1:
        add_property_method = getattr(_quippy, 'f90wrap_dictionary_set_value_{}_a'.format(fortran_type_name))
    elif dim == 2:
        add_property_method = getattr(_quippy, 'f90wrap_atoms_add_property_{}_2da'.format(fortran_type_name))
        value = value.T
    else:
        raise ValueError(
            'unsupported dimension ({}) of attribute in conversion from ase to quip atoms objects'.format(dim))
    add_property_method(this=quip_atoms.params._handle, key=name, value=value)


def add_property_array(quip_atoms, name, value):
    """
    Adds property to a quip atoms object

    Supports:
    =============== ===== ==== =======================================================================
    type            1D    2D   Note
    =============== ===== ==== =======================================================================
    real            yes   yes  scalars kept as arrays of one dimension with one element
    int             yes   yes
    bool            yes   NO   zero dim. kept as array of one dimension with one element
    =============== ===== ==== =======================================================================

    Note: logical arrays are kept as ones and zeros, and normally retrieved as int32.

    :param quip_atoms:
    :param name:
    :param value:
    :return:
    """
    # to make sure it is a numpy array, so we can use the dtype and shape of it
    if not isinstance(value, np.ndarray):
        value = np.array(value)

    # add the value, 1d/2d array
    dim = len(value.shape)
    arr_dtype_kind = value.dtype.kind

    # decide the fortran type
    if arr_dtype_kind == 'b':
        if dim < 2:
            fortran_type_name = 'logical'
        else:
            raise TypeError('2d logical array is not supported')
    elif arr_dtype_kind in ['u', 'i']:
        fortran_type_name = 'int'
    elif arr_dtype_kind == 'f':
        fortran_type_name = 'real'
    # elif arr_dtype_kind in ['S', 'U']:
    #     fortran_type_name = 'str'
    else:
        # so it is one of:
        # c complex floating - point
        # m timedelta
        # M datetime
        # O object
        # V void
        # strings not supported yet, f90wrap_atoms_add_property_str needs some
        raise TypeError('given dtype ({}) is not supported'.format(arr_dtype_kind))

    # decide dim
    if dim == 1:
        add_property_method = getattr(_quippy, 'f90wrap_atoms_add_property_{}_a'.format(fortran_type_name))
        add_property_method(this=quip_atoms._handle, name=name, value=value)
    elif dim == 2:
        add_property_method = getattr(_quippy, 'f90wrap_atoms_add_property_{}_2da'.format(fortran_type_name))
        add_property_method(this=quip_atoms._handle, name=name, value=value.T)
    else:
        raise ValueError(
            'unsupported dimension ({}) of attribute in conversion from ase to quip atoms objects'.format(dim))


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

    if not isinstance(desc_data_mono, quippy.descriptors_module.descriptor_data_mono):
        raise TypeError('Not descriptor_data_mono given')

    out_data_dict = dict()

    def take_value(key):
        """
        Take the arg if it exists
        """

        try:
            out_data_dict[key] = getattr(desc_data_mono, key)
        except AttributeError:
            pass
        except ValueError:
            pass

    # fixme: only take the ones actually needed, this is good for debugging now though
    for key in ['has_grad_data', 'ii', 'pos', 'grad_covariance_cutoff', 'covariance_cutoff', 'data', 'has_data',
                'grad_data', 'ci']:
        take_value(key)

    return out_data_dict


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
        try:  # this is an unsufficient temporary fix
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
        indices = [i for i, line in enumerate(lines) if line == section]
        if len(indices) == 1:
            break
    else:
        indices = [len(lines) - 1]  # insert at end

    index, = indices
    doc = '\n'.join([line.rstrip() for line in lines[:index] + extra.split('\n') + lines[index:]])
    doc = doc.replace('\n\n\n', '\n\n')

    return doc

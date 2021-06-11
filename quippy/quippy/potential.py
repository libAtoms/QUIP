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
ASE-compatible Calculator from a quip-potential object

"""

from copy import deepcopy as cp

import ase
import ase.calculators.calculator
import numpy as np
import quippy
from ase.io.extxyz import key_val_dict_to_str
from quippy.convert import set_doc

__all__ = ['Potential']


@set_doc(quippy.potential_module.__doc__,
"""
Pythonic interface to auto-generated quippy.potential_module.Potential class
""")
class Potential(ase.calculators.calculator.Calculator):
    callback_map = {}

    # this list only includes standard ASE properties. Additional non-standard
    # results are returned in the `extra_results` dictionary.
    implemented_properties = ['energy', 'free_energy', 'forces', 'stress',
                              'stresses', 'energies']
    
    extra_properties = ['virial', 'local_virial', 'local_energy']
    
    @set_doc(quippy.potential_module.Potential.__init__.__doc__,
    """
    ------------------------------------------------------------------------
    from old quippy arguments:

    not implemented yet:
        bulk_scale=None
        mpi_obj=None
        callback=None
        finalise=True
    """)
    def __init__(self, args_str="",
                 pot1=None, pot2=None,
                 param_str=None,
                 param_filename=None,
                 atoms=None,
                 calculation_always_required=False, calc_args=None,
                 add_arrays=None, add_info=None, **kwargs):
        quippy.potential_module.Potential.__init__.__doc__

        self._default_properties = ['energy', 'forces']
        self.calculation_always_required = calculation_always_required

        ase.calculators.calculator.Calculator.__init__(self, restart=None, label=None,
                                                       atoms=atoms, directory='.', **kwargs)
        # init the quip potential
        if param_filename is not None and isinstance(param_filename, str):
            # from a param filename
            self._quip_potential = quippy.potential_module.Potential.filename_initialise(args_str=args_str,
                                                                                         param_filename=param_filename)
        elif pot1 is not None and pot2 is not None:
            # from sum of two potentials
            # noinspection PyProtectedMember
            self._quip_potential = quippy.potential_module.Potential(
                args_str=args_str,
                pot1=(pot1._quip_potential if isinstance(pot1, quippy.potential.Potential) else pot1),
                pot2=(pot2._quip_potential if isinstance(pot2, quippy.potential.Potential) else pot2))
        else:
            # from a param string
            self._quip_potential = quippy.potential_module.Potential(args_str=args_str, param_str=param_str)

        # init the quip atoms as None, to have the variable
        self._quip_atoms = None
        # init the info and array keys that need to be added when converting atoms objects
        self.add_arrays = add_arrays
        self.add_info = add_info

        # from old
        if atoms is not None:
            atoms.calc = self
        self.name = args_str
        if isinstance(calc_args, dict):
            calc_args = key_val_dict_to_str(calc_args)
        elif calc_args is None:
            calc_args = ""
        self.calc_args = calc_args
        
        # storage of non-standard results
        self.extra_results = {'config': {},
                              'atoms': {}}

    @set_doc(quippy.potential_module.Potential.calc.__doc__,
    """
    Pythonic wrapper to `quippy.potential_module.Potential.calc()`

    atoms: ase.atoms.Atoms object
        Atoms object with which to calculate. This is converted to
        a `quippy.atoms_module.Atoms` object to pass to Fortran.
    properties: list of str
        List of what needs to be calculated.  Can be any combination
        of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
        and 'magmoms'.
    system_changes: list of str
        List of what has changed since last calculation.  Can be
        any combination of these six: 'positions', 'numbers', 'cell',
        'pbc', 'initial_charges' and 'initial_magmoms'.
    calc_args: argument string to pass to Fortran calc() routine.
        This is appended to `self.calc_args`, if this is set.
    add_arrays: list, str or True
        Arrays to carry over from `atoms.arrays` to the Fortran atoms object when calling calculate()
        If the same argument is filled in calculate, then it will overwrite this.
        String arrays not supported.
        possible types:
            None - only defaults: `pos`, `Z`, `cell`, `pbc`, `momenta` (as velo, if exists)
            str  - single key
            list - all list elements as keys
            True - all of the arrays from `atoms.arrays`
    add_info: list, str or True
        Info to carry over from `atoms.info` to the Fortran atoms object when calling calculate()
        If the same argument is filled in calculate, then it will overwrite this.
        possible types:
            None - only defaults: `pos`, `Z`, `cell`, `pbc`, `momenta` (as velo, if exists)
            str  - single key
            list - all list elements as keys
            True - all of the values from `atoms.info`

    Arrays can also be passed directly to the Fortran routine using
    the `forces`, `virial`, `local_energy`, `local_virial`
    arguments. These should be pre-allocated as Fortran-contigous arrays,
    e.g. `forces = np.zeros((len(atoms), 3), order='F')`.

    If present, `vol_per_atom` is used to convert from local_virial to per-atom
    stresses; this can either be a scalar or the name of an array in
    `atoms.arrays`. If not present the average volume per atom is used.

    Additional keyword arguments are appended to `calc_args`.
    
    Non-standard calculation results (such as `gap_local_variance`) are stored in 
    `self.extra_results`.
    
    """)
    def calculate(self, atoms=None, properties=None, system_changes=None,
                  forces=None, virial=None, local_energy=None,
                  local_virial=None, vol_per_atom=None,
                  calc_args=None, add_arrays=None,
                  add_info=None, **kwargs):

        # handling the property inputs
        if properties is None:
            properties = self.get_default_properties()
        else:
            properties = list(set(self.get_default_properties() + properties))

        if len(properties) == 0:
            raise RuntimeError('Nothing to calculate')

        for prop in properties:
            if prop not in self.implemented_properties + self.extra_properties:
                raise RuntimeError("Don't know how to calculate property '%s'" % prop)

        # initialise dictionary to arguments to be passed to calculator
        _dict_args = {}
        val = _check_arg(forces)
        if val == 'y':
            properties += ['force']
        elif val == 'add':
            properties += ['force']
            _dict_args['force'] = forces

        val = _check_arg(virial)
        if val == 'y':
            properties += ['virial']
        elif val == 'add':
            properties += ['virial']
            _dict_args['virial'] = virial

        val = _check_arg(local_energy)
        if val == 'y':
            properties += ['local_energy']
        elif val == 'add':
            properties += ['local_energy']
            _dict_args['local_energy'] = local_energy

        val = _check_arg(local_virial)
        if val == 'y':
            properties += ['local_virial']
        elif val == 'add':
            properties += ['local_virial']
            _dict_args['local_virial'] = local_virial

        # check if a calculation is needed - this must be done before we update self.atoms
        if not self.calculation_always_required and not self.calculation_required(atoms, properties):
            return
        
        # call the base class method, which updates self.atoms to atoms
        ase.calculators.calculator.Calculator.calculate(self, atoms, properties, system_changes)

        self.extra_results = {'config': {},
                              'atoms': {}} # reset all extra results on each new calculation       

        # construct the quip atoms object which we will use to calculate on
        # if add_arrays/add_info given to this object is not None, then OVERWRITES the value set in __init__
        self._quip_atoms = quippy.convert.ase_to_quip(self.atoms,
                                                      add_arrays=add_arrays if add_arrays is not None else self.add_arrays,
                                                      add_info=add_info if add_info is not None else self.add_info)

        # constructing args_string with automatically aliasing the calculateable non-quippy properties
        # calc_args string to be passed to Fortran code
        args_str = self.calc_args
        if calc_args is not None:
            if isinstance(calc_args, dict):
                calc_args = key_val_dict_to_str(calc_args)
            args_str += ' ' + calc_args
        if kwargs is not None:
            args_str += ' ' + key_val_dict_to_str(kwargs)

        args_str += ' energy'
        # no need to add logic to energy, it is calculated anyways (returned when potential called)
        if 'virial' in properties or 'stress' in properties:
            args_str += ' virial'
        if 'local_virial' in properties or 'stresses' in properties:
            args_str += ' local_virial'
        if 'energies' in properties or 'local_energy' in properties:
            args_str += ' local_energy'
        if 'forces' in properties:
            args_str += ' force'
        # TODO: implement 'elastic_constants', 'unrelaxed_elastic_constants', 'numeric_forces'

        # fixme: workaround to get the calculated energy, because the wrapped dictionary is not handling that float well
        ener_dummy = np.zeros(1, dtype=float)

        # the calculation itself
        # print('Calling QUIP Potential.calc() with args_str "{}"'.format(args_str))
        self._quip_potential.calc(self._quip_atoms, args_str=args_str, energy=ener_dummy, **_dict_args)

        # retrieve data from _quip_atoms.properties and _quip_atoms.params
        _quip_properties = quippy.convert.get_dict_arrays(self._quip_atoms.properties)
        _quip_params = quippy.convert.get_dict_arrays(self._quip_atoms.params)

        self.results['energy'] = ener_dummy[0]
        self.results['free_energy'] = self.results['energy']

        # process potential output to ase.properties
        # not handling energy here, because that is always returned by the potential above
        if 'virial' in _quip_params.keys():
            stress = -_quip_params['virial'].copy() / self.atoms.get_volume()
            # convert to 6-element array in Voigt order
            self.results['stress'] = np.array([stress[0, 0], stress[1, 1], stress[2, 2],
                                               stress[1, 2], stress[0, 2], stress[0, 1]])
            self.extra_results['config']['virial'] = _quip_params['virial'].copy()

        if 'force' in _quip_properties.keys():
            self.results['forces'] = np.copy(_quip_properties['force'].T)

        if 'local_energy' in _quip_properties.keys():
            self.results['energies'] = np.copy(_quip_properties['local_energy'])
            self.extra_results['atoms']['local_energy'] = np.copy(_quip_properties['local_energy'])

        if 'local_virial' in _quip_properties.keys():
            self.extra_results['atoms']['local_virial'] = np.copy(_quip_properties['local_virial'])

        if 'stresses' in properties:
            # use the correct atomic volume
            if vol_per_atom is not None:
                if vol_per_atom in self.atoms.arrays.keys():
                    # case of reference to a column in atoms.arrays
                    _v_atom = self.atoms.arrays[vol_per_atom]
                else:
                    # try for case of a given volume
                    try:
                        _v_atom = float(vol_per_atom)
                    except ValueError:
                        # cannot convert to float, so wrong
                        raise ValueError('volume_per_atom: not found in atoms.arrays.keys() and cannot utilise value '
                                         'as given atomic volume')
            else:
                # just use average
                _v_atom = self.atoms.get_volume() / self._quip_atoms.n
            self.results['stresses'] = -np.copy(_quip_properties['local_virial']).T.reshape((self._quip_atoms.n, 3, 3),
                                                                                            order='F') / _v_atom

        # all non-standard results now go in self.extra_results
        _skip_keys = set(list(self.results.keys()) + ['Z', 'pos', 'species',
                                                      'map_shift', 'n_neighb',
                                                      'force', 'local_energy',
                                                      'local_virial', 'velo'])
        # any other params (per-config properties)
        for param, val in _quip_params.items():
            if param not in _skip_keys:
                self.extra_results['config'][param] = cp(val)

        # any other arrays (per-atom properties)
        for prop, val in _quip_properties.items():
            if prop not in _skip_keys:
                # transpose before copying because of setting `order=C` here; issue#151
                self.extra_results['atoms'][prop] = np.copy(val.T, order='C')

    def get_virial(self, atoms=None):
        self.get_stress(atoms)
        return self.extra_results['config']['virial']

    def get_local_virial(self, atoms=None):
        self.get_stresses(atoms)
        return self.extra_results['atoms']['local_virial']

    def get_local_energy(self, atoms=None):
        return self.get_energies(atoms)

    def get_energies(self, atoms=None):
        return self.get_property('energies', atoms)

    def get_stresses(self, atoms=None):
        return self.get_property('stresses', atoms)

    def get_default_properties(self):
        """Get the list of properties to be calculated by default"""
        return self._default_properties[:]

    def set_default_properties(self, properties):
        """Set the list of properties to be calculated by default"""
        self._default_properties = properties[:]


def _check_arg(arg):
    """Checks if the argument is True bool or string meaning True"""

    true_strings = ('True', 'true', 'T', 't', '.true.', '.True.')

    if arg is None:
        return 'n'
    else:
        if isinstance(arg, bool):
            if arg:
                return 'y'
            else:
                return 'n'
        if isinstance(arg, str):
            if arg in true_strings:
                return 'y'
            else:
                return 'n'
        return 'add'

# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
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
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import quippy
import numpy as np
from quippy.util import quip_xml_parameters
from quippy.elastic import stress_matrix

def calculator_callback_factory(calculator):
    """Return a Python function which can be used as a quippy
       Potential callback routine for the given ASE calculator."""

    def callback(at):
        at.set_calculator(calculator)
        if at.calc_energy:
            at.params['energy'] = at.get_potential_energy()
        if at.calc_local_e:
            at.add_property('local_e', 0.0, overwrite=True)
            at.local_e[:] = at.get_potential_energies()
        if at.calc_force:
            at.add_property('force', 0.0, n_cols=3, overwrite=True)
            at.force[:] = at.get_forces().T
        if at.calc_virial:
            stress = at.get_stress()
            virial = quippy.fzeros((3,3))
            virial[:,:] = stress_matrix(-stress*at.get_volume()/quippy.GPA)
            at.params['virial'] = virial
        if at.calc_local_virial:
            stresses = at.get_stresses()
            at.add_property('local_virial', 0.0, n_cols=9, overwrite=True)
            lv = at.local_virial.view(np.ndarray)
            vol = at.get_volume()
            for i, stress in enumerate(stresses):
                lv[:,i] = -(stress*vol/quippy.GPA).reshape((9,), order='F')

    return callback

class Potential(quippy.FortranPotential):
    """Potential class which abstracts all QUIP interatomic potentials.

       Provides interface to all energy/force/virial
       calculating schemes, including actual calculations, as well as
       abstract hybrid schemes such as LOTF, Force Mixing, ONIOM,
       and the local energy scheme.

       It can also be used to geometry optimise an Atoms structure,
       using the 'minim' method."""

    __doc__ = quippy.FortranPotential.__doc__
    callback_map = {}

    def __init__(self, args_str=None, pot1=None, pot2=None, param_str=None,
                 param_filename=None, bulk_scale=None, mpi_obj=None,
                 callback=None, calculator=None,
                 fortran_indexing=True, fpointer=None, finalise=True,
                 error=None):
        """Typically a Potential is constructed from an initialisation
        args_str and an XML parameter file, e.g.::

            pot = Potential('IP SW', param_filename='SW.xml')

        creates a Stillinger-Weber potential using the paramters from
        the file 'SW.xml'.

        The Potential class also implements the ASE calculator
        interface: see the get_forces(), get_stress(), get_stresses(),
        get_potential_energy(), get_potential_energies() methods.

        Moreover, it's possible to do this the other way round and use
        an ASE calculator as a QUIP potential by passing it as the
        `calculator` argument to the Potential constructor, e.g.::

            from ase.calculators.morse import MorsePotential
            pot = Potential(calculator=MorsePotential)
        """

        self.atoms = None
        self.energy = None
        self.energies = None
        self.forces = None
        self.stress = None
        self.stresses = None
        self.numeric_forces = None

        if callback is not None or calculator is not None:
            if args_str is None:
                args_str = 'callbackpot'

        if param_filename is not None:
            param_str = open(param_filename).read()

        if args_str is None and param_str is None:
            raise ValueError('Need one of args_str,param_str,param_filename')

        if args_str is None:
            xml_label = param_str[:param_str.index('\n')].translate(None,'<>')
            args_str = 'xml_label=%s' % xml_label

        if args_str.lower().startswith('callbackpot'):
            if not 'label' in args_str:
                args_str = args_str + ' label=%d' % id(self)
        else:
            # if param_str missing, try to find default set of QUIP params
            if param_str is None and pot1 is None and pot2 is None:
                param_str = quip_xml_parameters(args_str)

        quippy.FortranPotential.__init__(self, args_str, pot1=pot1, pot2=pot2,
                                         param_str=param_str,
                                         bulk_scale=bulk_scale,
                                         mpi_obj=mpi_obj,
                                         fortran_indexing=fortran_indexing,
                                         fpointer=fpointer, finalise=finalise,
                                         error=error)

        if args_str.lower().startswith('callbackpot'):
            quippy.FortranPotential.set_callback(self, Potential.callback)

            if callback is not None:
                self.set_callback(callback)

            if calculator is not None:
                self.set_callback(calculator_callback_factory(calculator))

    __init__.__doc__ = quippy.FortranPotential.__init__.__doc__

    @staticmethod
    def callback(at_ptr):
        from quippy import Atoms
        at = Atoms(fpointer=at_ptr, finalise=False)
        if at.params['label'] not in Potential.callback_map:
            raise ValueError('Unknown Callback label %s' % at.params['label'])
        Potential.callback_map[at.params['label']](at)

    def set_callback(self, callback):
        Potential.callback_map[str(id(self))] = callback

    def update(self, atoms):
        if self.atoms is not None and self.atoms.equivalent(atoms):
            return

        # store a copy of `atoms` as a quippy.Atoms instance
        self.atoms = quippy.Atoms(atoms)

        # build neighbour list
        self.atoms.set_cutoff(self.cutoff())
        self.atoms.calc_connect()

        # mark all quantities as needing to be recalculated
        self.energy = None
        self.energies = None
        self.forces = None
        self.stress = None
        self.stresses = None
        self.numeric_forces = None

    # Synonyms for `update` for compatibility with ASE calculator interface
    initialize = update
    set_atoms  = update

    def calculation_required(self, atoms, quantities):
        self.update(atoms)
        for quantity in quantities:
            if getattr(self, quantity) is None:
                return True
        return False

    def calculate(self, atoms, quantities=None):
        if quantities is None:
            quantities = ['energy', 'forces', 'stress']

        if len(quantities) == 0:
            raise RuntimeError('Nothing to calculate')

        if not self.calculation_required(atoms, quantities):
            return

        args_map = {'energy':          'energy',
                    'energies':        'local_energy',
                    'forces':          'force',
                    'stress':          'virial',
                    'numeric_forces':  'force=numeric_force force_using_fd=T force_fd_delta=1.0e-5',
                    'stresses':        'local_virial'}

        args_str = ' '.join([args_map[quantity] for quantity in quantities])
        print 'calling Potential.calc with args_str %s' % args_str
        self.calc(self.atoms, args_str=args_str)

        if 'energy' in quantities:
            self.energy = float(self.atoms.energy)
        if 'energies' in quantities:
            self.energies = self.atoms.local_energy.view(np.ndarray)
        if 'forces' in quantities:
            self.forces = self.atoms.force.view(np.ndarray).T
        if 'numeric_forces' in quantities:
            self.numeric_forces = self.atoms.numeric_force.view(np.ndarray).T
        if 'stress' in quantities:
            self.stress = -(self.atoms.virial.view(np.ndarray)*
                            quippy.GPA/self.atoms.get_volume())
        if 'stresses' in quantities:
            self.stresses = np.zeros((len(self.atoms), 3, 3))

            lv = self.atoms.local_virial.view(np.ndarray)
            vol = self.atoms.get_volume()
            for i in range(len(self.atoms)):
                self.stresses[i,:,:] = -(lv[:,i].reshape((3,3),order='F')*
                                         quippy.GPA/vol)

    def get_potential_energy(self, atoms):
        self.calculate(atoms, ['energy'])
        return self.energy

    def get_potential_energies(self, atoms):
        self.calculate(atoms, ['energies'])
        return self.energies.copy()

    def get_forces(self, atoms):
        self.calculate(atoms, ['forces'])
        return self.forces.copy()

    def get_numeric_forces(self, atoms):
        self.calculate(atoms, ['numeric_forces'])
        return self.numeric_forces.copy()

    def get_stress(self, atoms):
        self.calculate(atoms, ['stress'])
        return self.stress.copy()

    def get_stresses(self, atoms):
        self.calculate(atoms, ['stresses'])
        return self.stresses.copy()

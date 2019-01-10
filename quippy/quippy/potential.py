"""
ASE-compatible Calculator from a fortran-potential object

"""

import ase
import ase.calculators.calculator
import numpy as np

import _quippy
import quippy


class potential(ase.calculators.calculator.Calculator):

    callback_map = {}

    implemented_properties = ['energy', 'forces']
    # earlier in quippy
    # ['energy', 'energies', 'forces', 'stress', 'stresses',
    #  'numeric_forces', 'elastic_constants',
    #  'unrelaxed_elastic_constants']

    def __init__(self, args_str, param_str, atoms, **kwargs):

        # update_docstring not implemented yet, it was oo_fortran.update_doc_string() in the earlier version

        """

        from fortran_docstring:

        args_str : str
            Valid arguments are 'Sum', 'ForceMixing', 'EVB', 'Local_E_Mix' and 'ONIOM', and any type of simple_potential
        param_str : str
            contents of xml parameter file for potential initializers, if needed

        -----------------------------------------------------------------------

        ase calculator has the following arguments for initialisation:

        let's not care about files for now, so just take them as None

        nay     restart: str
                    Prefix for restart file.  May contain a directory.  Default
                    is None: don't restart.
        nay     ignore_bad_restart_file: bool
                    Ignore broken or missing restart file.  By default, it is an
                    error if the restart file is missing or broken.
        nay     label: str
                    Name used for all files.  May contain a directory.
                atoms: Atoms object
                    Optional Atoms object to which the calculator will be
                    attached.  When restarting, atoms will get its positions and
                    unit-cell updated from file.

        ------------------------------------------------------------------------
        from old quippy arguments:

        used:
            init_args=None, param_str=None, atoms=None

        not used:
            calculator=None
            fpointer=None
            error=None

        not implemented yet:
            pot1=None, pot2=None
            param_filename=None
            bulk_scale=None
            mpi_obj=None
            callback=None
            calculation_always_required=False
            finalise=True
                """

        ase.calculators.calculator.Calculator.__init__(self, restart=None, ignore_bad_restart_file=False, label=None,
                                                       atoms=atoms, **kwargs)
        # init the fortran potential
        self._fortran_potential = quippy.potential_module.Potential(args_str=args_str, param_str=param_str)
        # init the fortran atoms as None, to have the variable
        self._fortran_atoms = None

        # from old
        if atoms is not None:
            atoms.set_calculator(self)
        self.name = args_str

        pass

    def calculate(self, atoms=None, properties=('energy', 'forces'),
                  system_changes=None):
        """Do the calculation.

        properties: list of str
            List of what needs to be calculated.  Can be any combination
            of 'energy', 'forces', 'stress', 'dipole', 'charges', 'magmom'
            and 'magmoms'.
        system_changes: list of str
            List of what has changed since last calculation.  Can be
            any combination of these six: 'positions', 'numbers', 'cell',
            'pbc', 'initial_charges' and 'initial_magmoms'.

        Subclasses need to implement this, but can ignore properties
        and system_changes if they want.  Calculated properties should
        be inserted into results dictionary like shown in this dummy
        example::

            self.results = {'energy': 0.0,
                            'forces': np.zeros((len(atoms), 3)),
                            'stress': np.zeros(6),
                            'dipole': np.zeros(3),
                            'charges': np.zeros(len(atoms)),
                            'magmom': 0.0,
                            'magmoms': np.zeros(len(atoms))}

        The subclass implementation should first call this
        implementation to set the atoms attribute.
        """

        if atoms is not None:
            self.atoms = atoms.copy()
        # construct the fortran atoms object which we will use to calculate on
        quippy.convert.ase_to_fortran(self.atoms, self._fortran_atoms)

        # construct adequate arrays to put the results into
        force = np.zeros((3, self._fortran_atoms.n), order='F')

        # perform the calculation
        energy, _ferror = self._fortran_potential.calc(self._fortran_atoms, force=force)

        # store the results according to ase's standards
        self.results = {'energy': energy,
                        'forces': np.copy(force.T)
                        }


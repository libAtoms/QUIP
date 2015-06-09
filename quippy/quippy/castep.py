# HQ
import sys, string, os, operator, itertools, logging, glob, re, shutil
import numpy as np

from quippy.atoms import Atoms, make_lattice, get_lattice_params
from quippy.io import AtomsReaders, AtomsWriters, atoms_reader
from quippy.dictionary import Dictionary
from quippy.units import AU_FS, HARTREE, BOHR, BOLTZMANN_K, GPA, DEBYE
from quippy.periodictable import atomic_number

from ordereddict import OrderedDict
from math import pi
import xml.dom.minidom

__all__ = ['CastepCell', 'CastepParam', 'CastepPotential',
           'check_pspots', 'run_castep', 'read_formatted_potential',
           'read_formatted_density']
           
castep_units = {
   'ang' : 1,
   'bohr' : BOHR
   }

castep_value_map = {
   'T': True,
   'true': True,
   'F': False,
   'false': False
   }

castep_output_map = {
   True: 'true',
   False: 'false'
   }

# Valid CELL and PARAMETER keywords. Generated using get_valid_keywords() in this module with CASTEP 8.0.

valid_cell_keywords = ['lattice_cart', 'lattice_abc', 'positions_frac', 'positions_abs',
                       'symmetry_generate', 'symmetry_ops', 'symmetry_tol', 'ionic_constraints',
                       'fix_com', 'cell_constraints', 'external_pressure', 'fix_all_ions',
                       'fix_all_cell', 'species_mass', 'species_pot', 'ionic_velocities',
                       'species_lcao_states', 'kpoints_list', 'kpoints_mp_grid', 'kpoints_mp_spacing',
                       'kpoints_mp_offset', 'kpoint_list', 'kpoint_mp_grid', 'kpoint_mp_spacing',
                       'kpoint_mp_offset', 'bs_kpoint_path', 'bs_kpoint_path_spacing',
                       'bs_kpoint_list', 'bs_kpoint_mp_grid', 'bs_kpoint_mp_spacing',
                       'bs_kpoint_mp_offset', 'bs_kpoints_path', 'bs_kpoints_path_spacing',
                       'bs_kpoints_list', 'bs_kpoints_mp_grid', 'bs_kpoints_mp_spacing',
                       'bs_kpoints_mp_offset', 'phonon_supercell_matrix', 'phonon_kpoint_path',
                       'phonon_kpoint_path_spacing', 'phonon_kpoint_list', 'phonon_kpoint_mp_grid',
                       'phonon_kpoint_mp_offset', 'phonon_kpoint_mp_spacing',
                       'phonon_gamma_directions', 'phonon_kpoints_path', 'phonon_kpoints_path_spacing',
                       'phonon_kpoints_list', 'phonon_fine_kpoint_list', 'phonon_fine_kpoint_path',
                       'phonon_fine_kpoint_path_spacing', 'phonon_fine_kpoint_mp_grid',
                       'phonon_fine_kpoint_mp_spacing', 'phonon_fine_kpoint_mp_offset',
                       'optics_kpoints_list', 'optics_kpoints_mp_grid', 'optics_kpoints_mp_spacing',
                       'optics_kpoints_mp_offset', 'optics_kpoint_list', 'optics_kpoint_mp_grid',
                       'optics_kpoint_mp_spacing', 'optics_kpoint_mp_offset', 'magres_kpoint_list',
                       'magres_kpoint_path', 'magres_kpoint_path_spacing', 'magres_kpoint_mp_grid',
                       'magres_kpoint_mp_spacing', 'magres_kpoint_mp_offset', 'positions_frac_product',
                       'positions_abs_product', 'positions_frac_intermediate',
                       'positions_abs_intermediate', 'fix_vol', 'species_gamma', 'species_q',
                       'supercell_kpoints_list', 'supercell_kpoints_mp_grid',
                       'supercell_kpoints_mp_spacing', 'supercell_kpoints_mp_offset',
                       'supercell_kpoint_list', 'supercell_kpoint_mp_grid',
                       'supercell_kpoint_mp_spacing', 'supercell_kpoint_mp_offset', 'supercell_matrix',
                       'nonlinear_constraints', 'external_efield', 'positions_noise', 'cell_noise',
                       'hubbard_u', 'hubbard_alpha', 'atomic_init', 'quantisation_axis',
                       'quantization_axis', 'jcoupling_site', 'chemical_potential',
                       'elnes_kpoint_list', 'elnes_kpoint_mp_grid', 'elnes_kpoint_mp_spacing',
                       'elnes_kpoint_mp_offset', 'snap_to_symmetry', 'spectral_kpoint_path',
                       'spectral_kpoint_path_spacing', 'spectral_kpoint_list',
                       'spectral_kpoint_mp_grid', 'spectral_kpoint_mp_spacing',
                       'spectral_kpoint_mp_offset', 'spectral_kpoints_path',
                       'spectral_kpoints_path_spacing', 'spectral_kpoints_list',
                       'spectral_kpoints_mp_grid', 'spectral_kpoints_mp_spacing',
                       'spectral_kpoints_mp_offset']

valid_parameters_keywords = ['comment', 'iprint', 'continuation', 'reuse', 'checkpoint', 'task',
                             'calculate_stress', 'calculate_densdiff', 'calculate_elf',
                             'calculate_hirshfeld', 'run_time', 'backup_interval', 'num_backup_iter',
                             'print_clock', 'print_memory_usage', 'write_formatted_potential',
                             'write_formatted_density', 'write_formatted_elf', 'write_orbitals',
                             'write_cif_structure', 'write_cell_structure', 'write_bib',
                             'calc_molecular_dipole', 'write_checkpoint', 'cml_output', 'cml_filename',
                             'length_unit', 'mass_unit', 'time_unit', 'charge_unit', 'spin_unit',
                             'energy_unit', 'force_unit', 'velocity_unit', 'pressure_unit',
                             'inv_length_unit', 'frequency_unit', 'force_constant_unit', 'volume_unit',
                             'ir_intensity_unit', 'dipole_unit', 'efield_unit', 'entropy_unit', 'page_wvfns',
                             'rand_seed', 'data_distribution', 'opt_strategy', 'opt_strategy_bias',
                             'num_farms', 'num_proc_in_smp', 'num_proc_in_smp_fine', 'message_size', 'stop',
                             'xc_functional', 'xc_vxc_deriv_epsilon', 'relativistic_treatment', 'sedc_apply',
                             'sedc_scheme', 'sedc_sr_ts', 'sedc_d_ts', 'sedc_s6_g06', 'sedc_d_g06',
                             'sedc_lambda_obs', 'sedc_n_obs', 'sedc_sr_jchs', 'sedc_s6_jchs', 'sedc_d_jchs',
                             'page_ex_pot', 'nlxc_page_ex_pot', 'ppd_integral', 'nlxc_ppd_integral',
                             'ppd_size_x', 'nlxc_ppd_size_x', 'ppd_size_y', 'nlxc_ppd_size_y', 'ppd_size_z',
                             'nlxc_ppd_size_z', 'impose_trs', 'nlxc_impose_trs', 'exchange_reflect_kpts',
                             'nlxc_exchange_reflect_kpts', 'k_scrn_den_function', 'nlxc_k_scrn_den_function',
                             'k_scrn_averaging_scheme', 'nlxc_k_scrn_averaging_scheme', 're_est_k_scrn',
                             'nlxc_re_est_k_scrn', 'nlxc_exchange_screening', 'nlxc_exchange_fraction',
                             'calc_full_ex_pot', 'nlxc_calc_full_ex_pot', 'nlxc_div_corr_on',
                             'nlxc_div_corr_s_width', 'nlxc_div_corr_tol', 'nlxc_div_corr_npts_step',
                             'pspot_nonlocal_type', 'pspot_beta_phi_type', 'basis_precision',
                             'cut_off_energy', 'grid_scale', 'fine_grid_scale', 'fine_gmax',
                             'fft_max_prime_factor', 'fixed_npw', 'finite_basis_corr', 'basis_de_dloge',
                             'finite_basis_npoints', 'finite_basis_spacing', 'nelectrons', 'charge', 'spin',
                             'nup', 'ndown', 'spin_polarized', 'spin_polarised', 'nspins', 'nextra_bands',
                             'perc_extra_bands', 'nbands', 'elec_temp', 'electronic_minimizer',
                             'max_sd_steps', 'max_cg_steps', 'max_diis_steps', 'elec_method',
                             'metals_method', 'elec_energy_tol', 'elec_eigenvalue_tol', 'elec_force_tol',
                             'elec_convergence_win', 'max_scf_cycles', 'spin_fix', 'fix_occupancy',
                             'smearing_scheme', 'smearing_width', 'efermi_tol', 'num_occ_cycles',
                             'dipole_correction', 'dipole_dir', 'elec_dump_file', 'num_dump_cycles',
                             'elec_restore_file', 'mixing_scheme', 'mix_history_length', 'mix_charge_amp',
                             'mix_charge_gmax', 'mix_spin_amp', 'mix_spin_gmax', 'mix_cut_off_energy',
                             'mix_metric_q', 'popn_calculate', 'popn_bond_cutoff', 'pdos_calculate_weights',
                             'bs_max_iter', 'bs_max_cg_steps', 'bs_nextra_bands', 'bs_perc_extra_bands',
                             'bs_nbands', 'bs_eigenvalue_tol', 'bs_xc_functional', 'bs_re_est_k_scrn',
                             'bs_write_eigenvalues', 'geom_method', 'geom_max_iter', 'geom_energy_tol',
                             'geom_force_tol', 'geom_disp_tol', 'geom_stress_tol', 'geom_convergence_win',
                             'geom_modulus_est', 'geom_frequency_est', 'geom_spin_fix', 'geom_use_linmin',
                             'geom_linmin_tol', 'geom_lbfgs_max_updates', 'geom_tpsd_init_stepsize',
                             'geom_tpsd_iterchange', 'md_num_iter', 'md_delta_t', 'md_ensemble',
                             'md_use_pathint', 'md_num_beads', 'md_pathint_init', 'md_pathint_staging',
                             'md_pathint_num_stages', 'md_temperature', 'md_thermostat', 'md_barostat',
                             'md_ion_t', 'md_cell_t', 'md_nhc_length', 'md_nose_t', 'md_langevin_t',
                             'md_extrap', 'md_extrap_fit', 'md_xlbomd', 'md_xlbomd_history',
                             'md_damping_scheme', 'md_damping_reset', 'md_opt_damped_delta_t',
                             'md_elec_energy_tol', 'md_elec_eigenvalue_tol', 'md_elec_force_tol',
                             'md_elec_convergence_win', 'md_sample_iter', 'md_eqm_method', 'md_eqm_ion_t',
                             'md_eqm_cell_t', 'md_eqm_t', 'md_use_plumed', 'optics_nextra_bands',
                             'optics_perc_extra_bands', 'optics_nbands', 'optics_xc_functional',
                             'tssearch_method', 'tssearch_lstqst_protocol', 'tssearch_qst_max_iter',
                             'tssearch_cg_max_iter', 'tssearch_max_path_points', 'tssearch_force_tol',
                             'tssearch_disp_tol', 'tssearch_energy_tol', 'phonon_const_basis',
                             'phonon_energy_tol', 'phonon_max_cg_steps', 'phonon_max_cycles',
                             'phonon_convergence_win', 'phonon_preconditioner', 'phonon_use_kpoint_symmetry',
                             'phonon_calculate_dos', 'phonon_dos_spacing', 'phonon_dos_limit',
                             'phonon_finite_disp', 'phonon_fine_cutoff_method',
                             'phonon_force_constant_cutoff', 'phonon_force_constant_cut_scale',
                             'phonon_force_constant_ellipsoid', 'phonon_fine_method', 'phonon_method',
                             'phonon_dfpt_method', 'secondd_method', 'phonon_calc_lo_to_splitting',
                             'phonon_sum_rule', 'phonon_sum_rule_method', 'calculate_born_charges',
                             'born_charge_sum_rule', 'calculate_raman', 'raman_range_low',
                             'raman_range_high', 'phonon_write_force_constants', 'phonon_write_dynamical',
                             'efield_dfpt_method', 'efield_max_cg_steps', 'efield_max_cycles',
                             'efield_convergence_win', 'efield_energy_tol', 'efield_calc_ion_permittivity',
                             'efield_ignore_molec_modes', 'efield_freq_spacing', 'efield_oscillator_q',
                             'thermo_calculate_helmholtz', 'thermo_t_start', 'thermo_t_stop',
                             'thermo_t_spacing', 'thermo_t_npoints', 'wannier_spread_tol',
                             'wannier_max_sd_steps', 'wannier_sd_step', 'wannier_print_cube',
                             'wannier_spread_type', 'wannier_min_algor', 'wannier_ion_moments',
                             'wannier_ion_rmax', 'wannier_ion_cut', 'wannier_ion_cut_fraction',
                             'wannier_restart', 'wannier_ion_cut_tol', 'wannier_ion_cmoments', 'magres_task',
                             'magres_method', 'magres_max_cg_steps', 'magres_convergence_win',
                             'magres_conv_tol', 'magres_xc_functional', 'magres_max_sc_cycles',
                             'magres_jcoupling_task', 'magres_write_response', 'elnes_nextra_bands',
                             'elnes_perc_extra_bands', 'elnes_nbands', 'elnes_xc_functional',
                             'elnes_eigenvalue_tol', 'spectral_theory', 'spectral_task', 'spectral_max_iter',
                             'spectral_max_steps_per_iter', 'spectral_nextra_bands',
                             'spectral_perc_extra_bands', 'spectral_nbands', 'spectral_eigenvalue_tol',
                             'spectral_xc_functional', 'spectral_re_est_k_scrn',
                             'spectral_write_eigenvalues', 'tddft_num_states', 'tddft_selected_state',
                             'tddft_eigenvalue_tol', 'tddft_convergence_win', 'tddft_max_iter',
                             'tddft_nextra_states', 'tddft_xc_functional', 'tddft_method',
                             'tddft_eigenvalue_method', 'tddft_approximation', 'ga_pop_size', 'ga_max_gens',
                             'ga_mutate_rate', 'ga_mutate_amp', 'ga_fixed_n', 'ga_bulk_slice',
                             'excited_state_scissors', 'devel_code']

class CastepCell(OrderedDict):
    """Class to wrap a CASTEP cell (.cell) file"""

    def __init__(self, cellfile=None,xml=None,atoms=None):
        OrderedDict.__init__(self)
        if cellfile is not None:
            self.read(cellfile)
        elif xml is not None:
            self.read_xml(xml)
        elif atoms is not None:
            self.update_from_atoms(atoms)

    def copy(self):
        new = CastepCell()
        new.update(self)
        return new

    def read(self, cellfile):
        "Read a CASTEP .cell file. cellfile can be a mapping type, filename or an open file"

        if operator.isMappingType(cellfile):
            self.update(cellfile)
            return

        if type(cellfile) == type(''):
            cellfile = open(cellfile,'r')

        current_block = None

        for line in cellfile:
            line = line.strip()

            # Skip comments and blank lines
            if line.startswith('#') or line.startswith('!') or line == '':
                continue

            # Start of block
            if line.upper().startswith('%BLOCK'):
                block, name = map(string.strip, line.split())
                if current_block is None:
                    current_block = name.upper()
                    self[current_block] = []
                else:
                    raise ValueError("Parse error in cell file: can't begin new block %s when already in block %s" \
                        % (name, current_block))

            # End of block
            elif line.upper().startswith('%ENDBLOCK'):
                block, name = map(string.strip, line.split())
                name = name.upper()
                if name == current_block:
                    current_block = None
                else:
                    raise ValueError('Parse error in cell file: endblock %s does not match %s' % (name, current_block))

            # Stand-alone line
            else:
                if current_block is not None:
                    self[current_block].append(line)
                else:
                    try:
                        for c in ':=':  # Remove delimeters
                            line = line.replace(c,'',1)

                        fields = line.split()
                        key = fields[0].lower()
                        if not key in valid_cell_keywords:
                            raise ValueError('Unknown CELL keyword %s'% key)
                        self[key] = ' '.join(fields[1:])
                    except ValueError:
                        raise ValueError('Error parsing cell file line: %s' % line)


    def read_xml(self, xml):
        els = xml.documentElement.getElementsByTagName('castep_cell')

        if els == []:
            raise ValueError('No <castep_cell> element in XML file')

        if len(els) != 1:
            raise ValueError('More than one <castep_cell> element in XML file')

        if len(els[0].childNodes) != 1 or els[0].firstChild.nodeType != xml.TEXT_NODE:
            raise ValueError('<castep_cell> element should have exactly one Text node')

        self.read(els[0].firstChild.data.split('\n'))


    def write(self, cellfile=sys.stdout):
        "Write CASTEP .cell file. cellfile can be a filename or an open file"

        if type(cellfile) == type(''):
            cellfile = open(cellfile,'w')

        for key, value in self.iteritems():
            if type(value) == type([]):
                cellfile.write('\n%BLOCK '+key+'\n')
                for line in value:
                    cellfile.write('  '+line+'\n')
                cellfile.write('%ENDBLOCK '+key+'\n\n')
            else:
                value = str(value)
                value = value.replace('[','')
                value = value.replace(']','')
                cellfile.write('%s  %s\n' % (key, value))


    def to_atoms(self):
        if not self.has_key('LATTICE_CART') and not self.has_key('LATTICE_ABC'):
            raise ValueError('cell is missing LATTICE_CART and LATTICE_ABC blocks')

        if self.has_key('LATTICE_CART') and self.has_key('LATTICE_ABC'):
            raise ValueError('cell has both LATTICE_CART and LATTICE_ABC')

        if not self.has_key('POSITIONS_ABS') and not self.has_key('POSITIONS_FRAC'):
            raise ValueError('cell is missing POSITIONS_ABS and POSITIONS_FRAC blocks')

        if self.has_key('POSITIONS_ABS') and self.has_key('POSITIONS_FRAC'):
            raise ValueError('cell has both POSITIONS_ABS and POSITIONS_FRAC')

        if self.has_key('POSITIONS_ABS'): pos_block = self['POSITIONS_ABS']
        else:                             pos_block = self['POSITIONS_FRAC']

        if self.has_key('LATTICE_CART'):
            block = self['LATTICE_CART']
            if len(block) == 4:
                unit = block[0].lower()
                if unit not in castep_units:
                    raise ValueError("Don't know how to convert from units of %s" % block[0])
                block = block[1:]
            else:
                unit = 'ang'

            lattice = farray([ [float(x) for x in row] for row in map(string.split, block) ]).T*castep_units[unit]
        else:
            block = self['LATTICE_ABC']
            if len(block) == 3:
                unit = block[0].lower()
                if unit not in castep_units:
                    raise ValueError("Don't know how to convert from units of %s" % block[0])
                block = block[1:]
            else:
                unit = 'ang'

            a, b, c = [float(x)*castep_units[unit] for x in block[0].split()]
            alpha, beta, gamma = [float(x)*pi/180.0 for x in block[1].split()]

            lattice = make_lattice(a,b,c,alpha,beta,gamma)

        # Check if first line in position block is a valid unit
        unit = pos_block[0].lower()
        if unit in castep_units:
            pos_block = pos_block[1:]
        else:
            unit = 'ang'

        atoms = Atoms(n=len(pos_block),lattice=lattice)

        field_list = [line.split() for line in pos_block]
        elements = map(operator.itemgetter(0), field_list)

        # Look up names of elements specified by atomic number
        elements = [ not el.isdigit() and atomic_number(el) or el for el in elements ]

        # Set the element and pos data
        atoms.set_atoms(elements)
        atoms.pos[:,:] = farray([ [float(x)*castep_units[unit] for x in row] \
                                  for row in [field[1:4] for field in field_list]]).T

        # Convert from fractional to real positions
        if self.has_key('POSITIONS_FRAC'):
            atoms.pos[:,:] = np.dot(atoms.lattice, atoms.pos)

        return atoms

    def update_from_atoms(self, at, frac_pos=False, lattice_abc=False):

        # Add lattice to cell
        if lattice_abc:
            self['LATTICE_ABC'] = []
            a, b, c, alpha, beta, gamma = get_lattice_params(at.lattice)
            self['LATTICE_ABC'].append('%f %f %f' % (a,b,c))
            self['LATTICE_ABC'].append('%f %f %f' % (alpha*180.0/pi,beta*180.0/pi,gamma*180.0/pi))

        else:
            self['LATTICE_CART'] = []
            for i in frange(3):
                self['LATTICE_CART'].append('%f %f %f' % tuple(at.lattice[:,i]))

        # Add atomic positions to cell
        if frac_pos:
            self['POSITIONS_FRAC'] = []
            for i in frange(at.n):
                self['POSITIONS_FRAC'].append(at.species[:,i].stripstrings() +' %f %f %f' % tuple(np.dot(at.g,at.pos[:,i])))
        else:
            self['POSITIONS_ABS'] = []
            for i in frange(at.n):
                self['POSITIONS_ABS'].append(at.species[:,i].stripstrings() +' %f %f %f' % tuple(at.pos[:,i]))

        for p in at.params:
            if p.lower() in valid_cell_keywords:
                self[p] = at.params[p]


    @staticmethod
    @atoms_reader('cell')
    def cellreader(source, format=None):
        self = CastepCell(source)
        yield self.to_atoms()

class CastepCellWriter(object):
    def __init__(self, dest):
        if dest == 'stdout': dest = sys.stdout
        self.dest = dest
        self.n = 0

    def write(self, at, cell_template=None, param_template=None):
        dest = self.dest
        if '%' in dest: dest = self.dest % self.n
        if cell_template is None:
            cell = CastepCell()
        else:
            cell = CastepCell(cell_template)
        cell.update_from_atoms(at)
        cell.write(dest)
        if param_template is not None:
            param = CastepParam(param_template)
            param.write(os.path.splitext(dest)[0]+'.param')
        self.n += 1

AtomsWriters['cell'] = CastepCellWriter

class CastepParam(OrderedDict):
    "Class to wrap a CASTEP parameter (.param) file"

    def __init__(self, paramfile=None, xml=None, atoms=None):
        OrderedDict.__init__(self)
        if paramfile is not None:
            self.read(paramfile)
        elif xml is not None:
            self.read_xml(xml)
        elif atoms is not None:
            self.update_from_atoms(atoms)

    def copy(self):
        new = CastepParam()
        new.update(self)
        return new

    def read(self, paramfile):
        "Read a CASTEP .param file. paramfile can be a filename or an open file"

        if operator.isMappingType(paramfile):
            self.update(paramfile)
            return

        if type(paramfile) == type(''):
            paramfile = open(paramfile,'r')

        for line in paramfile:
            line = line.strip()

            # Skip comments and blank lines
            if line.startswith('#') or line.startswith('!') or line == '':
                continue

            line = re.compile('[:=]').sub(' ', line, 1)
            fields = line.split()
            key = fields[0].lower()
            if not key in valid_parameters_keywords:
                raise ValueError('Unknown PARAMETERS keyword %s' % key)
            value = ' '.join(fields[1:])
            value = castep_value_map.get(value, value)
            self[key] = value

    def read_xml(self, xml):
        els = xml.documentElement.getElementsByTagName('castep_param')

        if els == []:
            raise ValueError('No <castep_param> element in XML file')

        if len(els) != 1:
            raise ValueError('More than one <castep_param> element in XML file')

        if len(els[0].childNodes) != 1 or els[0].firstChild.nodeType != xml.TEXT_NODE:
            raise ValueError('<castep_param> element should have exactly one Text node')

        self.read(els[0].firstChild.data.split('\n'))

    def read_from_castep_output(self, castep_output):
        "Read user parameters from .castep output. Input should be filename, file-like object or list of lines"

        if type(castep_output) == type(''):
            f = open(castep_output, 'r')
            castep_output = f.readlines()
            f.close()
        elif hasattr(castep_output, 'read'):
            castep_output = castep_output.readlines()

        # Remove newlines from end of each line in castep_output
        castep_output = [line.rstrip() for line in castep_output]

        # Bizarrelly, the parameter names in the output file have the underscores
        # removed from them. Here we construct a mapping from the .castep names to
        # the true .param file
        param_lookup = dict(zip(map(lambda s: s.replace('_',''),valid_parameters_keywords),
             valid_parameters_keywords))

        # Find the user parameters section
        try:
            user_param_start = castep_output.index(' ******************************* User Parameters *******************************')
        except ValueError:
            raise ValueError('No user parameters found in castep output')

        param_lines = []
        i = user_param_start + 1
        while castep_output[i].strip():
            line = castep_output[i]
            key, value = map(string.strip, line[:line.index('#')].split(':',1))
            if not key in param_lookup:
                raise ValueError('Unknown parameter %s in castep output file' % key)
            param_lines.append('%s = %s' % (param_lookup[key], value))
            i += 1

        self.read(param_lines)

    def update_from_atoms(self, at):
        for p in at.params:
            if p.lower() in valid_parameters_keywords:
                self[p] = at.params[p]

    def write(self,paramfile=sys.stdout):
        "Write CASTEP .param file"

        if type(paramfile) == type(''):
            paramfile = open(paramfile,'w')

        for key, value in self.iteritems():
            value = castep_output_map.get(value, value)
            paramfile.write('%s = %s\n' % (key, value))

@atoms_reader('geom')
@atoms_reader('md')
def CastepGeomMDReader(source, atoms_ref=None, format=None):
    """Generator to read frames from CASTEP .geom and .md files"""

    if type(source) == type(''):
        source = open(source, 'r')
    source = iter(source)

    if atoms_ref is not None and not atoms_ref.has_property('frac_pos'):
        atoms_ref.add_property('frac_pos', 0.0, n_cols=3)
        atoms_ref.frac_pos[:] = np.dot(atoms_ref.g, atoms_ref.pos)

    while True:
        lines = []
        line = source.next().strip() # raises StopIteration at end of file

        # Skip header if present
        if line.startswith('BEGIN header'):
            while not line.startswith('END header'):
                line = source.next().strip()
            source.next() # skip blank line
        else:
            lines.append(line)

        # Read until next blank line
        while line != '':
            line = source.next().strip()
            if line != '': lines.append(line)

        if len(lines) <= 1:
            raise StopIteration

        params = Dictionary()

        # First line is the time/step in a.u - we convert to fs
        params['time'] = float(lines[0])*AU_FS

        # Then the energy, in Hartree - we convert to eV
        energy_lines = filter(lambda s: s.endswith('<-- E'), lines)
        if len(energy_lines) != 1:
            raise ValueError('Number of energy lines should be exactly one. Got %r' % energy_lines)

        params['energy'], params['hamiltonian'] = \
                               [float(x)*HARTREE for x in energy_lines[0].split()[0:2]]

        # Temperature, in atomic units - we convert to Kelvin
        temperature_lines = filter(lambda s: s.endswith('<-- T'), lines)
        if temperature_lines:
            if len(temperature_lines) != 1:
                raise ValueError('Number of temperature lines should be exactly one. Got %r' % temperature_lines)
            params['temperature'] = float(temperature_lines[0].split()[0])*HARTREE/BOLTZMANN_K

        # Pressure, in atomic units - we convert to ev/A**3
        pressure_lines = filter(lambda s: s.endswith('<-- P'), lines)
        if pressure_lines:
            if len(pressure_lines) != 1:
                raise ValueError('Number of pressure lines should be exactly one. Got %r' % pressure_lines)
            params['pressure'] = float(pressure_lines[0].split()[0])*HARTREE/BOHR**3

        # Lattice is next, in units of Bohr
        lattice_lines = filter(lambda s: s.endswith('<-- h'), lines)
        if lattice_lines:
            lattice = farray([ [float(x)* BOHR for x in row[0:3]]
                               for row in map(string.split, lattice_lines) ]).T
        else:
            if atoms_ref is None:
                raise ValueError('No lattice in .geom file and atoms_ref not present')
            else:
                lattice = atoms_ref.lattice[:]

        # Then optionally virial tensor - convert stress tensor to eV/A**3
        # then multiply by cell volume to get virial in eV
        stress_lines  = filter(lambda s: s.endswith('<-- S'), lines)
        if stress_lines:
            virial = farray([ [float(x)*(HARTREE/(BOHR**3)) for x in row[0:3]]
                              for row in map(string.split, stress_lines) ]).T

        # Find positions and forces
        poslines   = filter(lambda s: s.endswith('<-- R'), lines)
        velolines  = filter(lambda s: s.endswith('<-- V'), lines)
        forcelines = filter(lambda s: s.endswith('<-- F'), lines)

        if poslines and forcelines and len(poslines) != len(forcelines):
            raise ValueError('Number of pos lines (%d) != force lines (%d)'
                             % (len(poslines), len(forcelines)))

        n_atoms = max(len(poslines), len(forcelines))
        if poslines == [] and forcelines == []:
            if atoms_ref is None:
                raise ValueError('No positions or forces in .geom file and atoms_ref not present')
            n_atoms = atoms_ref.n

        if atoms_ref is not None:
            # If we were passed in an Atoms object, construct mapping from
            # CASTEP (element, number) to original atom index
            assert n_atoms == atoms_ref.n
            at = atoms_ref.copy()

            # remove things which shouldn't be inherited from atoms_ref
            for p in ['energy', 'hamiltonian', 'temperature', 'pressure', 'virial']:
                if p in at.params: del at.params[p]
            for p in ['force', 'velo']:
                if p in at.properties: del at.properties[p]

            at.set_lattice(lattice, scale_positions=False)
            at.params.update(params)
            species_count = {}
            lookup = {}
            for i in frange(at.n):
                el = at.species[i].stripstrings()
                if species_count.has_key(el):
                    species_count[el] += 1
                else:
                    species_count[el] = 1
                lookup[(el,species_count[el])] = i
        else:
            # Otherwise we make a new, empty Atoms object. Atoms will
            # be ordered as they are in .castep file.
            lookup = {}
            at = Atoms(n=n_atoms, lattice=lattice, params=params)

        # Now parse the positions, converting from units of Bohr
        if poslines:
            for i, line in fenumerate(poslines):
                el, num, x, y, z, arrow, label = line.split()
                num = int(num)
                if not (el,num) in lookup:
                    lookup[(el,num)] = i
                at.z[lookup[(el,num)]] = atomic_number(el)
                at.pos[:,lookup[(el,num)]] = [ float(r)* BOHR for r in (x,y,z)]

            at.set_atoms(at.z) # set at.species property to match at.z
        else:
            at.pos[:] = np.dot(at.lattice, at.frac_pos)

        # Velocities, if this is an MD file, from atomic units to A/fs
        if velolines:
            at.add_property('velo', 0.0, n_cols=3)
            for i, line in fenumerate(velolines):
                el, num, vx, vy, vz, arrow, label = line.split()
                num = int(num)
                at.velo[:,lookup[(el,num)]] = [ float(f)*BOHR/AU_FS for f in (vx, vy, vz) ]

        # And finally the forces, which are in units of Hartree/Bohr
        if forcelines:
            at.add_property('force', 0.0, n_cols=3)
            for i, line in fenumerate(forcelines):
                el, num, fx, fy, fz, arrow, label = line.split()
                num = int(num)
                at.force[:,lookup[(el,num)]] = [ float(f)*HARTREE/BOHR for f in (fx, fy, fz) ]

        if stress_lines:
            at.params['virial'] = -virial*at.cell_volume()

        if atoms_ref is None:
            atoms_ref = at.copy()

        yield at

# Synonyms
CastepGeomReader = CastepMDReader = CastepGeomMDReader

@atoms_reader('castep')
@atoms_reader('castep_log')
def CastepOutputReader(castep_file, atoms_ref=None, abort=False, format=None):
    """Parse .castep file, and return Atoms object with positions,
       energy, forces, and possibly stress and atomic populations as
       well"""

    if type(castep_file) == type(''):
        castep_file_name = castep_file
        castep_file = open(castep_file,'r')
    else:
        castep_file_name = '<open file>'
    castep_file = iter(castep_file)

    param = CastepParam()

    if atoms_ref is not None and not atoms_ref.has_property('frac_pos'):
        atoms_ref.add_property('frac_pos', 0.0, n_cols=3)
        atoms_ref.frac_pos[:] = np.dot(atoms_ref.g, atoms_ref.pos)

    got_header = False
    eof = False
    while True:

        castep_output = []
        while True:
            try:
                line = castep_file.next()
            except StopIteration:
                eof = True
                break
            castep_output.append(line)

            if line.startswith(' Atomic calculation performed'):
                if got_header:
                    got_param = True
                else:
                    got_header = True

            if line.startswith(' Starting MD iteration'):
                break

            if line.startswith(' Starting BFGS iteration'):
                break

            if line.startswith(' BFGS: improving iteration'):
                break

        if castep_output == []:
            break

        # Remove newlines from end of each line in castep_output
        castep_output = [line.rstrip() for line in castep_output]

        # NB: CASTEP doesn't always print 'Total time'
        run_time = None
        total_time = filter(lambda s: s.startswith('Total time'), castep_output)
        if total_time == []:
            has_converged = filter(lambda s: s.startswith('Total energy has converged'), castep_output)
            if abort and has_converged == []:
                raise ValueError("castep didn't complete")
        else:
            run_time = float(total_time[0].split()[3])

        # Now we should have contents of a valid .castep file in castep_output

        # First let's read the user parameters for this run from top of file
        new_param = CastepParam()
        try:
            new_param.read_from_castep_output(castep_output)
            param.update(new_param)
        except ValueError:
            # If iprint < 2, params are not in .castep, file so let's look for .param file in same place
            param_file = os.path.splitext(castep_file_name)[0]+'.param'
            if os.path.exists(param_file):
                new_param.read(param_file)
                param.update(new_param)
            elif abort:
                raise

        # Next let's extract the lattice and atomic positions
        lattice_lines = [i for (i,x) in enumerate(castep_output) if x == '                                      Unit Cell']

        if lattice_lines == []:
            if atoms_ref is None or abort:
                raise ValueError('No unit cell found in castep file - try passing atoms_ref')
            else:
                lattice = atoms_ref.lattice.copy()
        else:
            lattice_line = lattice_lines[-1] # last lattice

            lattice_lines = castep_output[lattice_line+3:lattice_line+6]
            lattice = fzeros((3,3))
            lattice[:,1] = map(float, lattice_lines[0].split()[0:3])
            lattice[:,2] = map(float, lattice_lines[1].split()[0:3])
            lattice[:,3] = map(float, lattice_lines[2].split()[0:3])

        cell_contents = [i for (i,x) in  enumerate(castep_output) if x == '                                     Cell Contents']
        if cell_contents == []:
            if atoms_ref is None or abort:
                raise ValueError('No cell contents found in castep file - try passing atoms_ref')
            else:
                logging.warning('No cell contents. If this is a variable cell geometry optimisation with fixed ions this is normal')
                n_atoms = atoms_ref.n

        if cell_contents:
            cell_first_line = cell_contents[-1] # last cell contents line

            try:
                n_atoms = int(castep_output[cell_first_line+2].split()[-1])
                offset = 10
            except IndexError:
                offset = 7
                n_atoms = 0
                while 'xxxxxxxxxxx' not in castep_output[cell_first_line+offset+n_atoms]:
                   n_atoms += 1

        if atoms_ref is not None:
            # If we were passed in an Atoms object, construct mapping from
            # CASTEP (element, number) to original atom index
            assert n_atoms == atoms_ref.n, 'Number of atoms must match atoms_ref'
            atoms = atoms_ref.copy()

            # remove things which shouldn't be inherited from atoms_ref
            for p in ['energy', 'enthalpy', 'virial']:
                if p in atoms.params: del atoms.params[p]
            for p in ['force', 'popn_s', 'popn_p', 'popn_d', 'popn_f', 'popn_total', 'popn_charge', 'popn_spin']:
                if p in atoms.properties: del atoms.properties[p]

            atoms.set_lattice(lattice, scale_positions=False)
            species_count = {}
            lookup = {}
            for i in frange(atoms.n):
                el = atoms.species[i].stripstrings()
                if species_count.has_key(el):
                    species_count[el] += 1
                else:
                    species_count[el] = 1
                lookup[(el,species_count[el])] = i
        else:
            # Otherwise we make a new, empty Atoms object. Atoms will
            # be ordered as they are in .castep file.
            lookup = {}
            atoms = Atoms(n=n_atoms,lattice=lattice)

        if cell_contents:
            cell_lines = castep_output[cell_first_line+offset:cell_first_line+offset+n_atoms]

            # Fill in species and fractional positions
            atoms.add_property('frac_pos',0.0,n_cols=3)
            for i, line in fenumerate(cell_lines):
                x1, el, num, u, v, w, x2 = line.split()
                num = int(num)
                if not (el,num) in lookup:
                    lookup[(el,num)] = i
                atoms.z[lookup[(el,num)]] = atomic_number(el)
                atoms.frac_pos[:,lookup[(el,num)]] = map(float, (u,v,w))

            atoms.set_atoms(atoms.z) # set at.species from at.z

        # Calculate cartesian postions from fractional positions
        # (correct if we're doing a variable cell geom. opt. with fixed ions)
        atoms.pos[:] = np.dot(atoms.lattice, atoms.frac_pos)

        energy_lines = filter(lambda s: s.startswith('Final energy') and not s.endswith('<- EDFT'), castep_output)

        # If we're using smearing, correct energy is 'free energy'
        energy_lines.extend(filter(lambda s: s.startswith('Final free energy (E-TS)'), castep_output))

        # Are we doing finite basis correction?
        energy_lines.extend(filter(lambda s: s.startswith(' Total energy corrected for finite basis set'), castep_output))

        if (len(energy_lines) == 0):
            if abort:
                raise ValueError('No total energy found in castep file')
        else:
            # Use last matching energy line in file
            # Has this energy been corrected for the finite basis set?
            energy_param_name = 'energy'
            if 'not corrected for finite basis set' in energy_lines[-1]: energy_param_name = 'energy_no_basis_corr'

            fields = energy_lines[-1].split()
            if 'eV' not in fields:
                if abort:
                    raise ValueError('No value found in energy line "%s"' % energy_lines[-1])
            else:
                atoms.params[energy_param_name] = float(fields[fields.index('eV')-1])

        # If we're doing geom-opt, look for enthalpy
        enthalpy_lines = [s for s in castep_output if s.startswith(' BFGS: finished iteration ') or
                          s.startswith(' BFGS: Final Enthalpy')]
        if enthalpy_lines:
            atoms.params['enthalpy'] = float(enthalpy_lines[-1].split()[-2])

        try:

            for fn in ('Forces', 'Symmetrised Forces'):
                force_start_lines = [i for (i,s) in enumerate(castep_output) if s.find('****** %s ******' % fn) != -1]
                if force_start_lines: break

            if force_start_lines == []:
                raise ValueError

            # Use last set of forces
            force_start = force_start_lines[-1]

            # Extract force lines from .castep file
            force_lines = castep_output[force_start+6:force_start+6+atoms.n]

            # remove "cons" tags from constrained degrees of freedom
            force_lines = [ s.replace("(cons' d)", "") for s in force_lines ]

            atoms.add_property('force',0.0,n_cols=3)

            # Fill in the forces
            for i, line in enumerate(force_lines):
                line = line.replace('*','') # Remove the *s
                el, num_str, fx, fy, fz = line.split()
                num = int(num_str)
                atoms.force[:,lookup[(el,num)]] = (fx,fy,fz)

        except ValueError, m:
            if abort:
                raise ValueError('No forces found in castep file %s: ' % m)

        # Individual contributions to total force
        for name,label in [('force_ewald', 'Ewald forces'),
                           ('force_locpot', 'Local potential forces'),
                           ('force_nlpot', 'Non-local potential forces'),
                           ('force_extpot', 'External potential forces')]:
            force_start_lines = [i for (i,s) in enumerate(castep_output) if s.find('****** %s ******' % label) != -1]
            if force_start_lines == []: continue

            # Use last set of forces
            force_start = force_start_lines[-1]

            # Extract force lines from .castep file
            force_lines = castep_output[force_start+6:force_start+6+atoms.n]

            # remove "cons" tags from constrained degrees of freedom
            force_lines = [ s.replace("(cons' d)", "") for s in force_lines ]

            atoms.add_property(name,0.0,overwrite=True,n_cols=3)

            # Fill in the forces
            for i, line in enumerate(force_lines):
                line = line.replace('*','') # Remove the *s
                el, num_str, fx, fy, fz = line.split()
                num = int(num_str)
                getattr(atoms, name)[:,lookup[(el,num)]] = (fx,fy,fz)

        # Have we calculated stress?
        got_virial = False
        for sn in ('Stress Tensor', 'Symmetrised Stress Tensor'):
            stress_start_lines = [i for i,s in enumerate(castep_output) if s.find('****** %s ******' % sn) != -1 ]
            if stress_start_lines: break

        if stress_start_lines:
            stress_start = stress_start_lines[-1]
            stress_lines = castep_output[stress_start+6:stress_start+9]
            virial = fzeros((3,3),float)
            for i, line in fenumerate(stress_lines):
                star1, label, vx, vy, vz, star2 = line.split()
                virial[:,i] = [-float(v) for v in (vx,vy,vz) ]
            got_virial = True

        spin_polarised = 'spin_polarised' in param and param['spin_polarised']

        # Have we calculated local populations and charges?
        if 'popn_calculate' in param and param['popn_calculate']:
            try:
                try:
                    popn_start = castep_output.index('     Atomic Populations')
                except ValueError:
                    try:
                        popn_start = castep_output.index('     Atomic Populations (Mulliken)')
                    except ValueError:
                        raise

                popn_lines = castep_output[popn_start+4:popn_start+4+atoms.n]

                atoms.add_property('popn_s',0.0)
                atoms.add_property('popn_p',0.0)
                atoms.add_property('popn_d',0.0)
                atoms.add_property('popn_f',0.0)
                atoms.add_property('popn_total',0.0)
                atoms.add_property('popn_charge',0.0)
                if spin_polarised:
                    atoms.add_property('popn_spin', 0.0)

                for line in popn_lines:
                    if spin_polarised:
                        el, num_str, s, p, d, f, tot, charge, spin = line.split()
                    else:
                        el, num_str, s, p, d, f, tot, charge = line.split()
                    num = int(num_str)
                    atoms.popn_s[lookup[(el,num)]] = float(s)
                    atoms.popn_p[lookup[(el,num)]] = float(p)
                    atoms.popn_d[lookup[(el,num)]] = float(d)
                    atoms.popn_f[lookup[(el,num)]] = float(f)
                    atoms.popn_total[lookup[(el,num)]] = float(tot)
                    atoms.popn_charge[lookup[(el,num)]] = float(charge)
                    if spin_polarised:
                        atoms.popn_spin[lookup[(el,num)]] = float(spin)

            except ValueError:
                if abort:
                    raise ValueError('No populations found in castep file')

        mod_param = param.copy()

        # check for calculation of dipole moment - one vector per spin component
        if ' calculate electric dipole moment of system     : on' in castep_output:
            dipole_magnitudes = [s for s in castep_output if s.startswith('  +  Magnitude of Dipole')]
            dipole_directions = [s for s in castep_output if s.startswith('  +  Direction of Dipole')]

            dipoles = []
            for mag_line, dir_line in zip(dipole_magnitudes, dipole_directions):
                magnitude = float(mag_line.split()[5])*DEBYE
                direction = farray([float(x) for x in dir_line.split()[5:8]])
                dipoles.append(magnitude*direction)

            for i,dipole in fenumerate(dipoles):
                mod_param['dipole%d' % i] = dipole


        # append K-point information
        try:
            kpoint_start = castep_output.index('                              k-Points For BZ Sampling')
            kp_mesh_line = castep_output[kpoint_start+2]
            fields = kp_mesh_line.split()
            mod_param['kpoints_mp_grid'] = [int(fields[-3]), int(fields[-2]), int(fields[-1])]
        except ValueError:
            pass

        mod_param['castep_file_name'] = castep_file_name

        if run_time is not None:
            mod_param['castep_run_time'] = run_time

        # Convert virial to libAtoms units and add to params
        if got_virial:
            mod_param['virial'] = virial*atoms.cell_volume()/GPA

        if atoms_ref is None:
            atoms_ref = atoms.copy()

        atoms.params.update(mod_param)

        yield atoms

        if eof:
            break

@atoms_reader('magres')
def MagresReader(source, atoms_ref=None, format=None):
    """Generator to read .magres files. The specification of the format my be found at http://www.ccpnc.ac.uk/"""

    if type(source) == type(''):
        source = open(source, 'r')
    source = iter(source)

    if atoms_ref is not None and not atoms_ref.has_property('frac_pos'):
        atoms_ref.add_property('frac_pos',0.0,n_cols=3)
        atoms_ref.frac_pos[:] = np.dot(atoms_ref.g,atoms_ref.pos)

    eof = False
    parsing_started = False
    while True:
       block = ""
       atom_line = []
       ms_line = []
       efg_line = []
       while True:
           try:
               line=source.next()
           except StopIteration:
               eof = True
               
           if re.match("#\$magres",line) or eof:
               if parsing_started:
                  break
               else:
                  parsing_started = True
               continue
           elif re.match("\[[^/].*\]",line):
               block = line.replace("[","").replace("]","")
               continue
           elif re.match("\[\/.*\]",line):
               block = ""
               continue

           if "atoms" in block:
               if line.startswith("lattice"):
                   lattice_line = line.split()[1:10]
               elif line.startswith("atom"):
                   atom_line.append(line.split())
           if "magres" in block:
               if line.startswith("ms"):
                   ms_line.append(line.split())
               if line.startswith("efg"):
                   efg_line.append(line.split())
       
       lattice = fzeros((3,3))
       lattice[:,:] = np.reshape(map(float,lattice_line),(3,3))

       n_atoms=len(atom_line)
       atoms = Atoms(n=n_atoms,lattice=lattice)
       lookup = {}

       for i, line in fenumerate(atom_line):
           label, atom_type, species_label, number_in_species_label, x, y, z = line
           if not (species_label, number_in_species_label) in lookup:
               lookup[(species_label, number_in_species_label)] = i
           atoms.z[lookup[(species_label, number_in_species_label)]] = atomic_number(atom_type)
           atoms.pos[:,lookup[(species_label, number_in_species_label)]] = map(float,(x,y,z))

       atoms.set_atoms(atoms.z) # set at.species from at.z

       if n_atoms == len(ms_line):
           atoms.add_property('ms',0.0,n_cols=9)
           for i, line in fenumerate(ms_line):
               (label, species_label, number_in_species_label,
                ms_xx, ms_xy, ms_xz,
                ms_yx, ms_yy, ms_yz,
                ms_zx, ms_zy, ms_zz) = line
               atoms.ms[:,lookup[(species_label, number_in_species_label)]] = map(float,(ms_xx, ms_xy, ms_xz,
                                                                                         ms_yx, ms_yy, ms_yz,
                                                                                         ms_zx, ms_zy, ms_zz))
       if n_atoms == len(efg_line):
           atoms.add_property('efg',0.0,n_cols=9)
           for i, line in fenumerate(efg_line):
               (label, species_label, number_in_species_label,
                efg_xx, efg_xy, efg_xz,
                efg_yx, efg_yy, efg_yz,
                efg_zx, efg_zy, efg_zz) = line
               atoms.efg[:,lookup[(species_label, number_in_species_label)]] = map(float,(efg_xx, efg_xy, efg_xz,
                                                                                         efg_yx, efg_yy, efg_yz,
                                                                                         efg_zx, efg_zy, efg_zz))

       yield atoms
       if eof:
          break

def get_valid_keywords(castep):
    """Determines valid cell and parameter keyword by invoking castep with -help parameter.
       Returns a tuple (valid_cell_keywords, valid_parameters_keywords)"""

    valid_cell_keywords = []
    valid_parameters_keywords = []

    if castep.find('%s') == -1:
        castep = castep + ' %s'

    lines = os.popen(castep % ('-help all')).readlines()
    try:
       cell_start = lines.index('Help information on CELL keywords:\n')
       param_start = lines.index('Help information on PARAMETERS keywords:\n')
    except ValueError:
       raise ValueError('Error parsing output of castep -help all')

    cell_lines = lines[cell_start+2:param_start-2]
    param_lines = lines[param_start+2:-1]

    for lines, keywords  in zip((cell_lines, param_lines),
          (valid_cell_keywords,valid_parameters_keywords)):
       lines = map(string.strip, lines)
       for line in lines:
          keywords.append(line.split()[0].lower())

    print 'valid_cell_keywords = %r\n' % valid_cell_keywords
    print 'valid_parameters_keywords = %r' % valid_parameters_keywords

    return (valid_cell_keywords, valid_parameters_keywords)


def check_pspots(cluster, cell, param, orig_dir):
    """Check pseudopotential files are present, and that we have one for
    each element present in cluster. Also tries to check spin polarisation
    of system matches that specified in parameters."""

    log = logging.getLogger('castep_driver')

    # Get list of unique elements in cluster
    pspot_dict = {}
    for el in cluster.species.stripstrings():
        pspot_dict[el] = None
    elements = pspot_dict.keys()
    log.info(str(elements))

    if 'SPECIES_POT' in cell:
        for line in cell['SPECIES_POT']:
            element, pspot = map(string.strip, line.split())
            pspot_dict[element] = pspot

    valence_charges = {}
    missing_pspots = False

    for el in elements:
        pspot = pspot_dict[el]
        if pspot is None:
            # Not fatal since on-the-fly pseudo-potential will be generated
            # however, we can't estimate then n_electrons here
            log.warn('Pseudopotential missing for element %s' % el)
            valence_charges[el] = '?'
            missing_pspots = True
        elif not pspot.startswith('/"'):
            shutil.copy(orig_dir+'/'+pspot, '.')
            pspot_lines = open(pspot,'r').readlines()

            # First check for old-style pspot report
            zlines = filter(lambda s: s.startswith('    |  z ='), pspot_lines)
            if len(zlines) == 1:
                # Extract number of valence electrons
                fields = zlines[0].split()
                zv = float(fields[fields.index('zv')+2])
            else:
                # Now try newer style for OTF generated pspots
                zlines = filter(lambda s: s.startswith('   | Element: '), pspot_lines)
                if len(zlines) != 1:
                    raise ValueError('Malformed pseudopotential file: does not contain pseudopotential report')
                fields = zlines[0].split()
                zv = float(fields[fields.index('charge:')+1])

            valence_charges[el] = zv
        else:
            # It's an on-the-fly pseudopotential
            fields = pspot[2:-2].split('|')
            valence_charges[el] = float(fields[0])

    log.info('unique elements and valence electrons: %s' % valence_charges)

    if not missing_pspots:
        n_electrons = reduce(operator.add, \
                             map(lambda el: valence_charges[el]*sum(cluster.species.stripstrings() == el), \
                                 elements))
        log.info('total electrons %.1f' % n_electrons)
        if (param.has_key('spin_polarised') and param['spin_polarised'] and \
            int(n_electrons) % 2 != 0):
            raise ValueError('spin polarised set to false but got odd number of electrons!')
        if (param.has_key('spin_polarised') and param['spin_polarised'] and \
            int(n_electrons) % 2 == 0):
            raise ValueError('spin polarised set to true but got even number of electrons!')


def run_castep(cell, param, stem, castep, castep_log=None, save_all_check_files=False, save_all_input_files=False, test_mode=False, copy_in_files=None, subdir=None):
    "Invoke castep and return True if it completed successfully"

    log = logging.getLogger('castep_driver')

    orig_dir = os.getcwd()
    try:
        if subdir is not None:
            if not os.path.exists(subdir):
                os.makedirs(subdir)
            os.chdir(subdir)

        # Write parameter file ...
        param.write(stem+'.param')
        # ... and cell file
        cell.write(stem+'.cell')

        if copy_in_files is not None:
            for pspot in copy_in_files:
                shutil.copyfile(pspot, os.path.join(os.getcwd(), os.path.basename(pspot)))

        if not '%s' in castep: castep = castep + ' %s'

        if test_mode:
            log.info('test mode: not running castep')

            ## Simulate making check file
            #check_file = open(stem+'.check','w')
            #check_file.close()

        else:
            # Remove old output file and error files
            try:
                os.remove(stem+'.castep')
            except:
                pass
            for f in glob.glob('%s*.err' % stem):
                os.remove(f)

            os.system(castep % stem)

        if save_all_check_files:
            if os.path.exists('%s.check' % stem):
                n = 0
                while os.path.exists('%s.check.%d' % (stem, n)):
                    n += 1
                shutil.copyfile('%s.check' % stem, '%s.check.%d' % (stem, n))

        if save_all_input_files:
            if os.path.exists('%s.cell' % stem):
                n = 0
                while os.path.exists('%s.cell.%d' % (stem, n)):
                    n += 1
                shutil.copyfile('%s.cell' % stem, '%s.cell.%d' % (stem, n))
            if os.path.exists('%s.param' % stem):
                n = 0
                while os.path.exists('%s.param.%d' % (stem, n)):
                    n += 1
                shutil.copyfile('%s.param' % stem, '%s.param.%d' % (stem, n))

        err_files = glob.glob('%s*.err' % stem)
        got_error = False
        if (len(err_files) > 0):
            for f in err_files:
                error_text = open(f).read().strip()
                if error_text != '':
                    got_error = True
                    log.error(error_text)

        # Write log file here so that it will always be written
        if castep_log is not None and os.path.exists('%s.castep' % stem):
            logf = open(castep_log, 'a')
            castep_output = open('%s.castep' % stem, 'r').readlines()
            logf.writelines(castep_output)
            logf.close()

    finally:
        os.chdir(orig_dir)

    return not got_error

from quippy import Potential

class CastepPotential(Potential):
    def __init__(self, cell=None, param=None, castep_exec='castep %s', stem='castep_callback', test_mode=False,
                 little_clusters=False, copy_in_files=None, subdir=None):

        Potential.__init__(self, 'CallbackPot little_clusters=%s' % ({True: 'T', False: 'F'}[little_clusters]))
        self.set_callback(self.run)

        if isinstance(cell, str):
            self.cell = CastepCell(cell)
        else:
            self.cell = cell

        if isinstance(param, str):
            self.param = CastepParam(param)
        else:
            self.param = param

        self.castep_exec = castep_exec
        self.stem = stem
        self.test_mode = test_mode
        self.copy_in_files = copy_in_files
        self.subdir = subdir

    def run(self, at):
        if self.cell is not None:
            cell = self.cell.copy()
        else:
            cell = CastepCell()
        cell.update_from_atoms(at)

        if self.param is not None:
            param = self.param.copy()
        else:
            param = CastepParam()
        param.update_from_atoms(at)

        stem = self.stem
        if os.path.exists(stem+'.castep'):
            n = 0
            while True:
                stem = '%s_%05d' % (self.stem, n)
                n += 1
                if not os.path.exists(stem+'.castep'): break

        run_castep(cell, param, stem, self.castep_exec, test_mode=self.test_mode, 
                   copy_in_files=self.copy_in_files, subdir=self.subdir)

        if not self.test_mode:
            print 'Reading from file %s' % (stem+'.castep')
            result = Atoms(stem+'.castep', atoms_ref=at)

            # Update params -- this includes contents of .cell and .param templates
            at.params.update(result.params)

            # Energy and force
            at.params['energy'] = float(result.energy)
            at.add_property('force', result.force, overwrite=True)

            # Virial stress tensor
            if hasattr(result, 'virial'):
                at.params['virial'] = result.virial
            else:
                at.params['virial'] = fzeros((3,3))

            # Add popn_calculate output and force components
            for k in result.properties.keys():
                if k.startswith('popn_') or k.startswith('force_'):
                    at.add_property(k, getattr(result, k), overwrite=True)
        else:
            print 'Test mode: setting energy, force and virial to zeros'
            at.params['energy'] = 0.
            at.add_property('force', 0., n_cols=3, overwrite=True)
            at.params['virial'] = fzeros((3,3))



def read_formatted_potential(filename):
    """Load a potential write by CASTEP pot_write_formatted() routine, and convert
    to a 3-dimensional array suitable for writing to a .cube file."""

    # Read header lines
    header = []
    fh = open(filename)
    if fh.readline().startswith('BEGIN header'):
        line = fh.readline()
        while not line.startswith('END header'):
            header.append(line)
            line = fh.readline()
    else:
        fh.seek(0)

    # nsp_lines = [line for line in header if line.strip().endswith('! nsp') ]
    # if len(nsp_lines) != 1:
    #    raise IOError('Badly formatted potential file - missing "nsp" header line')
    # nsp = int(nsp_lines[0].split()[0])

    # num_ions_in_species_lines = [line for line in header if line.find('num_ions_in_species') != -1 ]
    # if len(num_ions_in_species_lines) != nsp:
    #    raise IOError('Badly formatted potential file - len(num_ions_in_species_lines) %d != nsp %d' % (len(num_ions_in_species_lines), nsp))
    # num_ions_in_species = [ int(line.split()[0]) for line in num_ions_in_species_lines ]

    # # One line per ion with charge, potential at ion and efield at ion
    # pot_ions = []
    # efield_ions = []
    # for sp in range(nsp):
    #    for i in range(num_ions_in_species[sp]):
    #       line = fh.readline()
    #       species, index, charge, potential, efield1, efield2, efield3 = line.split()
    #       pot_ions.append(float(potential))
    #       efield_ions.append([float(efield1), float(efield2), float(efield3)])

    # pot_ions = farray(pot_ions)
    # efield_ions = farray(efield_ions).T

    # Rest of file is potential data grid
    pot = np.loadtxt(fh)
    nx, ny, nz = pot[:,0].max(), pot[:,1].max(), pot[:,2].max()
    data = np.zeros((nx,ny,nz))
    for (i,j,k,value) in pot:
        data[int(i)-1,int(j)-1,int(k)-1] = value

    return data

def read_formatted_density(filename):
    """Load a potential write by CASTEP pot_write_formatted() routine, and convert
    to a 3-dimensional array suitable for writing to a .cube file."""

    den = np.loadtxt(filename, skiprows=11)
    nx, ny, nz = den[:,0].max(), den[:,1].max(), den[:,2].max()
    data = np.zeros((2,nx,ny,nz))
    for (i,j,k,charge,spin) in den:
        data[:,int(i)-1,int(j)-1,int(k)-1] = charge, spin
    return data

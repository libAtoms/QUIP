import sys, string, os, operator, itertools, logging
from ordereddict import OrderedDict
from farray import *
from quippy import Atoms, Dictionary, HARTREE, BOHR, GPA, atomic_number_from_symbol
from quippy import AtomsReaders, AtomsWriters, atoms_reader

import xml.dom.minidom

# Valid CELL and PARAMETER keywords. Generated using get_valid_keywords() in this module with CASTEP 4.3.

valid_cell_keywords = ['lattice_cart', 'lattice_abc', 'positions_frac', 'positions_abs', 
                       'symmetry_generate', 'symmetry_tol', 'ionic_constraints', 'fix_com', 
                       'cell_constraints', 'external_pressure', 'fix_all_ions', 'fix_all_cell',
                       'species_mass', 'species_pot', 'ionic_velocities', 'species_lcao_states',
                       'kpoints_list', 'kpoints_mp_grid', 'kpoints_mp_spacing', 'kpoints_mp_offset',
                       'kpoint_list', 'kpoint_mp_grid', 'kpoint_mp_spacing', 'kpoint_mp_offset', 
                       'bs_kpoint_path', 'bs_kpoint_path_spacing', 'bs_kpoint_list', 'bs_kpoint_mp_grid',
                       'bs_kpoint_mp_spacing', 'bs_kpoint_mp_offset', 'bs_kpoints_path', 'bs_kpoints_path_spacing',
                       'bs_kpoints_list', 'bs_kpoints_mp_grid', 'bs_kpoints_mp_spacing', 'bs_kpoints_mp_offset',
                       'phonon_supercell_matrix', 'phonon_kpoint_path', 'phonon_kpoint_path_spacing',
                       'phonon_kpoint_list', 'phonon_kpoint_mp_grid', 'phonon_kpoint_mp_offset', 
                       'phonon_kpoint_mp_spacing', 'phonon_gamma_directions', 'phonon_kpoints_path',
                       'phonon_kpoints_path_spacing', 'phonon_kpoints_list', 'phonon_fine_kpoint_list', 
                       'phonon_fine_kpoint_path', 'phonon_fine_kpoint_path_spacing', 'phonon_fine_kpoint_mp_grid',
                       'phonon_fine_kpoint_mp_spacing', 'phonon_fine_kpoint_mp_offset', 'optics_kpoints_list', 
                       'optics_kpoints_mp_grid', 'optics_kpoints_mp_spacing', 'optics_kpoints_mp_offset',
                       'optics_kpoint_list', 'optics_kpoint_mp_grid', 'optics_kpoint_mp_spacing',
                       'optics_kpoint_mp_offset', 'magres_kpoint_list', 'magres_kpoint_path', 
                       'magres_kpoint_path_spacing', 'magres_kpoint_mp_grid', 'magres_kpoint_mp_spacing',
                       'magres_kpoint_mp_offset', 'positions_frac_product', 'positions_abs_product', 
                       'positions_frac_intermediate', 'positions_abs_intermediate', 'fix_vol',
                       'species_gamma', 'species_q', 'supercell_kpoints_list', 'supercell_kpoints_mp_grid', 
                       'supercell_kpoints_mp_spacing', 'supercell_kpoints_mp_offset', 'supercell_kpoint_list',
                       'supercell_kpoint_mp_grid', 'supercell_kpoint_mp_spacing', 'supercell_kpoint_mp_offset', 
                       'supercell_matrix', 'nonlinear_constraints', 'external_efield', 'positions_noise',
                       'cell_noise', 'hubbard_u', 'hubbard_alpha', 'quantisation_axis', 'quantization_axis']

valid_parameters_keywords = ['comment', 'iprint', 'continuation', 'reuse', 'checkpoint', 'task',
                             'calculate_stress', 'calculate_elf', 'num_backup_iter', 'print_clock',
                             'print_memory_usage', 'write_formatted_potential', 'write_formatted_density', 
                             'write_formatted_elf', 'write_orbitals', 'calc_molecular_dipole', 'cml_output',
                             'cml_filename', 'stop', 'xc_functional', 'ppd_integral', 'nlxc_ppd_integral',
                             'nlxc_div_corr_on', 'pspot_nonlocal_type', 'basis_precision', 'fixed_npw', 
                             'finite_basis_corr', 'nelectrons', 'charge', 'spin', 'nup', 'ndown',
                             'spin_polarized', 'spin_polarised', 'nbands', 'electronic_minimizer',
                             'elec_method', 'metals_method', 'elec_energy_tol', 'elec_eigenvalue_tol',
                             'elec_force_tol', 'fix_occupancy', 'elec_dump_file', 'num_dump_cycles', 
                             'elec_restore_file', 'mixing_scheme', 'popn_calculate', 'popn_bond_cutoff',
                             'pdos_calculate_weights', 'bs_max_iter', 'bs_nbands', 'bs_eigenvalue_tol',
                             'bs_xc_functional', 'geom_method', 'geom_max_iter', 'geom_energy_tol', 
                             'geom_force_tol', 'geom_disp_tol', 'geom_stress_tol', 'geom_modulus_est',
                             'geom_frequency_est', 'md_num_iter', 'md_delta_t', 'md_ensemble', 'md_use_pathint', 
                             'md_num_beads', 'md_pathint_staging', 'md_pathint_num_stages', 'md_temperature', 
                             'md_thermostat', 'md_barostat', 'md_cell_t', 'md_langevin_t', 'md_extrap',
                             'md_extrap_fit', 'md_damping_scheme', 'md_opt_damped_delta_t', 'md_elec_force_tol', 
                             'md_sample_iter', 'md_eqm_method', 'md_eqm_ion_t', 'md_eqm_cell_t', 
                             'md_eqm_t', 'optics_xc_functional', 'tssearch_method', 'tssearch_lstqst_protocol', 
                             'tssearch_force_tol', 'tssearch_disp_tol', 'phonon_const_basis', 'phonon_energy_tol', 
                             'phonon_preconditioner', 'phonon_use_kpoint_symmetry', 'phonon_calculate_dos', 
                             'phonon_dos_spacing', 'phonon_dos_limit', 'phonon_finite_disp', 
                             'phonon_force_constant_cutoff', 'phonon_fine_method', 'phonon_method', 
                             'secondd_method', 'efield_energy_tol', 'thermo_t_start', 
                             'thermo_t_stop', 'thermo_t_spacing', 'wannier_spread_tol', 
                             'wannier_sd_step', 'wannier_spread_type', 'wannier_min_algor', 
                             'wannier_ion_rmax', 'wannier_ion_cut_fraction', 'wannier_restart', 
                             'wannier_ion_cut_tol', 'magres_task', 'magres_method', 
                             'magres_conv_tol', 'magres_xc_functional', 'magres_jcoupling_task', 
                             'ga_pop_size', 'ga_max_gens', 'ga_mutate_rate', 
                             'ga_mutate_amp', 'ga_fixed_n', 'ga_bulk_slice', 
                             'calculate_densdiff', 'run_time', 'backup_interval', 
                             'length_unit', 'mass_unit', 'time_unit', 
                             'charge_unit', 'energy_unit', 'force_unit', 
                             'velocity_unit', 'pressure_unit', 'inv_length_unit', 
                             'frequency_unit', 'force_constant_unit', 'volume_unit', 
                             'ir_intensity_unit', 'dipole_unit', 'efield_unit', 
                             'page_wvfns', 'data_distribution', 'opt_strategy', 
                             'opt_strategy_bias', 'num_farms', 'num_proc_in_smp', 
                             'page_ex_pot', 'nlxc_page_ex_pot', 'ppd_size_x', 
                             'nlxc_ppd_size_x', 'ppd_size_y', 'nlxc_ppd_size_y', 
                             'ppd_size_z', 'nlxc_ppd_size_z', 'impose_trs', 
                             'nlxc_impose_trs', 'exchange_reflect_kpts', 'nlxc_exchange_reflect_kpts', 
                             'k_scrn_den_function', 'nlxc_k_scrn_den_function', 'k_scrn_averaging_scheme', 
                             'nlxc_k_scrn_averaging_scheme', 're_est_k_scrn', 'nlxc_re_est_k_scrn', 
                             'calc_full_ex_pot', 'nlxc_calc_full_ex_pot', 'cut_off_energy', 
                             'basis_de_dloge', 'finite_basis_npoints', 'finite_basis_spacing', 
                             'nspins', 'nextra_bands', 'perc_extra_bands', 
                             'elec_temp', 'max_sd_steps', 'max_cg_steps', 
                             'max_diis_steps', 'elec_convergence_win', 'max_scf_cycles', 
                             'spin_fix', 'smearing_scheme', 'smearing_width', 
                             'efermi_tol', 'num_occ_cycles', 'mix_history_length', 
                             'mix_charge_amp', 'mix_spin_amp', 'mix_cut_off_energy', 
                             'mix_metric_q', 'bs_max_cg_steps', 'bs_nextra_bands', 
                             'bs_perc_extra_bands', 'bs_re_est_k_scrn', 'geom_convergence_win', 
                             'geom_spin_fix', 'geom_linmin_tol', 'md_ion_t', 
                             'md_nhc_length', 'md_nose_t', 'md_damping_reset', 
                             'md_elec_energy_tol', 'md_elec_eigenvalue_tol', 'md_elec_convergence_win', 
                             'optics_nextra_bands', 'optics_perc_extra_bands', 'optics_nbands', 
                             'tssearch_qst_max_iter', 'tssearch_cg_max_iter', 'phonon_max_cg_steps', 
                             'phonon_max_cycles', 'phonon_convergence_win', 'phonon_calc_lo_to_splitting', 
                             'phonon_sum_rule', 'calculate_born_charges', 'born_charge_sum_rule', 
                             'efield_max_cg_steps', 'efield_max_cycles', 'efield_convergence_win', 
                             'efield_calc_ion_permittivity', 'efield_ignore_molec_modes', 'efield_freq_spacing', 
                             'efield_oscillator_q', 'thermo_calculate_helmholtz', 'thermo_t_npoints', 
                             'wannier_max_sd_steps', 'wannier_print_cube', 'wannier_ion_moments', 
                             'wannier_ion_cut', 'wannier_ion_cmoments', 'magres_max_cg_steps', 
                             'magres_convergence_win', 'magres_max_sc_cycles', 'magres_write_response', 
                             'excited_state_scissors', 'rand_seed', 'num_proc_in_smp_fine', 
                             'message_size', 'xc_vxc_deriv_epsilon', 'nlxc_div_corr_s_width', 
                             'nlxc_div_corr_tol', 'nlxc_div_corr_npts_step', 'pspot_beta_phi_type', 
                             'grid_scale', 'fine_grid_scale', 'fine_gmax', 
                             'mix_charge_gmax', 'mix_spin_gmax', 'devel_code',
                             'max_scf_cycles_dm', 'max_scf_cycles_edft']

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
      "Read a CASTEP .cell file. cellfile can be a filename or an open file"

      if type(cellfile) == type(''):
         cellfile = open(cellfile,'r')

      current_block = None

      for line in cellfile:
         line = line.strip()

         # Skip comments and blank lines
         if line.startswith('#') or line == '':
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
            cellfile.write('%s  %s\n' % (key, value))


   def to_atoms(self):
      if not self.has_key('LATTICE_CART'):
         raise ValueError('self is missing LATTICE_CART block')

      if not self.has_key('POSITIONS_ABS'):
         raise ValueError('cell is missing POSITIONS_ABS block')
      
      atoms = Atoms(n=len(self['POSITIONS_ABS']),lattice = 
                farray([ [float(x) for x in row] for row in map(string.split, self['LATTICE_CART'])]))
      
      field_list = [line.split() for line in self['POSITIONS_ABS']]

      elements = map(operator.itemgetter(0), field_list)

      # Look up names of elements specified by atomic number
      elements = [ not el.isdigit() and atomic_number_from_symbol(el) or el for el in elements ]

      # Set the element and pos data
      atoms.set_atoms(elements)
      atoms.pos[:,:] = farray([ [float(x) for x in row] \
                                for row in [field[1:4] for field in field_list]])

      return atoms

   def update_from_atoms(self, at):
      
      # Add lattice to cell
      self['LATTICE_CART'] = []
      for i in frange(3):
         self['LATTICE_CART'].append('%f %f %f' % tuple(at.lattice[:,i]))

      # Add atomic positions to cell
      self['POSITIONS_ABS'] = []
      for i in frange(at.n):
         self['POSITIONS_ABS'].append(at.species[:,i].stripstrings() +' %f %f %f' % tuple(at.pos[:,i]))


   @staticmethod
   @atoms_reader('cell', False)
   def cellreader(source):
      self = CastepCell(source)
      yield self.to_atoms()

class CastepCellWriter(object):
   def __init__(self, dest):
      self.dest = dest

   def write(self, at):
      CastepCell(atoms=at).write(self.dest)

AtomsWriters['cell'] = CastepCellWriter

class CastepParam(OrderedDict):
   "Class to wrap a CASTEP parameter (.param) file"

   def __init__(self, paramfile=None, xml=None):
      OrderedDict.__init__(self)
      if paramfile is not None:
         self.read(paramfile)
      elif xml is not None:
         self.read_xml(xml)

   def copy(self):
      new = CastepParam()
      new.update(self)
      return new
   
   def read(self, paramfile):
      "Read a CASTEP .param file. paramfile can be a filename or an open file"

      if type(paramfile) == type(''):
         paramfile = open(paramfile,'r')

      for line in paramfile:
         line = line.strip()

         # Skip comments and blank lines
         if line.startswith('#') or line == '':
            continue

         fields = line.split()
         fields = [ f for f in fields if not f in (':','=')]
         
         key = fields[0].lower()
         if not key in valid_parameters_keywords:
            raise ValueError('Unknown PARAMETERS keyword %s' % key)
         self[key] = ' '.join(fields[1:])

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

      value_map = {'T': 'true', 'F': 'false'}

      if type(castep_output) == type(''):
         castep_output = open(castep_output, 'r')
         castep_output = f.readlines()
         f.close()
      elif hasattr(castep_output, 'read'):
         castep_output = castep_output.readlines()
      
      # Bizarrelly, the parameter names in the output file have the underscores
      # removed from them. Here we construct a mapping from the .castep names to
      # the true .param file
      param_lookup = dict(zip(map(lambda s: s.replace('_',''),valid_parameters_keywords),
           valid_parameters_keywords))
      
      # Find the user parameters section
      try:
         user_param_start = castep_output.index(' ******************************* User Parameters *******************************\n')
      except ValueError:
         raise ValueError('No user parameters found in castep output')
      
      param_lines = []
      i = user_param_start + 1
      while castep_output[i].strip():
         line = castep_output[i]
         key, value = map(string.strip, line[:line.index('#')].split(':',1))
         value = value_map.get(value, value)
         if not key in param_lookup:
            raise ValueError('Unknown parameter %s in castep output file' % key)
         param_lines.append('%s = %s' % (param_lookup[key], value))
         i += 1

      self.read(param_lines)

   def write(self,paramfile=sys.stdout):
      "Write CASTEP .param file"

      if type(paramfile) == type(''):
         paramfile = open(paramfile,'w')
         
      for key, value in self.iteritems():
         paramfile.write('%s = %s\n' % (key, value))


@atoms_reader('geom', False)
def CastepGeomReader(source, atoms_ref=None):
   """Generator to read frames from CASTEP .geom file"""

   if type(source) == type(''):
      source = open(source, 'r')
   source = iter(source)

   if atoms_ref is not None and not atoms_ref.has_property('frac_pos'):
      atoms_ref.add_property('frac_pos', 0.0, n_cols=3)
      atoms_ref.frac_pos[:] = numpy.dot(atoms_ref.g, atoms_ref.pos)

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

      # First line is the time/step
      params['time'] = float(lines[0])

      # Then the energy, in Hartree
      energy_lines = filter(lambda s: s.endswith('<-- E'), lines)
      if len(energy_lines) != 1:
         raise ValueError('Number of energy lines should be exactly one.')

      params['energy'], params['hamiltonian'] = \
                             [float(x)*HARTREE for x in energy_lines[0].split()[0:2]]

      # Lattice is next, in units of Bohr
      lattice_lines = filter(lambda s: s.endswith('<-- h'), lines)
      if lattice_lines != []:
         lattice = farray([ [float(x)* BOHR for x in row[0:3]]
                            for row in map(string.split, lattice_lines) ])
      else:
         if atoms_ref is None:
            raise ValueError('No lattice in .geom file and atoms_ref not present')
         else:
            lattice = atoms_ref.lattice[:]

      # Then optionally virial tensor - convert stress tensor to libAtoms units, then
      # divide by cell volume
      stress_lines  = filter(lambda s: s.endswith('<-- S'), lines)
      if stress_lines:
         virial = farray([ [float(x)*(HARTREE/(BOHR**3)) for x in row[0:3]]
                           for row in map(string.split, stress_lines) ])

      # Find positions and forces
      poslines   = filter(lambda s: s.endswith('<-- R'), lines)
      forcelines = filter(lambda s: s.endswith('<-- F'), lines)

      if poslines != [] and forcelines != [] and len(poslines) != len(forcelines):
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
         at.set_lattice(lattice)
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
      if poslines != []:
         for i, line in fenumerate(poslines):
            el, num, x, y, z, arrow, label = line.split()
            num = int(num)
            if not (el,num) in lookup:
               lookup[(el,num)] = i
            at.z[lookup[(el,num)]] = atomic_number_from_symbol(el)
            at.pos[:,lookup[(el,num)]] = [ float(r)* BOHR for r in (x,y,z)]

         at.set_atoms(at.z) # set at.species property to match at.z
      else:
         at.pos[:] = numpy.dot(at.lattice, at.frac_pos)
                
      # And finally the forces, which are in units of hartree/bohr
      if forcelines != []:
         at.add_property('force', 0.0, n_cols=3)
         for i, line in fenumerate(forcelines):
            el, num, fx, fy, fz, arrow, label = line.split()
            num = int(num)
            at.force[:,lookup[(el,num)]] = [ float(f)*HARTREE/BOHR for f in (fx, fy, fz) ]

      if stress_lines:
         at.params['virial'] = -virial*at.cell_volume()

      yield at

@atoms_reader('castep', False)
@atoms_reader('castep_log', False)
def CastepOutputReader(castep_file, atoms_ref=None, abort=True, save_params=False):
   """Parse .castep file, and return Atoms object with positions,
      energy, forces, and possibly stress and atomic populations as
      well"""

   if type(castep_file) == type(''):
      castep_file = open(castep_file,'r')
   castep_file = iter(castep_file)

   param = CastepParam()

   if atoms_ref is not None and not atoms_ref.has_property('frac_pos'):
      atoms_ref.add_property('frac_pos', 0.0, n_cols=3)
      atoms_ref.frac_pos[:] = numpy.dot(atoms_ref.g, atoms_ref.pos)

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

         if line == ' |      CCC   AA    SSS  TTTTT  EEEEE  PPPP        |\n':
            if got_header:
               got_param = True
               break
            else:
               got_header = True

         if line.startswith(' Starting BFGS iteration'):
            break

         if line.startswith(' BFGS: improving iteration'):
            break
         
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
         pass
         #if abort:
         #   raise

      # Next let's extract the lattice and atomic positions
      lattice_lines = [i for (i,x) in enumerate(castep_output) if x == '                                      Unit Cell\n']
      
      if lattice_lines == []:
         if abort:
            raise ValueError('No unit cell found in castep file')
         else:
            continue
         
      lattice_line = lattice_lines[-1] # last lattice
      
      lattice_lines = castep_output[lattice_line+3:lattice_line+6]
      lattice = fzeros((3,3))
      lattice[:,1] = map(float, lattice_lines[0].split()[0:3])
      lattice[:,2] = map(float, lattice_lines[1].split()[0:3])
      lattice[:,3] = map(float, lattice_lines[2].split()[0:3])

      cell_contents = [i for (i,x) in  enumerate(castep_output) if x == '                                     Cell Contents\n']
      if cell_contents == []:
         if atoms_ref is None or abort:
            raise ValueError('No cell contents found in castep file - try passing atoms_ref')
         else:
            logging.warning('No cell contents. If this is a variable cell geometry optimisation with fixed ions this is normal')
            n_atoms = atoms_ref.n

      if cell_contents != []:
         cell_first_line = cell_contents[-1] # last cell contents line

         try:
            n_atoms = int(castep_output[cell_first_line+2].split()[-1])
            offset = 10
         except IndexError:
            offset = 7

      if atoms_ref is not None:
         # If we were passed in an Atoms object, construct mapping from
         # CASTEP (element, number) to original atom index
         assert n_atoms == atoms_ref.n, 'Number of atoms must match atoms_ref'
         atoms = atoms_ref.copy()
         atoms.set_lattice(lattice)
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

      if cell_contents != []:
         cell_lines = castep_output[cell_first_line+offset:cell_first_line+offset+n_atoms]

         # Fill in species and fractional positions
         atoms.add_property('frac_pos',0.0,n_cols=3)
         for i, line in fenumerate(cell_lines):
            x1, el, num, u, v, w, x2 = line.split()
            num = int(num)
            if not (el,num) in lookup:
               lookup[(el,num)] = i
            atoms.z[lookup[(el,num)]] = atomic_number_from_symbol(el)
            atoms.frac_pos[:,lookup[(el,num)]] = map(float, (u,v,w))

         atoms.set_atoms(atoms.z) # set at.species from at.z

      # Calculate cartesian postions from fractional positions
      # (correct if we're doing a variable cell geom. opt. with fixed ions)
      atoms.pos[:] = numpy.dot(atoms.lattice, atoms.frac_pos)

      if param.has_key('finite_basis_corr') and param['finite_basis_corr'].lower() != 'none':
         energy_lines = filter(lambda s: s.startswith(' Total energy corrected for finite basis set'), \
                               castep_output)
      elif param.has_key('task') and (param['task'].lower() == 'geometryoptimization' or
                                      param['task'].lower() == 'geometryoptimisation'):
         energy_lines = filter(lambda s: s.startswith(' BFGS: Final Enthalpy'), castep_output)
      else:
         energy_lines = filter(lambda s: s.startswith('Final energy') and not s.endswith('<- EDFT\n'), castep_output)

      if (len(energy_lines) == 0):
         if abort:
            raise ValueError('No total energy found in castep file')
      else:
         # Energy is second to last field on line (last is "eV")
         # Use last matching energy line in file
         atoms.params['energy'] = float(energy_lines[-1].split()[-2])

      try:

         for fn in ('Forces', 'Symmetrised Forces'):
            force_start_lines = [i for (i,s) in enumerate(castep_output) if s.find('****** %s ******' % fn) != -1]
            if force_start_lines != []: break

         if force_start_lines == []:
            raise ValueError

         # Use last set of forces
         force_start = force_start_lines[-1]

         # Extract force lines from .castep file
         force_lines = castep_output[force_start+6:force_start+6+atoms.n]

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

      # Have we calculated stress?
      got_virial = False
      if 'calculate_stress' in param and param['calculate_stress'].lower() == 'true':
         try:
            for sn in ('Stress Tensor', 'Symmetrised Stress Tensor'):
               stress_start_lines = [i for i,s in enumerate(castep_output) if s.find('****** %s ******' % sn) != -1 ]
               if stress_start_lines != []: break

            if stress_start_lines == []:
               raise ValueError

            stress_start = stress_start_lines[-1]
            stress_lines = castep_output[stress_start+6:stress_start+9]
            virial = fzeros((3,3),float)
            for i, line in fenumerate(stress_lines):
               star1, label, vx, vy, vz, star2 = line.split()
               virial[:,i] = [-float(v) for v in (vx,vy,vz) ]

            got_virial = True

         except ValueError:
            if abort:
               raise ValueError('No stress tensor found in .castep file')

      # Have we calculated local populations and charges?
      if 'popn_calculate' in param and param['popn_calculate'].lower() == 'true':
         try:
            popn_start = castep_output.index('     Atomic Populations\n')

            popn_lines = castep_output[popn_start+4:popn_start+4+atoms.n]

            atoms.add_property('popn_s',0.0)
            atoms.add_property('popn_p',0.0)
            atoms.add_property('popn_d',0.0)
            atoms.add_property('popn_f',0.0)
            atoms.add_property('popn_total',0.0)
            atoms.add_property('popn_charge',0.0)

            for line in popn_lines:
               el, num_str, s, p, d, f, tot, charge = line.split()
               num = int(num_str)
               atoms.popn_s[lookup[(el,num)]] = float(s)
               atoms.popn_p[lookup[(el,num)]] = float(p)
               atoms.popn_d[lookup[(el,num)]] = float(d)
               atoms.popn_f[lookup[(el,num)]] = float(f)
               atoms.popn_total[lookup[(el,num)]] = float(tot)
               atoms.popn_charge[lookup[(el,num)]] = float(charge)

         except ValueError:
            if abort:
               raise ValueError('No populations found in castep file')

      if save_params:
         atoms.params.update(param)

      if run_time is not None:
         atoms.params['castep_run_time'] = run_time

      # Convert virial to libAtoms units and add to atoms.params
      if got_virial:
         atoms.params['virial'] = virial*atoms.cell_volume()/GPA

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

   for level in ('basic', 'inter', 'expert'):
      lines = os.popen(castep % ('-help %s' % level)).readlines()
      try:
         cell_start = lines.index('Help information on CELL keywords:\n')
         param_start = lines.index('Help information on PARAMETERS keywords:\n')
      except ValueError:
         raise ValueError('Error parsing output of castep -help %s' % level)

      cell_lines = lines[cell_start+2:param_start-2]
      param_lines = lines[param_start+2:-1]

      for lines, keywords  in zip((cell_lines, param_lines),
                                  (valid_cell_keywords,valid_parameters_keywords)):
         lines = map(string.strip, lines)
         lines = filter(lambda s: not s.startswith('*!'), lines)
         for line in lines:
            commentpos = line.find(' *!')
            if commentpos != -1:
               line = line[:commentpos]
            keywords.append(line.lower())

   print 'valid_cell_keywords = %r\n' % valid_cell_keywords
   print 'valid_parameters_keywords = %r' % valid_parameters_keywords
   
   return (valid_cell_keywords, valid_parameters_keywords) 

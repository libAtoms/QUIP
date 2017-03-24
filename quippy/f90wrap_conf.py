# mapping from Fortran kinds to C types
kind_map  = {
 'real':     {'8':   'double',
              'dp':  'double',
	      'DP':  'double',
              '16':  'long_double',
              'qp':  'double'}, # qp is usually double, not quad...
 'complex' : {'8':   'complex_double',
              'dp':  'complex_double'},
 'integer' : {'':       'int',
 	           '8':      'long_long',
	           'dp':     'long_long',
	           'c_intptr_t': 'long_long'}
}

# special cases for types that need initialisation in Fortran
init_lines = {
  'atoms': (('atoms_types_module', ('atoms_repoint',) ),
            ('if (associated(%(PTR)s)) call atoms_repoint(%(PTR)s)',
             'if (present(%(ARG)s)) call atoms_repoint(%(PTR)s)'))
}

# Mapping of Fortran type names to Python classes
class_names = {
}

py_mod_names = {
  # libAtoms modules (low-level)
  'system_module': 'system',
  'units_module': 'units',
  'linearalgebra_module': 'linearalgebra',
  'mpi_context_module': 'mpi_context',
  'quaternions_module': 'quaternions',
  'connection_module': 'connection',
  'clusters_module': 'clusters',
  'structures_module': 'structures',
  'domaindecomposition_module': 'domaindecomposition',
  'paramreader_module': 'paramreader',
  'spline_module': 'spline',
  'frametools_module': 'frametools',
  'topology_module': 'topology',
  'find_surface_atoms_module': 'find_surface_atoms',
  'ringstat_module': 'ringstat',
  'angular_functions_module': 'angular_functions',
  'steinhardt_nelson_qw_module': 'steinhard_nelson_qw',
  'nye_tensor_module': 'nye_tensor',

  # modules which are extended in Python
  'periodictable_module': '_periodictable',
  'table_module': '_table',
  'potential_module': '_potential',
  'dictionary_module': '_dictionary',
  'dynamicalsystem_module': '_dynamicalsystem',
  'cinoutput_module': '_cinoutput',
  'atoms_module': '_atoms',
  'atoms_types_module': 'atoms_types',
  'extendable_str_module': '_extendable_str',
  'structures_module': '_structures',
  'elasticity_module': '_elasticity'
}


# dictionary mapping full names of derived types to abbreviations used when naming methods
short_names = {
  'dynamicalsystem':'ds',
  'potential': 'pot'
}

joint_modules = {
  'atoms_types_module': ('atoms_module', 'connection_module')
}

remove_optional_arguments = ['error']

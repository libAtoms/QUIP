link_quip = [
    libAtoms,
    fox,
    GAP,
    Potentials,
    Utils,
    f90wrap_stub_lib,  # Provide f90wrap_abort for standalone programs
    ]

gap_fit = library('gap_fit',
  '../GAP/gap_fit_module.F90',
  dependencies:[
    blas_dep,
    ],
  link_with : link_quip,
)

gap_fit_exe = executable('gap_fit','../GAP/gap_fit.F90',
  link_with : link_quip+gap_fit,
  dependencies: [
    blas_dep,
    mpi_dep,
    ],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)

quip = executable('quip', 'quip.F90',
  link_with : link_quip,
  dependencies: [
    blas_dep,
    mpi_dep,
    ],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)

DFTB_to_xml = executable('DFTB_to_xml', 'DFTB_to_xml.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
FC_to_xml = executable('FC_to_xml', 'FC_to_xml.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
NRL_TB_to_xml = executable('NRL_TB_to_xml', 'NRL_TB_to_xml.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
align = executable('align', 'align.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
analyze_md_phonons = executable('analyze_md_phonons', 'analyze_md_phonons.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
calc_n_poles = executable('calc_n_poles', 'calc_n_poles.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
descriptors_wrapper_example = executable('descriptors_wrapper_example', 'descriptors_wrapper_example.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
fix_traj_latest = executable('fix_traj_latest', 'fix_traj_latest.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
get_qw = executable('get_qw', 'get_qw.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
local_E_fd = executable('local_E_fd', 'local_E_fd.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
local_random_search = executable('local_random_search', 'local_random_search.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
md = executable('md', 'md.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
md_gid = executable('md_gid', 'md_gid.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
metapot_test = executable('metapot_test', 'metapot_test.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
order_atoms_as_molecules = executable('order_atoms_as_molecules', 'order_atoms_as_molecules.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
quip_wrapper_example = executable('quip_wrapper_example', 'quip_wrapper_example.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
quip_wrapper_simple_example = executable('quip_wrapper_simple_example', 'quip_wrapper_simple_example.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
quip_wrapper_simple_example_C = executable('quip_wrapper_simple_example_C', 'quip_wrapper_simple_example_C.c',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
randomise_calc = executable('randomise_calc', 'randomise_calc.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
slice_sample = executable('slice_sample', 'slice_sample.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
tabletest = executable('tabletest', 'tabletest.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
test_angular = executable('test_angular', 'test_angular.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
test_cova = executable('test_cova', 'test_cova.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
test_grad_sphericals = executable('test_grad_sphericals', 'test_grad_sphericals.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
test_task_manager = executable('test_task_manager', 'test_task_manager.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
ts_calculation = executable('ts_calculation', 'ts_calculation.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
wrap = executable('wrap', 'wrap.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
xyz_to_SPECSYM_forces = executable('xyz_to_SPECSYM_forces', 'xyz_to_SPECSYM_forces.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)
xyz_to_SPECSYM_start = executable('xyz_to_SPECSYM_start', 'xyz_to_SPECSYM_start.F90',
  link_with : link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],  # Allow undefined f90wrap_abort in libAtoms
)

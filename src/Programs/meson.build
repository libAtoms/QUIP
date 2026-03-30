# Programs meson.build - Executable definitions for QUIP

link_quip = [
    libAtoms,
    fox,
    GAP,
    Potentials,
    Utils,
    f90wrap_stub_lib,  # Provide f90wrap_abort for standalone programs
]

# gap_fit library (special case - needed before gap_fit_exe)
gap_fit_lib = library('gap_fit',
  '../GAP/gap_fit_module.F90',
  dependencies: [blas_dep, mpi_dep],
  link_with: link_quip,
)

# Standard Fortran programs (22 programs, blas only, no MPI)
standard_programs = [
  'DFTB_to_xml',
  'FC_to_xml',
  'NRL_TB_to_xml',
  'align',
  'analyze_md_phonons',
  'calc_n_poles',
  'descriptors_wrapper_example',
  'fix_traj_latest',
  'get_qw',
  'local_E_fd',
  'local_random_search',
  'md',
  'md_gid',
  'metapot_test',
  'order_atoms_as_molecules',
  'quip_wrapper_example',
  'quip_wrapper_simple_example',
  'randomise_calc',
  'slice_sample',
  'tabletest',
  'test_angular',
  'test_cova',
  'test_grad_sphericals',
  'test_task_manager',
  'ts_calculation',
  'wrap',
  'xyz_to_SPECSYM_forces',
  'xyz_to_SPECSYM_start',
]

foreach prog : standard_programs
  executable(prog, prog + '.F90',
    link_with: link_quip,
    dependencies: [blas_dep],
    override_options: ['b_lundef=false'],
  )
endforeach

# MPI-enabled programs
executable('quip', 'quip.F90',
  link_with: link_quip,
  dependencies: [blas_dep, mpi_dep],
  override_options: ['b_lundef=false'],
)

executable('gap_fit', '../GAP/gap_fit.F90',
  link_with: link_quip + [gap_fit_lib],
  dependencies: [blas_dep, mpi_dep],
  override_options: ['b_lundef=false'],
)

# C program
executable('quip_wrapper_simple_example_C', 'quip_wrapper_simple_example_C.c',
  link_with: link_quip,
  dependencies: [blas_dep],
  override_options: ['b_lundef=false'],
)

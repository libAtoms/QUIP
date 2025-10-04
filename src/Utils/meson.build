Utils_F90_sources = [
  'crackparams.F90',
  'cracktools.F90',
  'elasticity.F90',
  'force_machine_learning_module.F90',
  'phonons.F90',
  'real_space_covariance.F90',
  'restraints_constraints_xml.F90',
  'structure_analysis_traj_routines.F90',
  'transition_state.F90',
  'ts_params.F90',
]
Utils = library('Utils',
  Utils_F90_sources,
  dependencies: [
    blas_dep,
    mpi_dep,
    ],
  link_with : [libAtoms,fox,GAP,Potentials],
  link_args: ['-lgomp'],
  )

libAtoms_F90_sources = [
  'f90wrap_stub.F90',  # Stub implementation for standalone builds
  'error.F90',
  'kind_module.F90',
  'System.F90',
  'LinkedList.F90',
  'ExtendableStr.F90',
  'MPI_context.F90',
  'Units.F90',
  'linearalgebra.F90',
  'Quaternions.F90',
  'statistics.F90',
  'Dictionary.F90',
  'Table.F90',
  'ParamReader.F90',
  'PeriodicTable.F90',
  'minimization.F90',
  'Atoms_types.F90',
  'Connection.F90',
  'DomainDecomposition.F90',
  'Atoms.F90',
  'RigidBody.F90',
  'Group.F90',
  'angular_functions.F90',
  'steinhardt_nelson_qw.F90',
  'gamma_functions.F90',
  'Constraints.F90',
  'Thermostat.F90',
  'Barostat.F90',
  'DynamicalSystem.F90',
  'Spline.F90',
  'Sparse.F90',
  'CInOutput.F90',
  'clusters.F90',
  'Structures.F90',
  'frametools.F90',
  'nye_tensor.F90',
  'Topology.F90',
  'Atoms_ll.F90',
  'ringstat.F90',
  'histogram1d.F90',
  'find_surface_atoms.F90',
  'k_means_clustering.F90',
  'SocketTools.F90',
  'partition.F90',
  'libAtoms.F90',
  'libAtoms_misc_utils.F90',
  'libAtoms_utils_no_module.F90',
  'ScaLAPACK.F90',
  'Matrix.F90',
  'task_manager.F90',
]
libAtoms_F77_sources = [
  'lbfgs.F'
]
libAtoms_c_sources = [
  'md5.c',
  'cutil.c',
  'netcdf.c',
  'xyz.c',
  'sockets.c',
]
# Platform-specific link args
if build_machine.system() == 'darwin'
  # macOS: allow undefined symbols (will be resolved from _quippy extension)
  libatoms_link_args = ['-lgomp', '-Wl,-undefined,dynamic_lookup']
else
  # Linux: use standard linking (weak symbols in stub will be overridden)
  libatoms_link_args = ['-lgomp']
endif

libAtoms = library('libAtoms',
  libAtoms_c_sources+libAtoms_F77_sources+libAtoms_F90_sources,
  dependencies: [
    blas_dep,
    mpi_dep,
    ],
  c_args: ['-DPROTOTYPES=1'],
  link_args: libatoms_link_args,
  )

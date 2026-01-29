project('quip','fortran','c',
  version : '0.1',
  default_options : ['warning_level=3'],
  meson_version : '>=1.1'
  )

fortran = meson.get_compiler('fortran')
sizeof_fortran_res = fortran.run('''
program sizeof_fortran_t

  type ptr_type
     type(ptr_type), pointer :: p => NULL()
  end type ptr_type
  type(ptr_type) :: ptr
  integer, allocatable, dimension(:) :: ptr_int

  write (*,*) size(transfer(ptr, ptr_int))

end program sizeof_fortran_t
  ''')
sizeof_fortran_t=sizeof_fortran_res.stdout().strip()

git_version_res = run_command('bin/gitversion',check:true)
add_global_arguments([
'-DGIT_VERSION="'+git_version_res.stdout().strip()+'"',
], language:'fortran')

if get_option('gap')
  gap_version_res = run_command('src/GAP/gapversion',check:true)
  add_global_arguments([
    '-DGAP_VERSION='+gap_version_res.stdout().strip(),
    '-DHAVE_GAP',
    '-DHAVE_TB',
    ], language:'fortran')
endif

# Determine QUIP_ARCH dynamically based on host machine
quip_system = host_machine.system()
quip_cpu = host_machine.cpu_family()
if quip_system == 'darwin'
  quip_os = 'darwin'
elif quip_system == 'linux'
  quip_os = 'linux'
else
  quip_os = quip_system
endif

quip_arch = quip_os + '_' + quip_cpu + '_gfortran'
if get_option('mpi')
  quip_arch += '_openmpi'
endif

add_global_arguments([
  '-ffree-line-length-none',
  '-fno-second-underscore',
  '-fPIC',
  '-fno-realloc-lhs',
  '-DFORTRAN_UNDERSCORE',
  '-DGETARG_F2003',
  '-DGETENV_F2003',
  '-DGFORTRAN',
  '-DQUIP_ARCH="' + quip_arch + '"',
  '-DSIZEOF_FORTRAN_T='+sizeof_fortran_t,
  ], language:'fortran')

if get_option('mpi')
  message('MPI ON')
  mpi_dep = [
    dependency('scalapack'),
    dependency('mpi', language: 'fortran'),
    dependency('mpi', language: 'c'),
  ]
  add_global_arguments([
    '-D_MPI',
    '-DSCALAPACK',
    ],language:'fortran')
else
  message('MPI OFF')
  mpi_dep = dependency('',required: false)
endif

# Shared OpenMP link args for all libraries
openmp_link_args = ['-lgomp']

# Get openblas dependency with fallback for systems without pkg-config support
openblas_dep = dependency('openblas', required: false)
if not openblas_dep.found()
  # Fallback: try to find openblas as a system library
  message('openblas not found via pkg-config, trying system library search')
  openblas_lib = fortran.find_library('openblas', required: true)
  # Assume standard include path for openblas
  openblas_dep = declare_dependency(
    dependencies: openblas_lib,
    compile_args: ['-I/usr/include', '-I/usr/include/openblas'],
  )
endif
c = meson.get_compiler('c')
if c.get_id() == 'clang'
  # Apple clang doesn't support -fopenmp, so we need to filter it out
  # For Fortran, gfortran supports it, so we keep it
  blas_dep = declare_dependency(
    compile_args: ['-I' + openblas_dep.get_variable(pkgconfig: 'includedir')],
    link_args: ['-L' + openblas_dep.get_variable(pkgconfig: 'libdir'), '-lopenblas'],
  )
else
  blas_dep = openblas_dep
endif

subdir('src')

# CMake package config generation for find_package(QUIP) support
cmake = import('cmake')

cmake_conf = configuration_data()
cmake_conf.set('QUIP_VERSION', meson.project_version())
cmake_conf.set('QUIP_HAS_GAP', get_option('gap') ? 'TRUE' : 'FALSE')

cmake.write_basic_package_version_file(
    name: 'QUIP',
    version: meson.project_version(),
    compatibility: 'SameMajorVersion',
    install_dir: get_option('libdir') / 'cmake' / 'QUIP'
)

cmake.configure_package_config_file(
    name: 'QUIP',
    input: 'cmake/QUIPConfig.cmake.in',
    configuration: cmake_conf,
    install_dir: get_option('libdir') / 'cmake' / 'QUIP'
)

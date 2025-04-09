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
    ], language:'fortran')
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
  '-DQUIP_ARCH="darwin_x86_64_gfortran"',
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

blas_dep = dependency('openblas')

subdir('src')


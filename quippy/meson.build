project('quippy', 'c', 'fortran',
  version: '0.9.0',
  meson_version: '>=1.1.0',
  default_options: [
    'buildtype=release',
  ],
)

py = import('python').find_installation(pure: false)
py_dep = py.dependency()

# Get compiler
fortran = meson.get_compiler('fortran')
gfortran = find_program('gfortran')

# Get sizeof_fortran_t
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
sizeof_fortran_t = sizeof_fortran_res.stdout().strip()

# Fortran preprocessor flags
fpp_flags = [
  '-DHAVE_GAP',
  '-DHAVE_TB',
  '-DGFORTRAN',
  '-DGETARG_F2003',
  '-DGETENV_F2003',
  '-DFORTRAN_UNDERSCORE',
  '-DQUIP_ARCH="darwin_arm64_gfortran"',
  '-DSIZEOF_FORTRAN_T=' + sizeof_fortran_t,
]

# Sources to wrap
libatoms_sources = [
  'Atoms_types',
  'Atoms',
  'System',
  'Dictionary',
  'DynamicalSystem',
  'nye_tensor',
]

pot_sources = [
  'Potential',
  'Potential_simple',
  'TB',
]

gap_sources = [
  'descriptors',
]

# Preprocess sources
fpp_files = []
foreach src : libatoms_sources
  fpp_file = custom_target(src + '_fpp',
    input: '../src/libAtoms/' + src + '.F90',
    output: src + '.fpp',
    command: [gfortran, '-E', '-cpp', fpp_flags, '@INPUT@'],
    capture: true,
  )
  fpp_files += fpp_file
endforeach

foreach src : pot_sources
  fpp_file = custom_target(src + '_fpp',
    input: '../src/Potentials/' + src + '.F90',
    output: src + '.fpp',
    command: [gfortran, '-E', '-cpp', fpp_flags, '@INPUT@'],
    capture: true,
  )
  fpp_files += fpp_file
endforeach

foreach src : gap_sources
  fpp_file = custom_target(src + '_fpp',
    input: '../src/GAP/' + src + '.F90',
    output: src + '.fpp',
    command: [gfortran, '-E', '-cpp', fpp_flags, '@INPUT@'],
    capture: true,
  )
  fpp_files += fpp_file
endforeach

# Run f90wrap to generate wrappers
f90wrap = find_program('f90wrap', required: true)

f90wrap_cmd = [
  f90wrap,
  '-m', 'quippy',
  '--f90-mod-name', 'quippy._quippy',
  '@INPUT@',
  '-P',
  '-a', 'c_error_abort',
  '--init-file', meson.current_source_dir() + '/init.py',
  '--move-methods',
  '--shorten-routine-names',
  '-k', meson.current_source_dir() + '/KIND_MAP',
  '-s', meson.current_source_dir() + '/STRING_LENGTHS',
  '-S', '10240',
  '--skip', 'atoms_shallowcopy', 'atoms_initialise_ptr', 'atoms_finalise_multi',
  'potential_initialise_inoutput', 'cplx_2d', 'cplx_2d_array1_finalise',
  'potential_local_e_mix_initialise', 'potential_local_e_mix_finalise',
  '--skip-types', 'spherical_harmonics_type', 'potential_local_e_mix',
  '--force-public', 'set_cutoff',
  '--documentation-plugin', meson.current_source_dir() + '/doc_plugin.py',
  '--py-max-line-length', '120',
  '--f90-max-line-length', '120',
]

f90wrap_output = custom_target('f90wrap',
  input: fpp_files,
  output: [
    'f90wrap_Atoms_types.f90',
    'f90wrap_Atoms.f90',
    'f90wrap_System.f90',
    'f90wrap_Dictionary.f90',
    'f90wrap_DynamicalSystem.f90',
    'f90wrap_nye_tensor.f90',
    'f90wrap_Potential.f90',
    'f90wrap_Potential_simple.f90',
    'f90wrap_TB.f90',
    'f90wrap_descriptors.f90',
  ],
  command: f90wrap_cmd,
  depend_files: [
    'init.py',
    'KIND_MAP',
    'STRING_LENGTHS',
    'doc_plugin.py',
  ],
)

# Patch f90wrap-generated Python files to fix overloaded interface shadowing bug
patch_script = find_program('patch_f90wrap_interfaces.py', required: true)
patched_py_files = custom_target('patch_f90wrap',
  input: f90wrap_output,
  output: 'f90wrap_patched.stamp',
  command: [
    patch_script,
    meson.current_build_dir() / 'quippy' / 'atoms_types_module.py',
    meson.current_build_dir() / 'quippy' / 'atoms_module.py',
    meson.current_build_dir() / 'quippy' / 'dictionary_module.py',
    meson.current_build_dir() / 'quippy' / 'system_module.py',
    meson.current_build_dir() / 'quippy' / 'dynamicalsystem_module.py',
    meson.current_build_dir() / 'quippy' / 'nye_tensor_module.py',
    meson.current_build_dir() / 'quippy' / 'potential_module.py',
    meson.current_build_dir() / 'quippy' / 'potential_simple_module.py',
    meson.current_build_dir() / 'quippy' / 'tb_module.py',
    meson.current_build_dir() / 'quippy' / 'descriptors_module.py',
    '&&', 'touch', '@OUTPUT@'
  ],
  depends: f90wrap_output,
  build_by_default: true,
)

# Create alias target to build all wrappers
alias_target('wrappers', patched_py_files)

# Now we need to compile the wrappers and create the Python extension
# Get paths to QUIP libraries
quip_root = meson.project_source_root() / '..'
quip_builddir = quip_root / 'builddir' / 'src'

# Try multiple methods to find OpenBLAS
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

# Get f2py-f90wrap
f2py_f90wrap = find_program('f2py-f90wrap', required: true)

# Get Python extension suffix
ext_suffix = run_command(py, ['-c',
  'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))'],
  check: true).stdout().strip()

# Find the Fortran module directories (they have platform-specific suffixes)
# We need to find directories containing .mod files
find_mod_dir = '''
import os
import sys
libdir = sys.argv[1]
# Look for subdirectories containing .mod files
for entry in os.listdir(libdir):
    entry_path = os.path.join(libdir, entry)
    if os.path.isdir(entry_path):
        # Check if this directory contains .mod files
        for f in os.listdir(entry_path):
            if f.endswith('.mod'):
                print(entry_path)
                sys.exit(0)
sys.exit(1)
'''

libatoms_moddir = run_command(py, ['-c', find_mod_dir, quip_builddir / 'libAtoms'], check: true).stdout().strip()
fox_moddir = run_command(py, ['-c', find_mod_dir, quip_builddir / 'fox'], check: true).stdout().strip()
gap_moddir = run_command(py, ['-c', find_mod_dir, quip_builddir / 'GAP'], check: true).stdout().strip()
potentials_moddir = run_command(py, ['-c', find_mod_dir, quip_builddir / 'Potentials'], check: true).stdout().strip()
utils_moddir = run_command(py, ['-c', find_mod_dir, quip_builddir / 'Utils'], check: true).stdout().strip()

# Create wrapper script for f2py-f90wrap
# This is needed because we need to move the output file after f2py builds it
# We also need to add the QUIP build directories to the Fortran module search path
f2py_wrapper_script = '''#!/bin/bash
set -e
@F2PY@ --build-dir "@BUILDDIR@/f2py_build" -c -m _quippy "$@" \
  -L"@QUIP_LIBATOMS@" \
  -L"@QUIP_FOX@" \
  -L"@QUIP_GAP@" \
  -L"@QUIP_POTENTIALS@" \
  -L"@QUIP_UTILS@" \
  -llibAtoms -lfox -lGAP -lPotentials -lUtils \
  -L"@OPENBLAS_LIBDIR@" -lopenblas -lgomp -lgfortran \
  --f90flags="-ffree-line-length-none -fPIC -I@LIBATOMS_MODDIR@ -I@FOX_MODDIR@ -I@GAP_MODDIR@ -I@POTENTIALS_MODDIR@ -I@UTILS_MODDIR@" \
  --fcompiler=gfortran

# Patch the C source to export f90wrap_abort_ symbol with visibility attribute
# This is needed so libAtoms can call it when errors occur
cd "@BUILDDIR@/f2py_build"
# Add visibility attribute to the function declarations and definitions
sed -i.bak1 's/void f90wrap_abort_(/__attribute__((visibility("default"))) void f90wrap_abort_(/g' _quippymodule.c
sed -i.bak2 's/void f90wrap_abort__(/__attribute__((visibility("default"))) void f90wrap_abort__(/g' _quippymodule.c

# Patch the meson.build file to add rpath settings
sed -i.bak3 "s|link_args.*: *\[|link_args: ['-Wl,-rpath,@QUIP_LIBATOMS@', '-Wl,-rpath,@QUIP_FOX@', '-Wl,-rpath,@QUIP_GAP@', '-Wl,-rpath,@QUIP_POTENTIALS@', '-Wl,-rpath,@QUIP_UTILS@', '-Wl,-rpath,@OPENBLAS_LIBDIR@', |" meson.build

# Build the extension
cd bbdir && meson compile
# f2py-f90wrap creates the .so in the f2py_build/bbdir directory, move it to output
mv "_quippy@EXT_SUFFIX@" "@OUTPUT@"
'''

f2py_wrapper_script = f2py_wrapper_script.replace('@F2PY@', f2py_f90wrap.full_path())
f2py_wrapper_script = f2py_wrapper_script.replace('@BUILDDIR@', meson.current_build_dir())
f2py_wrapper_script = f2py_wrapper_script.replace('@QUIP_LIBATOMS@', quip_builddir / 'libAtoms')
f2py_wrapper_script = f2py_wrapper_script.replace('@QUIP_FOX@', quip_builddir / 'fox')
f2py_wrapper_script = f2py_wrapper_script.replace('@QUIP_GAP@', quip_builddir / 'GAP')
f2py_wrapper_script = f2py_wrapper_script.replace('@QUIP_POTENTIALS@', quip_builddir / 'Potentials')
f2py_wrapper_script = f2py_wrapper_script.replace('@QUIP_UTILS@', quip_builddir / 'Utils')
f2py_wrapper_script = f2py_wrapper_script.replace('@LIBATOMS_MODDIR@', libatoms_moddir)
f2py_wrapper_script = f2py_wrapper_script.replace('@FOX_MODDIR@', fox_moddir)
f2py_wrapper_script = f2py_wrapper_script.replace('@GAP_MODDIR@', gap_moddir)
f2py_wrapper_script = f2py_wrapper_script.replace('@POTENTIALS_MODDIR@', potentials_moddir)
f2py_wrapper_script = f2py_wrapper_script.replace('@UTILS_MODDIR@', utils_moddir)
f2py_wrapper_script = f2py_wrapper_script.replace('@OPENBLAS_LIBDIR@', openblas_dep.get_variable(pkgconfig: 'libdir'))
f2py_wrapper_script = f2py_wrapper_script.replace('@EXT_SUFFIX@', ext_suffix)
f2py_wrapper_script = f2py_wrapper_script.replace('@OUTPUT@', meson.current_build_dir() / ('_quippy' + ext_suffix))

f2py_wrapper = configure_file(
  output: 'f2py_wrapper.sh',
  command: ['echo', f2py_wrapper_script],
  capture: true,
)

# Build the extension module using f2py-f90wrap via the wrapper script
quippy_ext = custom_target('quippy_extension',
  input: f90wrap_output,
  output: '_quippy' + ext_suffix,
  command: ['bash', f2py_wrapper, '@INPUT@'],
  depends: [f90wrap_output, patched_py_files],
  build_by_default: true,
  install: true,
  install_dir: py.get_install_dir() / 'quippy',
)

# Install the quippy Python source files
py.install_sources(
  [
    'quippy/cli.py',
    'quippy/convert.py',
    'quippy/descriptors.py',
    'quippy/dynamicalsystem.py',
    'quippy/gap_tools.py',
    'quippy/nye_tensor.py',
    'quippy/potential.py',
  ],
  subdir: 'quippy',
)

# Install f90wrap-generated Python wrapper files
# Use custom_target to create a dummy target that depends on f90wrap and patches
# This ensures the files exist and are patched before install
f90wrap_py_files = custom_target('f90wrap_py_install',
  input: f90wrap_output,
  output: 'f90wrap_installed.stamp',
  command: ['touch', '@OUTPUT@'],
  depends: patched_py_files,
  build_by_default: true,
  install: false,
)

# Now install the generated Python files
# These will be picked up by the editable loader
install_subdir(
  meson.current_build_dir() / 'quippy',
  install_dir: py.get_install_dir(),
  exclude_files: ['__pycache__'],
  strip_directory: false,
)

message('f90wrap wrapper generation and extension build configured')

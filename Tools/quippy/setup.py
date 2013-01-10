#!/usr/bin/env python
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import f90doc, f2py_wrapper_gen, patch_f2py, numpy
import sys, os, cPickle, glob, stat, subprocess, re, string, StringIO

from numpy.distutils.core import setup, Extension
from numpy.distutils.ccompiler import gen_preprocess_options
from numpy import get_include
from numpy.distutils.system_info import get_info
from distutils.file_util import copy_file
from distutils.dep_util import newer, newer_group
from custom_commands import *
from distutils.util import get_platform

major, minor = sys.version_info[0:2]
if (major, minor) < (2, 4):
    sys.stderr.write('Python 2.4 or later is needed to use this package\n')
    sys.exit(1)

def find_quip_root_and_arch():
    """Find QUIP root directory."""
    quip_root = os.path.abspath(os.path.join(os.getcwd(), '../../'))
    os.environ['QUIP_ROOT'] = quip_root # export to enviroment for Makefile variables

    if not 'QUIP_ARCH' in os.environ:
        raise ValueError('QUIP_ARCH environment variable not set')

    quip_arch = os.environ['QUIP_ARCH']
    if 'QUIP_ARCH_SUFFIX' in os.environ:
      quip_arch = quip_arch+os.environ['QUIP_ARCH_SUFFIX']

    return (quip_root, quip_arch)


def F90WrapperBuilder(modname, wrap_sources, targets, cpp, sizeof_fortran_t,
                      dep_type_maps=[], kindlines=[], short_names={},
                      initlines={}, filtertypes=None, prefix='',callback_routines={}):
    """Build a Fortran 90 wrapper for the given F95 source files
    that is suitable for use with f2py. Derived types are wrapped to 
    give access to methods and instance data."""

    def func(extension, build_dir):

        # first ensure libraries are up to date
        for (dir, target) in targets:
            command = "cd %s && make %s" % (dir, target)
            print 'Rebuilding target %s with command "%s"' % (target, command)
            if os.system(command) != 0:
                raise SystemExit('Command "%s" failed' % command)

        # Have any of wrap_sources changed since we last scanned source files?
        f90doc_file = os.path.join(build_dir, '../../%s.f90doc' % modname)
        if newer_group(wrap_sources, f90doc_file):

            # Rebuild .f90doc file containing interfaces of wrapped routines
            tmp_wrap_sources = []
            cpp_opt = ' '.join(gen_preprocess_options(macros, include_dirs))                
            for src in wrap_sources:
                tmp_file = os.path.join(build_dir.replace('src', 'temp'), os.path.basename(src))
                if not os.path.exists(os.path.dirname(tmp_file)): os.makedirs(os.path.dirname(tmp_file))
                command = "%s %s %s | grep -v '^#' > %s" % (' '.join(cpp), cpp_opt, src, tmp_file)
                print 'Executing command %s' % command
                os.system(command)
                if os.path.exists(src[:-4]+'.s'): os.remove(src[:-4]+'.s')
                tmp_wrap_sources.append(tmp_file)

            programs, modules, functs, subts = f90doc.read_files(tmp_wrap_sources)
            cPickle.dump((programs, modules, functs, subts), open(f90doc_file, 'w'))
        else:
            # Read previous .f90doc file
            (programs, modules, functs, subts) = cPickle.load(open(f90doc_file))

        # Update map from type names to module in which they are defined
        for mod, name in modules:
            for n in [t.name for t in mod.types]:
                type_map[n.lower()] = mod.name
            
        for item in dep_type_maps:
            if hasattr(item, '__getitem__') and hasattr(item, 'keys'): # dictionary
                type_map.update(item)
            else: # assume it's a string
                type_map.update(cPickle.load(open('%s.type' % item)))

        # Try to load previous .spec file
        res = []
        fortran_spec = {}
        spec_file = os.path.join(build_dir, '../../%s.spec' % modname)
        if os.path.exists(spec_file):
            fortran_spec = cPickle.load(open(spec_file))

        # Write new wrapper files and update .spec file
        wrap_modules = []
        for file in wrap_sources:

            for mod, name in modules:
                if os.path.basename(name) == os.path.basename(file):
                    break
            else:
                raise ValueError("Can't find Fortran module corresponding to file %s" % file)

            wrap_mod_name = mod.name.lower()[:-7]
            wrap_modules.append(wrap_mod_name)

            wrapper = '%s/%s_%s_wrap.f90' % (build_dir, modname, wrap_mod_name)

            if not newer(name, wrapper):
                res.append(wrapper)
                continue

            public_symbols = f2py_wrapper_gen.find_public_symbols(file)

            print 'public_symbols = ', public_symbols

            tmpf = StringIO.StringIO()
            new_spec = f2py_wrapper_gen.wrap_mod(mod, type_map, tmpf, kindlines=kindlines, initlines=initlines,
                                                 filtertypes=filtertypes, prefix=prefix, callback_routines=callback_routines,
                                                 public_symbols=public_symbols, sizeof_fortran_t=sizeof_fortran_t)

            if not os.path.exists(wrapper) or (new_spec[wrap_mod_name] != fortran_spec.get(wrap_mod_name, None)):
                print 'Interface for module %s has changed. Rewriting wrapper file' % mod.name
                fortran_spec.update(new_spec)
                wrapperf = open(wrapper, 'w')
                wrapperf.write(tmpf.getvalue())
                wrapperf.close()
            else:
                print 'Interface for module %s unchanged' % mod.name

            tmpf.close()
            res.append(wrapper)

        fortran_spec['wrap_modules'] = wrap_modules
        fortran_spec['short_names'] = short_names
        fortran_spec['quip_root'] = quip_root
        fortran_spec['quip_arch'] = quip_arch
        fortran_spec['quip_makefile'] = makefile
        cPickle.dump(fortran_spec, open(os.path.join(build_dir, '../../%s.spec' % modname), 'w'))

        import pprint
        spec_py_name = '%s/spec.py' % build_dir
        spec_py = open(spec_py_name, 'w')
        spec_py.write('spec = %s\n' % pprint.pformat(fortran_spec))
        spec_py.close()
        res.append(spec_py_name)

        return res

    return func

def expand_addsuffix(s):
    add_suffix = re.compile(r'\$[\(\{]addsuffix (.*?),(.*?)[\)\}]')
    try:
        m = add_suffix.search(s)
    except TypeError:
        return s
    if m:
        while m is not None:
            suffix, files = m.groups()
            s = add_suffix.sub(' '.join([f + suffix for f in files.split()]),s,1).strip()
            m = add_suffix.search(s)
    return s

def expand_addprefix(s):
    add_prefix =  re.compile(r'\$[\(\{]addprefix (.*?),(.*?)[\}\)]')
    try:
        m = add_prefix.search(s)
    except TypeError:
        return s
    if m:
        while m is not None:
            prefix, files = m.groups()
            s = add_prefix.sub(' '.join([prefix + f for f in files.split()]),s,1).strip()
            m = add_prefix.search(s)
    return s


def read_arch_makefiles_and_environment(quip_root, quip_arch):

    # Write a Makefile which simply includes Makefile.inc, Makefile.rules and Makefile.${QUIP_ARCH}
    f = open('Makefile.quippy', 'w')
    f.write("""ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif
include Makefile.inc
include Makefile.rules""")
    f.close()

    # Dump make database to file
    os.system('make -f Makefile.quippy -I %s -I %s -I %s -p > make.dump' %
              (os.path.join(quip_root, 'Makefiles'),
               os.path.join(quip_root, 'build.%s' % quip_arch),
               quip_root))

    # Parse dumped file
    makefile = parse_makefile('make.dump')

    # Tidy up
    os.remove('make.dump')
    os.remove('Makefile.quippy')

    # Allow environment variables to override contents of Makefiles
    makefile.update(os.environ)

    # Expand prefices and suffices
    for k, v in makefile.iteritems():
        v = expand_addprefix(v)
        v = expand_addsuffix(v)
        makefile[k] = v

    return makefile


def find_wrap_sources(makefile, quip_root):
    source_dirs = []
    wrap_sources = []
    wrap_types   = []
    libraries = []
    targets = []

    libatoms_dir   = os.path.join(quip_root, 'libAtoms/')
    wrap_sources += [os.path.join(libatoms_dir, s) for s in
                     ['System.f95', 'ExtendableStr.f95', 'MPI_context.f95', 'Units.f95', 'linearalgebra.f95', 'Quaternions.f95', 
                     'Dictionary.f95', 'Table.f95', 'PeriodicTable.f95', 'Atoms_types.f95', 'Atoms.f95', 'Connection.f95', 'DynamicalSystem.f95',
                     'clusters.f95','Structures.f95', 'DomainDecomposition.f95', 'CInOutput.f95', 'ParamReader.f95', 'Spline.f95',
		     'frametools.f95', 'Topology.f95', 'find_surface_atoms.f95']]
    wrap_types += ['inoutput', 'mpi_context', 'dictionary', 'table', 'atoms', 'connection', 'quaternion',
                   'dynamicalsystem', 'domaindecomposition', 'cinoutput', 'extendable_str', 'spline']
    source_dirs.append(libatoms_dir)
    libraries.append('atoms')
    targets.append((quip_root, 'libAtoms/libatoms.a'))

    if 'HAVE_GP_PREDICT' in makefile and int(makefile['HAVE_GP_PREDICT']) == 1:
        gp_dir = os.path.join(quip_root, 'GAP_predict')
        source_dirs.append(gp_dir)
        libraries.append('gap_predict')
        targets.extend([(quip_root, 'GAP_predict/libgap_predict.a')])

    if 'HAVE_GP_TEACH' in makefile and int(makefile['HAVE_GP_TEACH']) == 1:
        gp_dir = os.path.join(quip_root, 'GAP_teach')
        source_dirs.append(gp_dir)
        libraries.append('gap_teach')
        targets.extend([(quip_root, 'GAP_teach/libgap_teach.a')])

    quip_core_dir = os.path.join(quip_root, 'QUIP_Core/')
    source_dirs.append(quip_core_dir)
    wrap_sources += [os.path.join(quip_core_dir, s) for s in ['Potential.f95', 'ElectrostaticEmbed.f95', 'AdjustablePotential.f95']]
    wrap_types += ['potential']
    libraries = ['quip_core', 'thirdparty'] + libraries
    targets.append((quip_root, 'QUIP_Core/libquip_core.a'))

    do_tools = not 'QUIPPY_NO_TOOLS' in makefile or ('QUIPPY_NO_TOOLS' in makefile and not int(makefile['QUIPPY_NO_TOOLS']))
    do_crack = not 'QUIPPY_NO_CRACK' in makefile or ('QUIPPY_NO_CRACK' in makefile and not int(makefile['QUIPPY_NO_CRACK']))
       
    if do_tools or do_crack:
        quip_utils_dir = os.path.join(quip_root, 'QUIP_Utils/')
        source_dirs.append(quip_utils_dir)
        libraries = ['quiputils'] + libraries
        targets.append((quip_root, 'QUIP_Utils'))

    if do_tools:
        wrap_sources += [os.path.join(quip_utils_dir, s) for s in ['elasticity.f95', 'real_space_covariance.f95']]
        wrap_types += ['realspacecovariance']

    if do_crack:
        wrap_sources += [os.path.join(quip_utils_dir,f) for f in ('crackparams.f95', 'cracktools.f95')]
        wrap_types += ['crackparams', 'crackmdparams']

    if 'HAVE_CP2K' in makefile and int(makefile['HAVE_CP2K']) == 1:
	quip_filepot_drivers_dir = os.path.join(quip_root, 'QUIP_FilePot_Drivers')
	source_dirs.append(quip_filepot_drivers_dir)
        targets.append((quip_root, 'QUIP_FilePot_Drivers'))
        wrap_sources.append(os.path.join(quip_filepot_drivers_dir, 'cp2k_driver_module.f95'))
	libraries = ['cp2k_driver'] + libraries
        
    return source_dirs, wrap_sources, wrap_types, libraries, targets


type_map = {}

quip_root, quip_arch = find_quip_root_and_arch()

if not os.path.isdir(os.path.join(quip_root, 'libAtoms')):
    raise ValueError('Cannot find libAtoms directory under %s - please set QUIP_ROOT env var' % quip_root)

makefile = read_arch_makefiles_and_environment(quip_root, quip_arch)
makefile_test = lambda var: var in makefile and int(makefile[var])

# Check for essential Makefile variables
for key in ('QUIPPY_FCOMPILER', 'QUIPPY_CPP'):
    if not key in makefile:
        raise ValueError('Mandatory variable %s must be specified in Makefile or environment' % key)

# C preprocessor
cpp = makefile.get('QUIPPY_CPP', 'cpp').split()

# extract include directories from INCLUDES Makefile variable
fields = makefile['INCLUDES'].split()
include_dirs  = [s[2:] for s in fields if s[:2] == '-I']

# extract libraries and library_dirs from SYSLIBS Makefile variable
fields = makefile['SYSLIBS'].split()
libraries = [s.startswith('-l')  and s[2:] or s  for s in fields
             if s.startswith('-l') or s == '-Bstatic' or s == '-Bdynamic' or s.startswith('-Wl') or s.endswith('.a') ]
library_dirs  = [s[2:] for s in fields if s[:2] == '-L']

# everything else in SYSLIBS is an extra link argument
extra_link_args = [s for s in fields if not s[2:] in libraries and not s in libraries and not s[2:] in library_dirs and not s[2:] in include_dirs]

if 'QUIPPY_LDFLAGS' in makefile:
    extra_link_args.extend(makefile['QUIPPY_LDFLAGS'].split())

# Preprocessor macros
macros = [('SVN_VERSION',r'\"%s\"' % os.popen('%s/utility_scripts/svnversion -n .' % quip_root).read())]
for defn in makefile['DEFINES'].split():
    if defn[:2] == '-D':
        if '=' in defn:
            n, v = defn[2:].split('=')
            macros.append((n,v))
        else:
            macros.append((defn[2:], None))
    elif defn[:2] == '-U':
        macros.append(defn[2:])

print 'include_dirs', include_dirs
print 'libraries', libraries
print 'lib_dirs', library_dirs
print 'macros', macros
print 'extra_link_args', extra_link_args

# Default distutils options -- will be overriden by command line options
# once setup() is invoked.
default_options= {
    'config': {
    },
    
    'config_fc': {
    'f90flags':  makefile.get('QUIPPY_F90FLAGS', '').split(),
    'f77flags':  makefile.get('QUIPPY_F77FLAGS', '').split(),
    'debug':     int(makefile.get('QUIPPY_DEBUG',1))
    },

    'build': {
    'build_base': 'build.%s' % quip_arch,
    'debug':     int(makefile.get('QUIPPY_DEBUG',1))
   
    },

    'build_ext': {
    'fcompiler': makefile['QUIPPY_FCOMPILER'],
    'dep_libs': True
    },

    'build_src': {
    'f2py_opts': makefile.get('QUIPPY_F2PY_OPTS', '--no-wrap-functions')
    },

    'clean':{
    'all': True
    },

    'test':{
    'verbosity': 2
    }
}

got_gfortran45 = False
if makefile['QUIPPY_FCOMPILER'] == 'gnu95':
    version = os.popen("gfortran --version | head -1 | sed 's/(.*)//' | awk '{ print $3 }'").read()
    print version
    version = [int(s) for s in version.split('.')[0:2]]
    got_gfortran45 = version[0] == 4 and version[1] == 5

if makefile_test('QUIPPY_DEBUG') or got_gfortran45:
    os.environ['FOPT'] = '-O0'
    os.environ['FARCH'] = ''
    macros.append(('DEBUG',None))

if 'QUIPPY_OPT' in makefile:
    default_options['config_fc']['opt'] = makefile['QUIPPY_OPT'].split()

sizeof_fortran_t = int(makefile['SIZEOF_FORTRAN_T'])

# Install options
if 'QUIPPY_INSTALL_OPTS' in makefile:
    install_opts = makefile['QUIPPY_INSTALL_OPTS'].split()
    default_options['install'] = {}
    for opt in install_opts:
        n, v = opt.split('=')
        n = n[2:] # remove --
        default_options['install'][n] = v

# Find Fortran source code files
source_dirs, wrap_sources, wrap_types, quip_libraries, quip_targets = find_wrap_sources(makefile, quip_root)
include_dirs.extend(source_dirs)
libraries = quip_libraries + libraries

# Add build.${QUIP_ARCH} to include and library paths
include_dirs.append(os.path.join(quip_root, 'build.%s' % quip_arch))
library_dirs.append(os.path.join(quip_root, 'build.%s' % quip_arch))

# arraydata extension module
f2py_info = get_info('f2py')
arraydata_ext = Extension(name='quippy.arraydata', 
                          sources=['arraydatamodule.c'] + f2py_info['sources'],
                          include_dirs=f2py_info['include_dirs']+include_dirs,
                          define_macros=[('SIZEOF_FORTRAN_T', sizeof_fortran_t)])


# _quippy extension module
quippy_ext = Extension(name='quippy._quippy',
                       sources=[ F90WrapperBuilder('quippy',
                                                   wrap_sources=wrap_sources,
                                                   cpp=cpp,
                                                   sizeof_fortran_t=sizeof_fortran_t,
                                                   targets=quip_targets,
                                                   dep_type_maps=[{'c_ptr': 'iso_c_binding',
                                                                   'dictionary_t':'FoX_sax'}], 
                                                   kindlines=['use system_module',
                                                               'use iso_c_binding'],
                                                   short_names={'dynamicalsystem':'ds',
                                                                'potential': 'pot'},
                                                   initlines={'atoms': ('atoms_module', ('if (associated(%(PTR)s)) call atoms_repoint(%(PTR)s)',
                                                                                         'if (present(%(ARG)s)) call atoms_repoint(%(PTR)s)'))},
                                                   filtertypes=wrap_types,
                                                   prefix='qp_',
                                                   callback_routines={'potential_set_callback':{'arglines': ['external :: qp_callback'],
                                                                                                'attributes':[],
                                                                                                'call':'qp_this'},
                                                                      'dynamicalsystem_run':{'arglines' : ['external :: qp_hook'],
                                                                                             'attributes': [],
                                                                                             'call':'qp_this'}})
                                                   
                                 ],
                       library_dirs=library_dirs,
                       include_dirs=include_dirs,
                       libraries=libraries,
                       define_macros= macros,
                       extra_link_args=extra_link_args)

exts = [arraydata_ext, quippy_ext]

# Finally, call setup() to run command

setup(name='quippy',
      packages = ['quippy'],
      py_modules=['qlab'],
      ext_modules = exts,
      scripts=['scripts/quippy'] + glob.glob('scripts/*.py'),
      cmdclass = {'clean': clean, 'test': test, 'build_ext': build_ext, 'interact': interact},
      version=os.popen('%s/utility_scripts/svnversion -n .' % quip_root).read(),
      description='Python bindings to QUIP code',
      author='James Kermode',
      author_email='james.kermode@kcl.ac.uk',
      url='http://www.jrkermode.co.uk/quippy',
      options=default_options)

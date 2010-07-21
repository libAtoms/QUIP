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
# from distutils.sysconfig import parse_makefile
from custom_commands import *
from distutils.util import get_platform

major, minor = sys.version_info[0:2]
if (major, minor) < (2, 4):
    sys.stderr.write('Python 2.4 or later is needed to use this package\n')
    sys.exit(1)

# parse_makefile() from Python 2.6.1 distutils.sys_config.

# Regexes needed for parsing Makefile (and similar syntaxes,
# like old-style Setup files).
_variable_rx = re.compile("([a-zA-Z][a-zA-Z0-9_]+)\s*=\s*(.*)")
_findvar1_rx = re.compile(r"\$\(([A-Za-z][A-Za-z0-9_]*)\)")
_findvar2_rx = re.compile(r"\${([A-Za-z][A-Za-z0-9_]*)}")

def parse_makefile(fn, g=None):
   """Parse a Makefile-style file.

   A dictionary containing name/value pairs is returned.  If an
   optional dictionary is passed in as the second argument, it is
   used instead of a new dictionary.
   """
   from distutils.text_file import TextFile
   fp = TextFile(fn, strip_comments=1, skip_blanks=1, join_lines=1)

   if g is None:
       g = {}
   done = {}
   notdone = {}

   while 1:
       line = fp.readline()
       if line is None:                # eof
           break
       m = _variable_rx.match(line)
       if m:
           n, v = m.group(1, 2)
           v = string.strip(v)
           if "$" in v:
               notdone[n] = v
           else:
               try: v = int(v)
               except ValueError: pass
               done[n] = v

   # do variable interpolation here
   while notdone:
       for name in notdone.keys():
           value = notdone[name]
           m = _findvar1_rx.search(value) or _findvar2_rx.search(value)
           if m:
               n = m.group(1)
               found = True
               if n in done:
                   item = str(done[n])
               elif n in notdone:
                   # get it on a subsequent round
                   found = False
               elif n in os.environ:
                   # do it like make: fall back to environment
                   item = os.environ[n]
               else:
                   done[n] = item = ""
               if found:
                   after = value[m.end():]
                   value = value[:m.start()] + item + after
                   if "$" in after:
                       notdone[name] = value
                   else:
                       try: value = int(value)
                       except ValueError:
                           done[name] = string.strip(value)
                       else:
                           done[name] = value
                       del notdone[name]
           else:
               # bogus variable reference; just drop it since we can't deal
               del notdone[name]

   fp.close()

   # save the results in the global dictionary
   g.update(done)
   return g


def find_quip_root_and_arch():
    """Find QUIP root directory."""
    if 'QUIP_ROOT' in os.environ:
        quip_root = os.environ['QUIP_ROOT']
    else:
        quip_root = os.path.abspath(os.path.join(os.getcwd(), '../../'))

    if not 'QUIP_ARCH' in os.environ:
        raise ValueError('QUIP_ARCH environment variable not set')
    
    quip_arch = os.environ['QUIP_ARCH']

    return (quip_root, quip_arch)


def SourceImporter(infile, defines, include_dirs, cpp):
    """Import source code from infile and copy to build_dir/package,
       passing through filter if it's an F95 source file. The filter
       converts all private variables to public and runs the cpp
       preprocessor."""

    def func(extension, build_dir):

        outfile = os.path.join(build_dir, os.path.split(infile)[1])

        if outfile.endswith('.f95'):
            outfile = outfile[:-4]+'.f90'

        cpp_opt = ' '.join(gen_preprocess_options(macros, include_dirs))

        if infile.endswith('.f95'):
            if newer(infile, outfile):
                print 'filtering %s to create %s' % (infile, outfile)
                # Munge source a little... this could be converted to pure python at some point

                os.system("""perl -e 'while(<>) {
    if (/^(double precision\s+)?\s*(recursive\s+)?(subroutine|function)\s+(\w*)/) {
	$last_subrt = $4;
	print;
    }
    elsif (/^\s*end\s+(subroutine|function)\s*$/) { 
	chomp();
	print $_." ".$last_subrt."\n";
    }
    else { s/^\s*private\s*$/!private/; s/^\s*private([ :])/!private$1/; print; }
}' %s > tmp.out; %s %s tmp.out | grep -v '^#' >  %s""" % (infile, ' '.join(cpp), cpp_opt, outfile))    
        else:
            copy_file(infile, outfile, update=True)

        return outfile
    
    return func


def F90WrapperBuilder(modname, all_sources, wrap_sources, dep_type_maps=[], kindlines=[], short_names={}, initlines={}, 
                      filtertypes=None, prefix='',callback_routines={}):
    """Build a Fortran 90 wrapper for the given F95 source files
    that is suitable for use with f2py. Derived types are wrapped to 
    give access to methods and instance data."""

    def func(extension, build_dir):
        in_files = ['%s/../%s' % (build_dir, f) for f in all_sources]
        wrap_files = ['%s/../%s' % (build_dir, f) for f in wrap_sources]

        if newer_group(wrap_files, '%s.f90doc' % modname):
            programs, modules, functs, subts = f90doc.read_files(in_files)
            cPickle.dump((programs, modules, functs, subts), open('%s.f90doc' % modname, 'w'))
        else:
            (programs, modules, functs, subts) = cPickle.load(open('%s.f90doc' % modname))

        for mod, name in modules:
            for n in [t.name for t in mod.types]:
                type_map[n.lower()] = mod.name
            
        for item in dep_type_maps:
            if hasattr(item, '__getitem__') and hasattr(item, 'keys'): # dictionary
                type_map.update(item)
            else: # assume it's a string
                type_map.update(cPickle.load(open('%s.type' % item)))

        res = []
        fortran_spec = {}
        if os.path.exists('%s.spec' % modname):
            fortran_spec = cPickle.load(open('%s.spec' % modname))

        wrap_modules = []
        for file in wrap_sources:

            for mod, name in modules:
                if os.path.basename(name) == os.path.basename(file):
                    break
            else:
                raise ValueError("Can't find module %s" % file)

            wrap_mod_name = mod.name.lower()[:-7]
            wrap_modules.append(wrap_mod_name)

            wrapper = '%s/%s_%s_wrap.f90' % (build_dir, modname, wrap_mod_name)

            if not newer(name, wrapper):
                res.append(wrapper)
                continue

            tmpf = StringIO.StringIO()
            new_spec = f2py_wrapper_gen.wrap_mod(mod, type_map, tmpf, kindlines=kindlines, initlines=initlines,
                                                 filtertypes=filtertypes, prefix=prefix, callback_routines=callback_routines)

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
        cPickle.dump(fortran_spec, open('%s.spec' % modname, 'w'))

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
    f.write("""include Makefile.inc
include Makefile.rules
ifeq (${QUIP_ARCH},)
  include Makefile.arch
else
  include Makefile.${QUIP_ARCH}
endif""")
    f.close()

    # Dump make database to file
    os.system('make -f Makefile.quippy -I %s -I %s -I %s -p > make.dump' %
              (quip_root,
               os.path.join(quip_root, 'Makefiles'),
               os.path.join(quip_root, 'build.%s' % quip_arch)))

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


def find_sources(makefile, quip_root):

    source_dirs  = []
    all_sources  = []
    wrap_sources = []
    wrap_types   = []

    libatoms_dir   = os.path.join(quip_root, 'libAtoms/')
    makefile_libatoms   = parse_makefile(os.path.join(libatoms_dir,'Makefile'))
    libatoms_sources = [os.path.join(libatoms_dir,f) for f in
                        (expand_addsuffix(makefile_libatoms['LIBATOMS_F77_SOURCES']).split() +
                         expand_addsuffix(makefile_libatoms['LIBATOMS_F95_SOURCES']).split() +
                         expand_addsuffix(makefile_libatoms['LIBATOMS_C_SOURCES']).split())]
    all_sources += libatoms_sources
    wrap_sources += ['System.f95', 'MPI_context.f95', 'Units.f95', 'linearalgebra.f95',
                     'Dictionary.f95', 'Table.f95', 'PeriodicTable.f95', 'Atoms.f95', 'DynamicalSystem.f95',
                     'clusters.f95','Structures.f95', 'CInOutput.f95', 'Topology.f95']
    wrap_types += ['inoutput', 'mpi_context', 'dictionary', 'table', 'atoms', 'connection', 'dynamicalsystem', 'cinoutput']
    source_dirs.append(libatoms_dir)

    if 'HAVE_GP' in makefile and int(makefile['HAVE_GP']) == 1:
        gp_dir = os.path.join(quip_root, 'gp')
        gp_makefile = parse_makefile(os.path.join(gp_dir, 'Makefile'))
        gp_sources = [os.path.join(gp_dir, f) for f in expand_addsuffix(gp_makefile['GP_F95_SOURCES']).split()]        
        all_sources += gp_sources
        # wrap_sources += ... # list of Fortran files to wrap
        # wrap_types += .... # list of types to wrap
        source_dirs.append(gp_dir)


    quip_core_dir = os.path.join(quip_root, 'QUIP_Core/')
    makefile_quip_core = parse_makefile(os.path.join(quip_core_dir, 'Makefile'))
    quip_core_sources = [os.path.join(quip_core_dir,f) for f in
                         (expand_addsuffix(makefile_quip_core['TB_F77_SOURCES']).split() +
                          expand_addsuffix(makefile_quip_core['ALL_F95_FILES']).split()  + 
                          expand_addsuffix(makefile_quip_core['POT_F95_SOURCES']).split())]
    all_sources += quip_core_sources
    wrap_sources += ['Potential.f95']
    wrap_types += ['potential']
    source_dirs.append(quip_core_dir)

    if (not 'QUIPPY_NO_TOOLS' in makefile or
        ('QUIPPY_NO_TOOLS' in makefile and not int(makefile['QUIPPY_NO_TOOLS']))):
        quip_utils_dir = os.path.join(quip_root, 'QUIP_Utils/')
        makefile_quip_utils = parse_makefile(os.path.join(quip_utils_dir, 'Makefile'))
        quip_utils_sources = [os.path.join(quip_utils_dir,f)+'.f95' for f in makefile_quip_utils['F95_FILES'].split()]
        all_sources += quip_utils_sources
        wrap_sources += ['elasticity.f95']
        source_dirs.append(quip_utils_dir)

    if (not 'QUIPPY_NO_CRACK' in makefile or
        ('QUIPPY_NO_CRACK' in makefile and not int(makefile['QUIPPY_NO_CRACK']))):
        crack_dir = os.path.join(quip_root, 'QUIP_Programs/')
        crack_sources = [os.path.join(crack_dir,f) for f in ('crackparams.f95', 'cracktools.f95')]
        all_sources += crack_sources
        wrap_sources += ['cracktools.f95', 'crackparams.f95']
        wrap_types += ['crackparams']
        source_dirs.append(crack_dir)

    return source_dirs, all_sources, wrap_sources, wrap_types


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
libraries = [s.startswith('-l') and s[2:] or s  for s in fields
             if s.startswith('-l') or s == '-Bstatic' or s == '-Bdynamic' or s.startswith('-Wl') ]
library_dirs  = [s[2:] for s in fields if s[:2] == '-L']

# everything else in SYSLIBS is an extra link argument
extra_link_args = [s for s in fields if not s[2:] in libraries and not s in libraries and not s[2:] in library_dirs and not s[2:] in include_dirs]

# Preprocessor macros
macros = [('HAVE_QUIPPY',None), ('SVN_VERSION',r'\"%s\"' % os.popen('svnversion -n .').read())]
for defn in makefile['DEFINES'].split():
    if defn[:2] == '-D':
        if '=' in defn:
            n, v = defn[2:].split('=')
            macros.append((n,v))
        else:
            macros.append((defn[2:], None))
    elif defn[:2] == '-U':
        macros.append(defn[2:])

macros.extend([
    ('NUMPY_PTR_SIZE', numpy.dtype('O').itemsize),
    ('NUMPY_INTEGER',  numpy.dtype('int32').num),
    ('NUMPY_REAL_DP',  numpy.dtype('d').num),
    ('NUMPY_LOGICAL',  numpy.dtype('int32').num),
    ('NUMPY_COMPLEX',  numpy.dtype('complex').num),
    ('NUMPY_CHAR',     numpy.dtype('S').num)
    ])

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

if makefile_test('QUIPPY_DEBUG'):
    os.environ['FOPT'] = '-O0'
    os.environ['FARCH'] = ''

if 'QUIPPY_OPT' in makefile:
    default_options['config_fc']['opt'] = makefile['QUIPPY_OPT'].split()

# Install options
if 'QUIPPY_INSTALL_OPTS' in makefile:
    install_opts = makefile['QUIPPY_INSTALL_OPTS'].split()
    default_options['install'] = {}
    for opt in install_opts:
        n, v = opt.split('=')
        n = n[2:] # remove --
        default_options['install'][n] = v

# Find Fortran source code files
source_dirs, all_sources, wrap_sources, wrap_types = find_sources(makefile, quip_root)
include_dirs.extend(source_dirs)

# arraydata extension module
f2py_info = get_info('f2py')
arraydata_ext = Extension(name='quippy.arraydata', 
                          sources=['arraydatamodule.c'] + f2py_info['sources'],
                          include_dirs=f2py_info['include_dirs'])

# underlying quippy library, libquippy.a
quippy_lib = ('quippy', {
    'sources': [ SourceImporter(f, macros, source_dirs, cpp) for f in all_sources ],
    'include_dirs':  include_dirs,
    'macros': macros,
    })

# _quippy extension module
quippy_ext = Extension(name='quippy._quippy',
                       sources=[ F90WrapperBuilder('quippy',
                                                   [os.path.basename(f)[:-4]+'.f90' for f in all_sources if f.endswith('.f95')],
                                                   [os.path.basename(f)[:-4]+'.f90' for f in wrap_sources],
                                                   dep_type_maps=[{'c_ptr': 'iso_c_binding',
                                                                   'dictionary_t':'FoX_sax'}], 
                                                   kindlines=['use system_module, only: dp, qp'],
                                                   short_names={'dynamicalsystem':'ds',
                                                                'metapotential': 'metapot'},
                                                   initlines={'atoms': ('atoms_module', ('call atoms_repoint(%(PTR)s)',
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
                       libraries=[quippy_lib] + libraries,
                       define_macros= macros,
                       extra_link_args=extra_link_args)

exts = [arraydata_ext, quippy_ext]

# Finally, call setup() to run command
setup(name='quippy',
      packages = ['quippy'],
      ext_modules = exts,
      data_files = [('quippy',['quippy.spec'])],
      scripts=glob.glob('scripts/*.py'),      
      cmdclass = {'clean': clean, 'test': test, 'build_ext': build_ext},
      version=os.popen('svnversion -n .').read(),
      description='Python bindings to QUIP code',
      author='James Kermode',
      author_email='james.kermode@kcl.ac.uk',
      url='http://www.jrkermode.co.uk/quippy',
      options=default_options)

#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
from numpy.distutils.system_info import get_info
from numpy.distutils.misc_util import dict_append
from numpy.distutils.ccompiler import new_compiler, gen_preprocess_options 
from distutils.file_util import copy_file
from distutils.dir_util import remove_tree
from distutils.dep_util import newer, newer_group
from distutils.util import get_platform
from numpy import get_include
from distutils.command.clean import clean as _clean
import sys, os, cPickle, glob, stat, subprocess
import f90doc

print sys.argv

major, minor = sys.version_info[0:2]
if (major, minor) < (2, 4):
    sys.stderr.write('Python 2.4 or later is needed to use this package\n')
    sys.exit(1)

from numpy.f2py.rules import f2py_version
if not f2py_version.endswith('patched_JRK'):
    sys.stderr.write('\nnumpy.f2py.rules must be patched to use this package; see README for details\n')
    sys.exit(1)

class clean(_clean):
    def run(self):
        _clean.run(self)

        for file in glob.glob('*.spec') + glob.glob('*.type') + ['sizeof_void_ptr', 'sizeof_void_ptr.o']:
            os.remove(file)

        remove_tree(build_base)


def SourceImporter(infile, defines, include_dirs, cpp):
    """Import source code from infile and copy to build_dir/package,
       passing through filter if it's an F95 source file. The filter
       converts all private variables to public and runs the cpp
       preprocessor."""

    def func(extension, build_dir):

        outfile = os.path.join(build_dir, os.path.split(infile)[1])

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
    else { s/^\s*private/!private/; print; }
}' %s > tmp.out; %s %s tmp.out -o %s""" % (infile, ' '.join(cpp), cpp_opt, outfile))    
#}' %s | %s %s - | perl -ne 'print if !/^$/' > %s""" % (infile, ' '.join(cpp), cpp_opt, outfile))
        else:
            #if newer(infile, outfile):
            #    os.system("cat %s | cpp %s > %s" % (infile, cpp_opt, outfile))
            copy_file(infile, outfile, update=True)

        return outfile
    
    return func


def F90WrapperBuilder(modname, f95_sources, cpp, dep_type_maps=[], donothing=False, kindlines=[]):
    """Build a Fortran 90 wrapper for the given F95 source files
    that is suitable for use with f2py. Derived types are wrapped to 
    give access to methods and instance data."""

    def func(extension, build_dir):

        wrapper = '%s/%swrap.f90' % (build_dir, modname)

        if donothing and os.path.exists(wrapper): return wrapper

        in_files = ['%s/../%s' % (build_dir, f) for f in f95_sources]

        if os.path.exists('%s.spec' % modname) and not newer_group(in_files + [f90doc.__file__], wrapper): return wrapper

        programs, modules, functs, subts = f90doc.read_files(in_files)

        type_map = {}
        for item in dep_type_maps:
            if hasattr(item, '__getitem__') and hasattr(item, 'keys'): # dictionary
                type_map.update(item)
            else: # assume it's a string
                type_map.update(cPickle.load(open('%s.type' % item)))

        for mod, name in modules:
            for n in [t.name for t in mod.types]:
                type_map[n.lower()] = mod.name
                
        cPickle.dump(type_map, open('%s.type' % modname, 'w'))

        fortran_spec = {}
        cpp_opt = gen_preprocess_options(macros, [])
        wrapperf = open(wrapper, 'w')
        cpp_process = subprocess.Popen(cpp + cpp_opt + ['-'], stdin=subprocess.PIPE, stdout=wrapperf)
        for mod, name in modules:
            mod.f2py(type_map, fortran_spec, cpp_process.stdin,
                     kindlines=kindlines)
        cpp_process.communicate()

        cPickle.dump(fortran_spec, open('%s.spec' % modname, 'w'))

        return wrapper

    return func

# Compile simple C program to find sizeof(void *) on this arch
print 'Calculating sizeof(void *)...', 
if not os.path.exists('./sizeof_void_ptr'):
    cc = new_compiler()
    cc.link_executable(cc.compile(['sizeof_void_ptr.c']),'sizeof_void_ptr')
sizeof_void_ptr = int(os.popen('./sizeof_void_ptr').read())
print sizeof_void_ptr, 'bytes.'

# Bit of a hack: we have to add directory which will contain .mod files
build_base = 'build'
if 'build' in sys.argv:
    for i, arg in enumerate(sys.argv[sys.argv.index('build'):]):
        if not arg.startswith('-'): break
        if arg == '-b':
            build_base = sys.argv[i+1]
            break
        if arg.startswith('--build-base'):
            build_base = arg.split('=')[1]
            break
print 'build_base directory is %s' % build_base
plat_specifier = ".%s-%s" % (get_platform(), sys.version[0:3])
mod_dir = os.path.join(build_base, 'temp'+plat_specifier)

type_map = {}

include_dirs = [os.path.expanduser(s[2:]) for s in sys.argv if s.startswith('-I')]
library_dirs = [os.path.expanduser(s[2:]) for s in sys.argv if s.startswith('-L')]
libraries = [s[2:] for s in sys.argv if s.startswith('-l')]
sys.argv = [ s for s in sys.argv if not s.startswith('-I') and not s.startswith('-L') and not s.startswith('-l')]

argfilt = filter(lambda s: s.startswith('-D'), sys.argv)
macros = []
for arg in argfilt:
    try:
        k, v = arg[2:].split('=')
        v = v.replace('"',r'\"')
    except ValueError:
        k, v = arg[2:], None
    macros.append((k,v))
    del sys.argv[sys.argv.index(arg)]

print 'macros = ', macros

do_quippy_extension = True

argfilt = filter(lambda s: s.startswith('--libatoms-dir'), sys.argv)
if argfilt:
    libatoms_dir = argfilt[0].split('=')[1]
    del sys.argv[sys.argv.index(argfilt[0])]
else:
    do_quippy_extension = False

argfilt = filter(lambda s: s.startswith('--libatoms-sources'), sys.argv)
if argfilt:
    libatoms_sources = argfilt[0].split('=')[1].split(':')
    del sys.argv[sys.argv.index(argfilt[0])]
else:
    do_quippy_extension = False

argfilt = filter(lambda s: s.startswith('--quip-core-dir'), sys.argv)
if argfilt:
    quip_core_dir = argfilt[0].split('=')[1]
    del sys.argv[sys.argv.index(argfilt[0])]
else:
    do_quippy_extension = False

argfilt = filter(lambda s: s.startswith('--quip-core-sources'), sys.argv)
if argfilt:
    quip_core_sources = argfilt[0].split('=')[1].split(':')
    del sys.argv[sys.argv.index(argfilt[0])]
else:
    do_quippy_extension = False

cpp = 'cpp'
argfilt = filter(lambda s: s.startswith('--cpp'), sys.argv)
if argfilt:
    cpp = argfilt[0].split('=')[1].split()
    del sys.argv[sys.argv.index(argfilt[0])]

do_atomeye = False
argfilt = filter(lambda s: s.startswith('--atomeye-dir'), sys.argv)
if argfilt:
    do_atomeye = True
    atomeye_dir = argfilt[0].split('=')[1]
    del sys.argv[sys.argv.index(argfilt[0])]

if do_quippy_extension:
    macros += [('SIZEOF_VOID_PTR', sizeof_void_ptr), ('HAVE_QUIPPY',None)]
    arraydata_ext = Extension(name='quippy.arraydata', 
                              sources=['arraydatamodule.c'],
                              include_dirs=[get_include()])


    libatoms_files = [ os.path.join(libatoms_dir, f) for f in libatoms_sources ]

    libatoms_lib = ('atoms', {
            'sources': [ SourceImporter(f, macros, [libatoms_dir], cpp) for f in libatoms_files ],
            'include_dirs': include_dirs + [libatoms_dir],
            'macros': macros
            })

    quip_core_files = [os.path.join(quip_core_dir, f) for f in quip_core_sources]

    quip_core_lib = ('quip_core', {
        'sources': [ SourceImporter(f, macros, [quip_core_dir], cpp) for f in quip_core_files ],
        'include_dirs': include_dirs + [quip_core_dir],
        'macros': macros
        })


    data_files = ['quippy.spec']

    ext_args = {'name': 'quippy._quippy',
                'sources': [ F90WrapperBuilder('quippy', filter(lambda f: f.endswith('.f95'), libatoms_sources + quip_core_sources),
                                               cpp, dep_type_maps=[{'c_ptr': 'iso_c_binding',
                                                                    'dictionary_t':'FoX_sax'}], 
                                               donothing=False,
                                               kindlines=['use system_module, only: dp',
                                                          'use iso_c_binding, only: C_SIZE_T'])
                           ],
                'library_dirs': library_dirs,
                'include_dirs': include_dirs + [mod_dir],
                'libraries':  ['quip_core', 'atoms'] + libraries,
                'define_macros': macros
                }

    lapack_opt = get_info('lapack_opt')
    if not lapack_opt:
        print 'No lapack_opt resources found in system'
        sys.exit(1)
    dict_append(ext_args,**lapack_opt)

    quippy_ext = Extension(**ext_args)

    build_libraries = [libatoms_lib, quip_core_lib]
    exts = [arraydata_ext, quippy_ext]

    if do_atomeye:
        ext_args = {'name': 'quippy._atomeye',
                    'sources': ['atomeyemodule.c'],
                    'library_dirs':  library_dirs + [os.path.join(atomeye_dir, 'lib'),'/usr/X11/lib'],
                    'libraries': ['AtomEye', 'AX', 'png', 'z', 'jpeg', 'Atoms', 'VecMat3', 'VecMat', 'IO', 'Scalar', 'Timer',
                                  'history', 'ncurses', 'm', 'Xpm', 'Xext', 'X11', 'readline'] + libraries,
                    'include_dirs': [os.path.join(atomeye_dir,'include')] + [libatoms_dir],
                    'depends': [os.path.join(atomeye_dir, 'lib/libAtomEye.a')],
                    'define_macros': macros
                    }

        dict_append(ext_args, **lapack_opt)
        atomeye_ext = Extension(**ext_args)
        exts.append(atomeye_ext)
else:
    build_libraries = []
    exts = []
    data_files = []


setup(name='quippy',
      libraries=build_libraries,
      packages = ['quippy'],
      ext_modules = exts,
      data_files = [('quippy',data_files)],
      cmdclass = {'clean':clean},
      version=os.popen('svnversion -n .').read(),
      description='Python binding to QUIP code',
      author='James Kermode',
      author_email='james.kermode@kcl.ac.uk',
      url='http://www.libatoms.org')

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

"""Overridden versions of distutils `clean` and `build_ext` commands,
and a new `test` command to find and run unittests."""

import sys, os, glob, re, string

from distutils.command.clean import clean as _clean
from numpy.distutils.command.build_ext import build_ext as _build_ext
from numpy.distutils.core import Command

from distutils import log
from distutils.file_util import copy_file
from distutils.dir_util import remove_tree
from distutils.dep_util import newer, newer_group
from numpy.distutils.misc_util import filter_sources, has_f_sources, \
     has_cxx_sources, get_ext_source_files, \
     get_numpy_include_dirs, is_sequence, get_build_architecture
from numpy.distutils.ccompiler import new_compiler
from unittest import TestLoader, TextTestRunner

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

class clean(_clean):
    def initialize_options(self):
        _clean.initialize_options(self)
        self.build_src = None

    def finalize_options(self):
        _clean.finalize_options(self)
        self.set_undefined_options('build_src', ('build_src', 'build_src'))

    def run(self):
        # remove the build/temp.<plat> directory (unless it's already
        # gone)
        if os.path.exists(self.build_temp):
            remove_tree(self.build_temp, dry_run=self.dry_run)
        else:
            log.debug("'%s' does not exist -- can't clean it",
                      self.build_temp)

        if self.all:
            # remove build directories
            for directory in (self.build_lib,
                              self.bdist_base,
                              self.build_scripts,
                              self.build_src):
                if os.path.exists(directory):
                    remove_tree(directory, dry_run=self.dry_run)
                else:
                    log.warn("'%s' does not exist -- can't clean it",
                             directory)

        if not self.dry_run and os.path.exists(self.build_base):
            remove_tree(self.build_base)
            log.info("removing '%s'", self.build_base)

class build_ext(_build_ext):
    """Subclass of numpy.distutils.command.build_ext.build_ext with two new features:

    (1) We add module_build_dir of included libraries to include_dirs for all extensions.
    (2) New build_ext dep_libs options which forces a re-link of the Extension .so library
        if any of the libraries with which it's being linked are newer than it."""

    user_options = _build_ext.user_options + [
        ('dep-libs', None, 're-link shared object if any libraries are newer than it')
        ]

    def initialize_options(self):
        _build_ext.initialize_options(self)
        self.dep_libs = None

    def finalize_options(self):
        _build_ext.finalize_options(self)
    
    def run(self):
        if not self.extensions:
            return

        # Make sure that extension sources are complete.
        self.run_command('build_src')

        if self.distribution.has_c_libraries():
            self.run_command('build_clib')
            build_clib = self.get_finalized_command('build_clib')
            self.library_dirs.append(build_clib.build_clib)
            ## JRK addition
            for ext in self.extensions:
                ext.include_dirs.append(build_clib.build_clib)
            ### end addition
        else:
            build_clib = None

        # Not including C libraries to the list of
        # extension libraries automatically to prevent
        # bogus linking commands. Extensions must
        # explicitly specify the C libraries that they use.

        from distutils.ccompiler import new_compiler
        from numpy.distutils.fcompiler import new_fcompiler

        compiler_type = self.compiler
        # Initialize C compiler:
        self.compiler = new_compiler(compiler=compiler_type,
                                     verbose=self.verbose,
                                     dry_run=self.dry_run,
                                     force=self.force)
        self.compiler.customize(self.distribution)
        self.compiler.customize_cmd(self)
        self.compiler.show_customization()

        # Create mapping of libraries built by build_clib:
        clibs = {}
        if build_clib is not None:
            for libname,build_info in build_clib.libraries or []:
                if libname in clibs and clibs[libname] != build_info:
                    log.warn('library %r defined more than once,'\
                             ' overwriting build_info\n%s... \nwith\n%s...' \
                             % (libname, `clibs[libname]`[:300], `build_info`[:300]))
                clibs[libname] = build_info
        # .. and distribution libraries:
        for libname,build_info in self.distribution.libraries or []:
            if libname in clibs:
                # build_clib libraries have a precedence before distribution ones
                continue
            clibs[libname] = build_info

        # Determine if C++/Fortran 77/Fortran 90 compilers are needed.
        # Update extension libraries, library_dirs, and macros.
        all_languages = set()
        for ext in self.extensions:
            ext_languages = set()
            c_libs = []
            c_lib_dirs = []
            macros = []
            for libname in ext.libraries:
                if libname in clibs:
                    binfo = clibs[libname]
                    c_libs += binfo.get('libraries',[])
                    c_lib_dirs += binfo.get('library_dirs',[])
                    for m in binfo.get('macros',[]):
                        if m not in macros:
                            macros.append(m)

                for l in clibs.get(libname,{}).get('source_languages',[]):
                    ext_languages.add(l)
            if c_libs:
                new_c_libs = ext.libraries + c_libs
                log.info('updating extension %r libraries from %r to %r'
                         % (ext.name, ext.libraries, new_c_libs))
                ext.libraries = new_c_libs
                ext.library_dirs = ext.library_dirs + c_lib_dirs
            if macros:
                log.info('extending extension %r defined_macros with %r'
                         % (ext.name, macros))
                ext.define_macros = ext.define_macros + macros

            # determine extension languages
            if has_f_sources(ext.sources):
                ext_languages.add('f77')
            if has_cxx_sources(ext.sources):
                ext_languages.add('c++')
            l = ext.language or self.compiler.detect_language(ext.sources)
            if l:
                ext_languages.add(l)
            # reset language attribute for choosing proper linker
            if 'c++' in ext_languages:
                ext_language = 'c++'
            elif 'f90' in ext_languages:
                ext_language = 'f90'
            elif 'f77' in ext_languages:
                ext_language = 'f77'
            else:
                ext_language = 'c' # default
            if l and l != ext_language and ext.language:
                log.warn('resetting extension %r language from %r to %r.' %
                         (ext.name,l,ext_language))
            ext.language = ext_language
            # global language
            all_languages.update(ext_languages)

        need_f90_compiler = 'f90' in all_languages
        need_f77_compiler = 'f77' in all_languages
        need_cxx_compiler = 'c++' in all_languages

        # Initialize C++ compiler:
        if need_cxx_compiler:
            self._cxx_compiler = new_compiler(compiler=compiler_type,
                                             verbose=self.verbose,
                                             dry_run=self.dry_run,
                                             force=self.force)
            compiler = self._cxx_compiler
            compiler.customize(self.distribution,need_cxx=need_cxx_compiler)
            compiler.customize_cmd(self)
            compiler.show_customization()
            self._cxx_compiler = compiler.cxx_compiler()
        else:
            self._cxx_compiler = None

        # Initialize Fortran 77 compiler:
        if need_f77_compiler:
            ctype = self.fcompiler
            self._f77_compiler = new_fcompiler(compiler=self.fcompiler,
                                               verbose=self.verbose,
                                               dry_run=self.dry_run,
                                               force=self.force,
                                               requiref90=False,
                                               c_compiler=self.compiler)
            fcompiler = self._f77_compiler
            if fcompiler:
                ctype = fcompiler.compiler_type
                fcompiler.customize(self.distribution)
            if fcompiler and fcompiler.get_version():
                fcompiler.customize_cmd(self)
                fcompiler.show_customization()
            else:
                self.warn('f77_compiler=%s is not available.' %
                          (ctype))
                self._f77_compiler = None
        else:
            self._f77_compiler = None

        # Initialize Fortran 90 compiler:
        if need_f90_compiler:
            ctype = self.fcompiler
            self._f90_compiler = new_fcompiler(compiler=self.fcompiler,
                                               verbose=self.verbose,
                                               dry_run=self.dry_run,
                                               force=self.force,
                                               requiref90=True,
                                               c_compiler = self.compiler)
            fcompiler = self._f90_compiler
            if fcompiler:
                ctype = fcompiler.compiler_type
                fcompiler.customize(self.distribution)
            if fcompiler and fcompiler.get_version():
                fcompiler.customize_cmd(self)
                fcompiler.show_customization()
            else:
                self.warn('f90_compiler=%s is not available.' %
                          (ctype))
                self._f90_compiler = None
        else:
            self._f90_compiler = None

        # Build extensions
        self.build_extensions()

        # Make sure that scons based extensions are complete.
        #if self.inplace:
        #    cmd = self.reinitialize_command('scons')
        #    cmd.inplace = 1
        #self.run_command('scons')


    def build_extension(self, ext):
        sources = ext.sources
        if sources is None or not is_sequence(sources):
            raise DistutilsSetupError(
                ("in 'ext_modules' option (extension '%s'), " +
                 "'sources' must be present and must be " +
                 "a list of source filenames") % ext.name)
        sources = list(sources)

        if not sources:
            return

        fullname = self.get_ext_fullname(ext.name)
        if self.inplace:
            modpath = fullname.split('.')
            package = '.'.join(modpath[0:-1])
            base = modpath[-1]
            build_py = self.get_finalized_command('build_py')
            package_dir = build_py.get_package_dir(package)
            ext_filename = os.path.join(package_dir,
                                        self.get_ext_filename(base))
        else:
            ext_filename = os.path.join(self.build_lib,
                                        self.get_ext_filename(fullname))
        depends = sources + ext.depends

        log.info("building '%s' extension", ext.name)

        extra_args = ext.extra_compile_args or []
        macros = ext.define_macros[:]
        for undef in ext.undef_macros:
            macros.append((undef,))

        c_sources, cxx_sources, f_sources, fmodule_sources = \
                   filter_sources(ext.sources)

        if self.compiler.compiler_type=='msvc':
            if cxx_sources:
                # Needed to compile kiva.agg._agg extension.
                extra_args.append('/Zm1000')
            # this hack works around the msvc compiler attributes
            # problem, msvc uses its own convention :(
            c_sources += cxx_sources
            cxx_sources = []

        # Set Fortran/C++ compilers for compilation and linking.
        if ext.language=='f90':
            fcompiler = self._f90_compiler
        elif ext.language=='f77':
            fcompiler = self._f77_compiler
        else: # in case ext.language is c++, for instance
            fcompiler = self._f90_compiler or self._f77_compiler
        cxx_compiler = self._cxx_compiler

        # check for the availability of required compilers
        if cxx_sources and cxx_compiler is None:
            raise DistutilsError, "extension %r has C++ sources" \
                  "but no C++ compiler found" % (ext.name)
        if (f_sources or fmodule_sources) and fcompiler is None:
            raise DistutilsError, "extension %r has Fortran sources " \
                  "but no Fortran compiler found" % (ext.name)
        if ext.language in ['f77','f90'] and fcompiler is None:
            self.warn("extension %r has Fortran libraries " \
                  "but no Fortran linker found, using default linker" % (ext.name))
        if ext.language=='c++' and cxx_compiler is None:
            self.warn("extension %r has C++ libraries " \
                  "but no C++ linker found, using default linker" % (ext.name))

        kws = {'depends':ext.depends}
        output_dir = self.build_temp

        include_dirs = ext.include_dirs + get_numpy_include_dirs()

        c_objects = []
        if c_sources:
            log.info("compiling C sources")
            c_objects = self.compiler.compile(c_sources,
                                              output_dir=output_dir,
                                              macros=macros,
                                              include_dirs=include_dirs,
                                              debug=self.debug,
                                              extra_postargs=extra_args,
                                              **kws)

        if cxx_sources:
            log.info("compiling C++ sources")
            c_objects += cxx_compiler.compile(cxx_sources,
                                              output_dir=output_dir,
                                              macros=macros,
                                              include_dirs=include_dirs,
                                              debug=self.debug,
                                              extra_postargs=extra_args,
                                              **kws)

        extra_postargs = []
        f_objects = []
        if fmodule_sources:
            log.info("compiling Fortran 90 module sources")
            module_dirs = ext.module_dirs[:]
            module_build_dir = os.path.join(
                self.build_temp,os.path.dirname(
                    self.get_ext_filename(fullname)))

            self.mkpath(module_build_dir)
            if fcompiler.module_dir_switch is None:
                existing_modules = glob('*.mod')
            extra_postargs += fcompiler.module_options(
                module_dirs,module_build_dir)
            f_objects += fcompiler.compile(fmodule_sources,
                                           output_dir=self.build_temp,
                                           macros=macros,
                                           include_dirs=include_dirs,
                                           debug=self.debug,
                                           extra_postargs=extra_postargs,
                                           depends=ext.depends)

            if fcompiler.module_dir_switch is None:
                for f in glob('*.mod'):
                    if f in existing_modules:
                        continue
                    t = os.path.join(module_build_dir, f)
                    if os.path.abspath(f)==os.path.abspath(t):
                        continue
                    if os.path.isfile(t):
                        os.remove(t)
                    try:
                        self.move_file(f, module_build_dir)
                    except DistutilsFileError:
                        log.warn('failed to move %r to %r' %
                                 (f, module_build_dir))
        if f_sources:
            log.info("compiling Fortran sources")
            f_objects += fcompiler.compile(f_sources,
                                           output_dir=self.build_temp,
                                           macros=macros,
                                           include_dirs=include_dirs,
                                           debug=self.debug,
                                           extra_postargs=extra_postargs,
                                           depends=ext.depends)

        objects = c_objects + f_objects

        if ext.extra_objects:
            objects.extend(ext.extra_objects)
        extra_args = ext.extra_link_args or []
        libraries = self.get_libraries(ext)[:]
        library_dirs = ext.library_dirs[:]

        linker = self.compiler.link_shared_object
        # Always use system linker when using MSVC compiler.
        if self.compiler.compiler_type=='msvc':
            # expand libraries with fcompiler libraries as we are
            # not using fcompiler linker
            self._libs_with_msvc_and_fortran(fcompiler, libraries, library_dirs)

        elif ext.language in ['f77','f90'] and fcompiler is not None:
            linker = fcompiler.link_shared_object
        if ext.language=='c++' and cxx_compiler is not None:
            linker = cxx_compiler.link_shared_object

        if sys.version[:3]>='2.3':
            kws = {'target_lang':ext.language}
        else:
            kws = {}

        if self.dep_libs is not None:
            dep_libs = [ self.compiler.find_library_file(library_dirs+self.library_dirs, lib) for lib in libraries ]
            dep_libs = [ libfile for libfile in dep_libs if libfile is not None ]

            if newer_group(dep_libs, ext_filename):
                self.compiler.force = 1
                fcompiler.force = 1

        linker(objects, ext_filename,
               libraries=libraries,
               library_dirs=library_dirs,
               runtime_library_dirs=ext.runtime_library_dirs,
               extra_postargs=extra_args,
               export_symbols=self.get_export_symbols(ext),
               debug=self.debug,
               build_temp=self.build_temp,**kws)

        if self.dep_libs is not None:
            self.compiler.force = 0
            fcompiler.force = 0

        


class test(Command):
    # Original version of this class posted
    # by Berthold H=F6llmann to distutils-sig@python.org
    description = "test the distribution prior to install"

    user_options = [
        ('test-dir=', None,
         "directory that contains the test definitions"),
        ('test-prefix=', None,
         "prefix to the testcase filename"),
        ('test-suffixes=', None,
         "a list of suffixes used to generate names of the testcases"),
        ('verbosity=', None,
         "verbosity for unit tests"),
        ("test=", None,
         "single test to run"),
        ("debug", None,
          "start debug ipython shell after running test")
        ]

    def initialize_options(self):
        self.build_base = 'build'
        # these are decided only after 'build_base' has its final value
        # (unless overridden by the user or client)
        self.test_dir = 'tests'
        self.test_prefix = 'test_'
        self.test_suffixes = None
        self.verbosity = None
        self.test = None
        self.debug = None

    def finalize_options(self):
        import os
        if self.test_suffixes is None:
            self.test_suffixes = []
            pref_len = len(self.test_prefix)
            for file in os.listdir(self.test_dir):
                if (file[-3:] == ".py" and
                    file[:pref_len]==self.test_prefix):
                    self.test_suffixes.append(file[pref_len:-3])

        if self.verbosity is None:
            self.verbosity = 1
        else:
            self.verbosity = int(self.verbosity)

        build = self.get_finalized_command('build')
        self.build_purelib = build.build_purelib
        self.build_platlib = build.build_platlib

    def run(self):
        import sys
        # Invoke the 'build' command to "build" pure Python modules
        # (ie. copy 'em into the build tree)
        self.run_command('build')

        # remember old sys.path to restore it afterwards
        old_path = sys.path[:]

        # extend sys.path
        sys.path.insert(0, self.build_purelib)
        sys.path.insert(0, self.build_platlib)
        sys.path.insert(0, self.test_dir)

        # run tests
        if self.test is not None:
            tests = TestLoader().loadTestsFromNames([self.test])
        else:
            tests = TestLoader().loadTestsFromNames([self.test_prefix+case for case in self.test_suffixes])


        if self.debug is not None:
            import sys
            from IPython.core import ultratb
            sys.excepthook = ultratb.FormattedTB(mode='Verbose',
                                                 color_scheme='LightBG', call_pdb=1)
            tests.debug()
        else:
           print 'Running tests with verbosity %d' % self.verbosity
           runner = TextTestRunner(verbosity=self.verbosity)
           result = runner.run(tests)

           # Exit script with exitcode 1 if any tests failed
           if result.failures or result.errors:
              sys.exit(1)

        # restore sys.path
        sys.path = old_path[:]

        


class interact(Command):
    description = "interact with the distribution prior to install"

    user_options = [('execute=', None, "code to execute on startup"),
                    ('execfile=', None, "file to execute on startup")]

    def initialize_options(self):
        self.build_base = 'build'
        self.execute = None
        self.execfile = None

    def finalize_options(self):
        build = self.get_finalized_command('build')
        self.build_purelib = build.build_purelib
        self.build_platlib = build.build_platlib

    def run(self):
        import sys
        # Invoke the 'build' command to "build" pure Python modules
        # (ie. copy 'em into the build tree)
        self.run_command('build')

        # remember old sys.path to restore it afterwards
        old_path = sys.path[:]

        # start embedded ipython shell
        from IPython import embed
        #ipshell = IPShellEmbed(['-pdb', '-pi1','interact In <\\#>: ','-pi2','interact    .\\D.: ',
        #                        '-po','interact Out<\\#>: ','-nosep'])

        # extend sys.path
        sys.path.insert(0, self.build_purelib)
        sys.path.insert(0, self.build_platlib)
        if os.getcwd() in sys.path: sys.path.remove(os.getcwd())
        if '' in sys.path: sys.path.remove('')
        
        from quippy import *
        if self.execute is not None:
           exec(self.execute)
        if self.execfile is not None:
           execfile(self.execfile)
        embed()

        # restore sys.path
        sys.path = old_path[:]

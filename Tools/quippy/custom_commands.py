"""Overridden versions of distutils `clean` and `build_ext` commands,
and a new `test` command to find and run unittests."""

import sys, os, glob

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

        # just for the heck of it, try to remove the base build directory:
        # we might have emptied it right now, but if not we don't care
        if not self.dry_run:
            try:
                os.rmdir(self.build_base)
                log.info("removing '%s'", self.build_base)
            except OSError:
                pass

class build_ext(_build_ext):
    """Subclass of numpy.distutils.command.build_ext.build_ext which adds module_build_dir of
       included libraries to ext.module_dirs for all ext in self.extensions."""
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
                #ext.module_dirs.append(build_clib.build_clib)
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
        if self.inplace:
            cmd = self.reinitialize_command('scons')
            cmd.inplace = 1
        self.run_command('scons')


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
         "verbosity for unit tests")
        ]

    def initialize_options(self):
        self.build_base = 'build'
        # these are decided only after 'build_base' has its final value
        # (unless overridden by the user or client)
        self.test_dir = 'tests'
        self.test_prefix = 'test_'
        self.test_suffixes = None
        self.verbosity = None

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
        tests = TestLoader().loadTestsFromNames([self.test_prefix+case for case in self.test_suffixes])
        print 'Running tests with verbosity %d' % self.verbosity
        TextTestRunner(verbosity=self.verbosity).run(tests)

        # restore sys.path
        sys.path = old_path[:]



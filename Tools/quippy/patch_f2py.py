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

"""This module patches :mod:`numpy.f2py` at runtime, to customise the C code that
is generated. We make several changes to f2py:

  1. Allow the Fortran :cfunc:`present` function to work correctly with optional arguments.
     If an argument to an f2py wrapped function is optional and is not given, replace it
     with ``NULL``.
     
  2. Allow Fortran routines to raise a :exc:`RuntimeError` exception
     with a message by calling an external function
     :cfunc:`quippy_error_abort`. This is implemented using a
     :cfunc:`setjmp`/ :cfunc:`longjmp` trap.

  3. Allow Fortran routines to be interrupted with :kbd:`Ctrl+C` by installing a custom
     interrupt handler before the call into Fortran is made. After the Fortran routine
     returns, the previous interrupt handler is restored.

"""

import numpy
if not tuple([int(x) for x in numpy.__version__.split('.')[0:3]]) >= (1,2,1):
   raise ImportError('patch_f2py only tested with numpy version 1.2.1 or later, found version %s' % numpy.__version__)


import numpy.f2py.auxfuncs

import numpy.f2py.capi_maps
numpy.f2py.capi_maps.cformat_map['long_long'] = '%Ld'

import numpy.f2py.rules, numpy.f2py.cb_rules

numpy.f2py.rules.module_rules['modulebody'] = numpy.f2py.rules.module_rules['modulebody'].replace('#includes0#\n', """#includes0#
#include <libatoms.h>""")

numpy.f2py.rules.routine_rules['body'] = numpy.f2py.rules.routine_rules['body'].replace('\tvolatile int f2py_success = 1;\n', """\tvolatile int f2py_success = 1;
\tint setjmpvalue; /* James Kermode - for setjmp */
""")

numpy.f2py.rules.routine_rules['body'] = numpy.f2py.rules.routine_rules['body'].replace('#callfortranroutine#\n', """/* setjmp() exception handling added by James Kermode */
PyOS_sighandler_t _npy_sig_save;
_npy_sig_save = PyOS_setsig(SIGINT, quippy_error_abort_int_handler);
setjmpvalue = setjmp(environment_buffer);
if (setjmpvalue != 0) {
  PyOS_setsig(SIGINT, _npy_sig_save);
  PyErr_SetString(PyExc_RuntimeError, abort_message);
} else {
 #callfortranroutine#
 PyOS_setsig(SIGINT, _npy_sig_save);
}
/* End addition */
""")

from numpy.f2py.auxfuncs import *

numpy.f2py.auxfuncs.options['persistant_callbacks'] = True

# Disable callback argument cleanup so that callbacks can be called after function returns.
# This will lead to a small memory leak every time a function with callback arguments is called.
def persistant_callbacks(var):
   return True

numpy.f2py.cb_rules.cb_routine_rules['body'] = numpy.f2py.cb_rules.cb_routine_rules['body'].replace('capi_longjmp_ok = 1', 'capi_longjmp_ok = 0')

numpy.f2py.rules.arg_rules[7]['cleanupfrompyobj'] = {l_not(persistant_callbacks): numpy.f2py.rules.arg_rules[7]['cleanupfrompyobj'],
                                                     persistant_callbacks: '}'}

numpy.f2py.rules.arg_rules[8]['callfortran'] = {isintent_c:'#varname#,',
                                                l_and(isoptional,l_not(isintent_c)):'#varname#_capi == Py_None ? NULL : &#varname#,',
                                                l_and(l_not(isoptional),l_not(isintent_c)):'&#varname#,'}


numpy.f2py.rules.arg_rules[14]['callfortran'] = {isintent_c:'#varname#,',l_and(isoptional,l_not(isintent_c)):'#varname#_capi == Py_None ? NULL : &#varname#,',
                                                 l_and(l_not(isoptional),l_not(isintent_c)):'&#varname#,'}

numpy.f2py.rules.arg_rules[21]['callfortran'] = {isintent_out:'#varname#,', l_and(isoptional, l_not(isintent_out)):'#varname#_capi == Py_None ? NULL : #varname#,',
                                                 l_and(l_not(isoptional), l_not(isintent_out)): '#varname#,'}


numpy.f2py.rules.arg_rules[26]['callfortran'] = {isintent_out:'#varname#,', l_and(isoptional, l_not(isintent_out)):'#varname#_capi == Py_None ? NULL : #varname#,',
                                                 l_and(l_not(isoptional), l_not(isintent_out)): '#varname#,'}

numpy.f2py.rules.arg_rules[33]['frompyobj'].insert(2, {isoptional: 'if (#varname#_capi != Py_None) {'})
numpy.f2py.rules.arg_rules[33]['frompyobj'].insert(5, {isoptional: '}'})

del numpy.f2py.rules.arg_rules[33]['frompyobj'][6]
numpy.f2py.rules.arg_rules[33]['frompyobj'].insert(6, {l_not(isoptional): """\
\tif (capi_#varname#_tmp == NULL) {
\t\tif (!PyErr_Occurred())
\t\t\tPyErr_SetString(#modulename#_error,\"failed in converting #nth# `#varname#\' of #pyname# to C/Fortran array\" );
\t} else {
\t\t#varname# = (#ctype# *)(capi_#varname#_tmp->data);
"""})

numpy.f2py.rules.arg_rules[33]['frompyobj'].insert(7, {isoptional:"""\
\tif (#varname#_capi != Py_None && capi_#varname#_tmp == NULL) {
\t\tif (!PyErr_Occurred())
\t\t\tPyErr_SetString(#modulename#_error,\"failed in converting #nth# `#varname#\' of #pyname# to C/Fortran array\" );
\t} else {
\t\tif (#varname#_capi != Py_None) #varname# = (#ctype# *)(capi_#varname#_tmp->data);
"""})

def library_option(self, lib):
   if lib[0] == '-' or lib.endswith('.a'):
      return lib
   else:
      return '-l' + lib

import numpy.distutils.fcompiler
import numpy.distutils.ccompiler
numpy.distutils.fcompiler.FCompiler.library_option = library_option
numpy.distutils.unixccompiler.UnixCCompiler.library_option = library_option

# Replace distutils.ccompiler.CCompiler setup_compile and _prep_compile methods 
# with versions from Python 2.6.4 to avoid excessive recompilation of .o files
# (see, e.g. http://stackoverflow.com/questions/3145933/python-distutils-builds-extensions-differently-on-different-machines)

from types import *
from copy import copy
from distutils.errors import *
from distutils.spawn import spawn
from distutils.file_util import move_file
from distutils.dir_util import mkpath
from distutils.dep_util import newer_pairwise, newer_group
from distutils.util import split_quoted, execute
from distutils import log
from distutils.ccompiler import gen_preprocess_options
import os

def old_setup_compile(self, outdir, macros, incdirs, sources, depends,
                      extra):
     """Process arguments and decide which source files to compile.

     Merges _fix_compile_args() and _prep_compile().
     """
     if outdir is None:
         outdir = self.output_dir
     elif type(outdir) is not StringType:
         raise TypeError, "'output_dir' must be a string or None"

     if macros is None:
         macros = self.macros
     elif type(macros) is ListType:
         macros = macros + (self.macros or [])
     else:
         raise TypeError, "'macros' (if supplied) must be a list of tuples"

     if incdirs is None:
         incdirs = self.include_dirs
     elif type(incdirs) in (ListType, TupleType):
         incdirs = list(incdirs) + (self.include_dirs or [])
     else:
         raise TypeError, \
               "'include_dirs' (if supplied) must be a list of strings"

     if extra is None:
         extra = []

     # Get the list of expected output (object) files
     objects = self.object_filenames(sources,
                                     strip_dir=0,
                                     output_dir=outdir)
     assert len(objects) == len(sources)

     # XXX should redo this code to eliminate skip_source entirely.
     # XXX instead create build and issue skip messages inline

     if self.force:
         skip_source = {}            # rebuild everything
         for source in sources:
             skip_source[source] = 0
     elif depends is None:
         # If depends is None, figure out which source files we
         # have to recompile according to a simplistic check. We
         # just compare the source and object file, no deep
         # dependency checking involving header files.
         skip_source = {}            # rebuild everything
         for source in sources:      # no wait, rebuild nothing
             skip_source[source] = 1

         n_sources, n_objects = newer_pairwise(sources, objects)
         for source in n_sources:    # no really, only rebuild what's
             skip_source[source] = 0 # out-of-date
     else:
         # If depends is a list of files, then do a different
         # simplistic check.  Assume that each object depends on
         # its source and all files in the depends list.
         skip_source = {}
         # L contains all the depends plus a spot at the end for a
         # particular source file
         L = depends[:] + [None]
         for i in range(len(objects)):
             source = sources[i]
             L[-1] = source
             if newer_group(L, objects[i]):
                 skip_source[source] = 0
             else:
                 skip_source[source] = 1

     pp_opts = gen_preprocess_options(macros, incdirs)

     build = {}
     for i in range(len(sources)):
         src = sources[i]
         obj = objects[i]
         ext = os.path.splitext(src)[1]
         self.mkpath(os.path.dirname(obj))
         if skip_source[src]:
             log.debug("skipping %s (%s up-to-date)", src, obj)
         else:
             build[obj] = src, ext

     return macros, objects, extra, pp_opts, build

def old_prep_compile(self, sources, output_dir, depends=None):
  """Decide which souce files must be recompiled.

  Determine the list of object files corresponding to 'sources',
  and figure out which ones really need to be recompiled.
  Return a list of all object files and a dictionary telling
  which source files can be skipped.
  """
  # Get the list of expected output (object) files
  objects = self.object_filenames(sources, output_dir=output_dir)
  assert len(objects) == len(sources)

  if self.force:
      skip_source = {}            # rebuild everything
      for source in sources:
          skip_source[source] = 0
  elif depends is None:
      # If depends is None, figure out which source files we
      # have to recompile according to a simplistic check. We
      # just compare the source and object file, no deep
      # dependency checking involving header files.
      skip_source = {}            # rebuild everything
      for source in sources:      # no wait, rebuild nothing
          skip_source[source] = 1

      n_sources, n_objects = newer_pairwise(sources, objects)
      for source in n_sources:    # no really, only rebuild what's
          skip_source[source] = 0 # out-of-date
  else:
      # If depends is a list of files, then do a different
      # simplistic check.  Assume that each object depends on
      # its source and all files in the depends list.
      skip_source = {}
      # L contains all the depends plus a spot at the end for a
      # particular source file
      L = depends[:] + [None]
      for i in range(len(objects)):
          source = sources[i]
          L[-1] = source
          if newer_group(L, objects[i]):
              skip_source[source] = 0
          else:
              skip_source[source] = 1

  return objects, skip_source

def find_library_file(self, dirs, lib, debug=0):
   return os.path.join(dirs[0], lib)

from distutils.ccompiler import CCompiler
CCompiler._setup_compile = old_setup_compile
CCompiler._prep_comile = old_prep_compile
CCompiler.find_library_file = find_library_file   

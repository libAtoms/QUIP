"""This module patches :mod:`numpy.f2py` at runtime, to customise the C code that
is generated. We make several changes to f2py:

  1. Allow the Fortran :cfunc:`present` function to work correctly with optional arguments.
     If an argument to an f2py wrapped function is optional and is not given, replace it
     with a ``NULL`` for scalars and arrays, or with a pointer to ``NULL`` for derived-types.

  2. Allow Fortran routines to raise a :exc:`RuntimeError` exception with a message
     by calling the :cfunc:`system_abort` function. This is implemented using a :cfunc:`setjmp`/
     :cfunc:`longjmp` trap.

  3. Allow Fortran routines to be interrupted with :kbd:`Ctrl+C` by installing a custom
     interrupt handler before the call into Fortran is made. After the Fortran routine
     returns, the previous interrupt handler is restored.

"""

import numpy
if not tuple([int(x) for x in numpy.__version__.split('.')]) >= (1,2,1):
   raise ImportError('patch_f2py only tested with numpy version 1.2.1 or later')

from numpy.f2py.rules import f2py_version
already_patched = 'patched_JRK' in f2py_version

if not already_patched:
   import numpy.f2py.auxfuncs

   numpy.f2py.auxfuncs.isderivedtypepointer = lambda var: 'typename' in var and 'attrspec' in var and 'pointer' in var['attrspec']

   import numpy.f2py.capi_maps
   numpy.f2py.capi_maps.cformat_map['long_long'] = '%Ld'

   import numpy.f2py.rules

   numpy.f2py.rules.module_rules['modulebody'] = numpy.f2py.rules.module_rules['modulebody'].replace('#includes0#\n', """#includes0#

   /* Added by James Kermode */
   #include <setjmp.h>
   jmp_buf environment_buffer;
   char abort_message[1024];

   void system_abort_(char *message, int len)
   {
     strncpy(abort_message, message, len);
     abort_message[len] = '\\0';
     longjmp(environment_buffer,0);
   }

   void system_abort_int_handler(int signum)
   {
     char message[] = "Interrupt occured";
     system_abort_(message, strlen(message));
   }
   /* end addition */
   """)

   numpy.f2py.rules.routine_rules['body'] = numpy.f2py.rules.routine_rules['body'].replace('\tvolatile int f2py_success = 1;\n', """\tvolatile int f2py_success = 1;
   \tint setjmpvalue; /* James Kermode - for setjmp */
   """)

   numpy.f2py.rules.routine_rules['body'] = numpy.f2py.rules.routine_rules['body'].replace('#callfortranroutine#\n', """/* setjmp() exception handling added by James Kermode */
   PyOS_sighandler_t _npy_sig_save;
   _npy_sig_save = PyOS_setsig(SIGINT, system_abort_int_handler);
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

   numpy.f2py.rules.arg_rules[8]['decl'] = {l_or(l_not(isoptional),
                                                 l_not(isderivedtypepointer)): '\t#ctype# #varname# = 0;',
                                            l_and(isoptional, isderivedtypepointer): '\t#ctype# #varname# = 0;\n\tvoid *#varname#_nullptr = NULL;'}


   numpy.f2py.rules.arg_rules[8]['callfortran'] = {isintent_c:'#varname#,',
                                                   l_and(isoptional,l_not(isintent_c),l_not(isderivedtypepointer)):'#varname#_capi == Py_None ? NULL : &#varname#,',
                                                   l_and(isoptional,l_not(isintent_c),isderivedtypepointer):'#varname#_capi == Py_None ? (#ctype#*)&#varname#_nullptr : &#varname#,',
                                                   l_and(l_not(isoptional),l_not(isintent_c)):'&#varname#,'}


   numpy.f2py.rules.arg_rules[14]['callfortran'] = {isintent_c:'#varname#,',l_and(isoptional,l_not(isintent_c)):'#varname#_capi == Py_None ? NULL : &#varname#,',
                                                    l_and(l_not(isoptional),l_not(isintent_c)):'&#varname#,'}

   numpy.f2py.rules.arg_rules[21]['callfortran'] = {isintent_out:'#varname#,', l_and(isoptional, l_not(isintent_out)):'#varname#_capi == Py_None ? NULL : #varname#,',
                                                    l_and(l_not(isoptional), l_not(isintent_out)): '#varname#,'}


   numpy.f2py.rules.arg_rules[26]['callfortran'] = {isintent_out:'#varname#,', l_and(isoptional, l_not(isintent_out)):'#varname#_capi == Py_None ? NULL : #varname#,',
                                                    l_and(l_not(isoptional), l_not(isintent_out)): '#varname#,'}

   del numpy.f2py.rules.arg_rules[33]['frompyobj'][4]
   numpy.f2py.rules.arg_rules[33]['frompyobj'].insert(4, {l_not(isoptional): """\
   \tif (capi_#varname#_tmp == NULL) {
   \t\tif (!PyErr_Occurred())
   \t\t\tPyErr_SetString(#modulename#_error,\"failed in converting #nth# `#varname#\' of #pyname# to C/Fortran array\" );
   \t} else {
   \t\t#varname# = (#ctype# *)(capi_#varname#_tmp->data);
   """})

   numpy.f2py.rules.arg_rules[33]['frompyobj'].insert(5, {isoptional:"""\
   {
   \t\tPyErr_Clear();
   if (capi_#varname#_tmp != NULL)
   \t\t#varname# = (#ctype# *)(capi_#varname#_tmp->data);
   """})

   def library_option(self, lib):
      if lib[0] == '-':
         return lib
      else:
         return '-l' + lib

   import numpy.distutils.fcompiler
   numpy.distutils.fcompiler.FCompiler.library_option = library_option

   

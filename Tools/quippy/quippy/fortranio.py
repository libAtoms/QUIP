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

import _quippy
from quippy import InOutput, print_

# Create InOutput objects associated with Fortran stdout and stderr                
mainlog_ptr, errorlog_ptr = _quippy.qp_get_mainlog_errorlog_ptr()
mainlog  = InOutput(fpointer=mainlog_ptr, finalise=False)
errorlog = InOutput(fpointer=errorlog_ptr, finalise=False)
del mainlog_ptr, errorlog_ptr

class FortranWriter:
   def __init__(self, fortran_file):
      self.fortran_file = fortran_file
      self.saved_prefix = None
      self.leading_space = ''
      
   def write(self, text):
      if text.startswith(' '):
         self.leading_space += ' '*(len(text) - len(text.lstrip()))
         text = text.lstrip()
      for line in text.splitlines(True):
         print_(self.leading_space+line.rstrip(), file=self.fortran_file, nocr=not line.endswith('\n'))
         self.leading_space = ''
         if self.saved_prefix is not None and line.endswith('\n'):
            self.fortran_file.prefix = self.saved_prefix
            self.saved_prefix = None
         if self.saved_prefix is None and not line.endswith('\n'):
            self.saved_prefix = self.fortran_file.prefix
            self.fortran_file.prefix = ''

# create file-like objects associated with Fortran mainlog and errorlog
fortran_stdout = FortranWriter(mainlog)
fortran_stderr = FortranWriter(errorlog)

# convenience wrapper functions which behave like python "print" statement
def fortran_print(*args):
   """Print `args` to mainlog using Fortran I/O routines. Respects prefix and verbosity settings"""
   print >>fortran_stdout, " ".join([str(a) for a in args])

def fortran_print_error(*args):
   """Print `args` to errorlog using Fortran I/O routines. Respects prefix and verbosity settings"""
   print >>fortran_stderr, " ".join([str(a) for a in args])

del _quippy

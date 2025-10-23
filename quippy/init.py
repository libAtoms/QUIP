# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Portions of this code were written by
# HQ X     Tamas K. Stenczel, James Kermode
# HQ X
# HQ X   Copyright 2019
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   https://warwick.ac.uk/fac/sci/eng/staff/jrk
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
"""quippy package

Maintained by James Kermode <j.r.kermode@warwick.ac.uk>

Contains python bindings to the libAtoms/QUIP Fortran 95 codes
<http://libatoms.github.org/QUIP>. """

# IMPORTANT: Set RTLD_GLOBAL before any imports to ensure symbol visibility
# This is critical for runtime symbol resolution: libAtoms.so has undefined
# references to f90wrap_abort_ which must be resolved from _quippy.so
# By default Python loads extensions with RTLD_LOCAL, hiding symbols
import sys
import os
_quippy_dlopen_flags = sys.getdlopenflags()
sys.setdlopenflags(os.RTLD_NOW | os.RTLD_GLOBAL)

import quippy.convert
import quippy.potential
import quippy.descriptors
import quippy.nye_tensor
import quippy.gap_tools
import atexit

# Reference values of .true. and .false. from Fortran
QUIPPY_TRUE = quippy.system_module.reference_true()
QUIPPY_FALSE = quippy.system_module.reference_false()


def quippy_cleanup():
    try:
        quippy.system_module.verbosity_pop()
        quippy.system_module.system_finalise()
    except AttributeError:
        pass


quippy.system_module.system_initialise(-1, quippy_running=QUIPPY_TRUE)
quippy.system_module.verbosity_push(0)
atexit.register(quippy_cleanup)

# Restore original dlopen flags after all quippy modules are loaded
sys.setdlopenflags(_quippy_dlopen_flags)

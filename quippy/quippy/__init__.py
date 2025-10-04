from __future__ import print_function, absolute_import, division
import quippy._quippy
import f90wrap.runtime
import logging
import numpy
import warnings
import quippy.tb_module
import quippy.dictionary_module
import quippy.potential_module
import quippy.descriptors_module
import quippy.potential_simple_module
import quippy.system_module
import quippy.atoms_types_module
import quippy.nye_tensor_module
import quippy.dynamicalsystem_module
import quippy.atoms_module

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

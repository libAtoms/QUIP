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

from quippy.atoms import *
from quippy.structures import *
from quippy.farray import fidentity, fzeros, frange, farray

import numpy as np
try:
    import atomeye
except ImportError:
    pass

__all__ = ['J_PER_M2', 'surface_energy']

from quippy.units import ELEM_CHARGE

J_PER_M2 = ELEM_CHARGE*1.0e20

def surface_energy(pot, bulk, surface, dir=2):

    if not hasattr(bulk, 'energy'):
        bulk.calc_connect()
        pot.calc(bulk, energy=True, force=True)

    bulk_energy_per_sio2 = bulk.energy/(bulk.z == 14).count()

    surface.calc_connect()
    pot.calc(surface, energy=True, force=True)

    area = surface.cell_volume()/surface.lattice[dir,dir]

    return (surface.energy - (surface.z == 14).count()*bulk_energy_per_sio2)/(2.0*area)


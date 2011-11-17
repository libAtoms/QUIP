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

from quippy.atoms import Atoms, atoms_reader
from quippy.farray import fzeros, s2a, farray
from quippy.table import TABLE_STRING_LENGTH
from quippy.units import MASSCONVERT

__all__ = ['CP2KOutputReader']

@atoms_reader('cp2k_output')
def CP2KOutputReader(fh, module='QUICKSTEP', type_map=None, kind_map=None):
    if type_map is None:
        type_map = {}
    if kind_map is None:
        kind_map = {}

    if module == 'QS': module = 'QUICKSTEP'
    if module not in ('FIST', 'QUICKSTEP', 'QM/MM'):
        raise ValueError('Unknown CP2K module %s' % module)
    
    opened = False
    filename = fh
    if type(fh) == type(''):
        fh = open(fh, 'r')
        opened = True
    lines = fh.readlines()
    if opened:
        fh.close()

    cell_lines = [i for i,line in enumerate(lines) if line.startswith(" CELL| Vector a")]
    if cell_lines == []:
        raise ValueError("Cannot find cell in file %s" % filename)

    if module == 'FIST':
        cell_line = cell_lines[0]
    elif module == 'QM/MM':
        try:
            cell_line = cell_lines[1]
        except IndexError:
            raise ValueError("Cannot file QMMM cell in file %s" % filename)
    elif module == 'QUICKSTEP':
        try:
            cell_line = cell_lines[2]
        except IndexError:
            raise ValueError("Cannot file QUICKSTEP cell in file %s" % filename)
    else:
        raise ValueError("Don't know how to find cell for module %s" % module)

    lattice = fzeros((3,3))
    for i in [0,1,2]:
        lattice[:,i+1] = [float(c) for c in lines[cell_line+i].split()[4:7]]
        
    try:
        start_line = lines.index(" MODULE %s:  ATOMIC COORDINATES IN angstrom\n" % module)
    except ValueError:
        raise ValueError("Cannot find atomic positions for module %s in file %s" % (module, filename))

    kinds = []
    species = []
    Zs = []
    pos = []
    masses = []
    Zeffs = []
    types = []
    qeffs = []
    for line in lines[start_line+4:]:
        if line.strip() == '':
            break
        if module == 'FIST':
            atom, kind, typ, x, y, z, qeff, mass = line.split()
            types.append(typ)
            Z = type_map.get(typ,0)
            kind = int(kind)
            if Z == 0:
                Z = kind_map.get(kind,0)
            Zs.append(Z)
            qeffs.append(float(qeff))
        else:
            atom, kind, sp, Z, x, y, z, Zeff, mass = line.split()
            species.append(sp)
            Zs.append(int(Z))
            Zeffs.append(float(Zeff))
        kinds.append(int(kind))
        pos.append([float(x),float(y),float(z)])
        masses.append(float(mass))

    at = Atoms(n=len(kinds), lattice=lattice)
    at.pos[...] = farray(pos).T
    at.set_atoms(Zs)
    at.add_property('mass', farray(masses)*MASSCONVERT)
    at.add_property('kind', kinds)
    if module == 'FIST':
        at.add_property('type', ' '*TABLE_STRING_LENGTH)
        at.add_property('qm', False)
        at.qm[:] = (at.type.stripstrings() == '_QM_') | (at.type.stripstrings() == '_LNK')
        at.type[...] = s2a(types, TABLE_STRING_LENGTH)
        at.add_property('qeff', qeffs)
    else:
        at.species[...] = s2a(species, TABLE_STRING_LENGTH)
        at.add_property('zeff', Zeffs)

    yield at

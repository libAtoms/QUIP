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

import logging, StringIO, sys
from math import pi
import numpy as np

from quippy.atoms import Atoms
from quippy.io import atoms_reader, AtomsReaders, AtomsWriters
from quippy.periodictable import ElementName
from quippy.dictionary import (Dictionary, PROPERTY_INT, PROPERTY_REAL,
                               PROPERTY_STR, PROPERTY_LOGICAL,
                               T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                               T_CHAR, T_LOGICAL, T_INTEGER_A,
                               T_REAL_A, T_COMPLEX_A, T_CHAR_A,
                               T_LOGICAL_A, T_INTEGER_A2, T_REAL_A2)
from quippy.table import TABLE_STRING_LENGTH
from quippy.io import *
from quippy.ordereddict import *
from quippy.farray import *


__all__ = ['PuPyXYZReader', 'PuPyXYZWriter']

def _props_dtype(props):
    "Return a record array dtype for the specified properties (default all)"

    result = []
    fmt_map = {'R':'d',
               'I': 'i',
               'S': 'S' + str(TABLE_STRING_LENGTH),
               'L': 'bool'}

    for prop in props:
        ptype, cols = props[prop]
        if cols == 1:
            result.append((prop,fmt_map[ptype]))
        else:
            for c in range(cols):
                result.append((prop+str(c),fmt_map[ptype]))

    return np.dtype(result)

def parse_properties(prop_str):
    properties = OrderedDict()
    if prop_str.startswith('pos:R:3'):
        prop_str = 'species:S:1:'+prop_str

    fields = prop_str.split(':')

    for name, ptype, cols in zip(fields[::3], fields[1::3], map(int,fields[2::3])):
        if ptype not in ('R', 'I', 'S', 'L'):
            raise ValueError('Unknown property type: '+ptype)
        properties[name] = (ptype, cols)

    return properties


@atoms_reader('pupyxyz')
def PuPyXYZReader(xyz, format=None):
    "Read from extended XYZ filename, string or open file."
    from quippy import Table

    def _getconv(dtype):
        typ = dtype.type
        if issubclass(typ, np.bool_):
            return lambda x: {'T' : True, 'F' : False, 'True' : True, 'False' : False}.get(x)
        if issubclass(typ, np.integer):
            return int
        elif issubclass(typ, np.floating):
            return float
        elif issubclass(typ, np.complex):
            return complex
        else:
            return str

    opened = False
    if type(xyz) == type(''):
        if '\n' in xyz and xyz[:xyz.index('\n')].isdigit():
            # string looks like an embedded XYZ file
            xyz = iter(xyz.split('\n'))
        else:
            xyz = open(xyz,'r')
            opened = True
    else:
        xyz = iter(xyz)

    while True:
        line = xyz.next()
        if not line: raise StopIteration

        n = int(line.strip())
        comment = (xyz.next()).strip()

        # Parse comment line
        params = Dictionary(comment)

        if not 'Properties' in params:
            raise ValueError('Properties missing from comment line')

        properties = parse_properties(params['Properties'])
        del params['Properties']

        # Get lattice
        if not 'Lattice' in params:
            raise ValueError('No lattice found in xyz file')

        lattice = params.get_value('Lattice') # make a copy
        del params['Lattice']

        # Get speccial params entries:
        #   nneightol, cutoff, cutoff_break, cutoff_factor, cutoff_break_factor, pbc
        special_params = {}
        for key in ['nneightol', 'cutoff', 'cutoff_break',
                    'cutoff_factor', 'cutoff_break_factor',
                    'pbc']:
            special_params[key] = params.get(key)
            if key in params:
                del params[key]

        props_dtype = _props_dtype(properties)

        converters = [_getconv(props_dtype.fields[name][0]) \
                      for name in props_dtype.names]

        X = []
        for i,line in enumerate(xyz):
            vals = line.split()
            row = tuple([converters[j](val) for j,val in enumerate(vals)])
            X.append(row)
            if i == n-1: break # Only read self.n lines

        try:
            data = np.array(X,props_dtype)
        except TypeError:
            raise IOError('Badly formatted data, or end of file reached before end of frame')

        # Empty dictionary passed for properties to avoid creating pos, species, z.
        at = Atoms(n=n,lattice=lattice,params=params,properties={})

        # Set the special attributes nneightol, cutoff, cutoff_break, use_uniform_cutoff
        if special_params.get('nneightol'):
            at.nneightol = special_params['nneightol']
        if special_params.get('cutoff'):
            at.set_cutoff(special_params['cutoff'], special_params['cutoff_break'])
        elif special_params.get('cutoff_factor'):
            at.set_cutoff_factor(special_params['cutoff_factor'], special_params['cutoff_break_factor'])
        if special_params.get('pbc'):
            at.pbc = special_params['pbc']

        for p in properties:
            ptype, cols = properties[p]
            if cols == 1:
                value = data[p]
                if value.dtype.kind == 'S':
                    value = s2a(value, TABLE_STRING_LENGTH).T
            else:
                value = np.vstack([data[p+str(c)] for c in range(cols) ])
            at.add_property(p, value, overwrite=True)

        if not at.has_property('Z')  and not at.has_property('species'):
            raise ValueError('Atoms read from XYZ has neither Z nor species')
        elif at.has_property('Z') and not at.has_property('species'):
            at.add_property('species', ' '*TABLE_STRING_LENGTH)
            at.set_atoms(at.z)
        elif at.has_property('species') and not at.has_property('z'):
            at.add_property('Z', 0)
            at.z[:] = [ElementName.index(sp) for sp in at.species.stripstrings()]

        yield at

    if opened: xyz.close()


class PuPyXYZWriter(object):
    "Write atoms in extended XYZ format. xyz can be a filename or open file"

    def __init__(self, xyz):
        self.opened = False
        self.string = False
        if type(xyz) == type(''):
            if xyz == 'stdout':
                self.xyz = sys.stdout
            elif xyz == 'string':
                self.xyz = StringIO.StringIO()
                self.string = True
            else:
                self.xyz = open(xyz, 'w')
                self.opened = True
        else:
            self.xyz = xyz

    def write(self, at, properties=None):

        def _getfmt(dtype):
            typ = dtype.type
            if issubclass(typ, np.bool_):
                return ' %.1s '
            if issubclass(typ, np.integer):
                return '%8d '
            elif issubclass(typ, np.floating):
                return '%16.8f '
            elif issubclass(typ, np.complex):
                return '(%f,%f) '
            else:
                return '%s '

        def ncols(props, name):
            t, l, (l1,l2) = props.get_type_and_size(name)
            if t in (T_INTEGER_A, T_REAL_A, T_LOGICAL_A, T_CHAR_A):
                return 1
            elif t in (T_INTEGER_A2, T_REAL_A2):
                return l1
            else:
                raise TypeError('bad property type %d' % t)

        def properties_comment(self, properties=None, params=None):
            if properties is None:
                props = self.properties.keys()
            else:
                props = properties

            if params is None:
                params = self.params

            pmap = { T_INTEGER_A: 'I',
                     T_REAL_A: 'R',
                     T_CHAR_A: 'S',
                     T_LOGICAL_A: 'L',
                     T_INTEGER_A2: 'I',
                     T_REAL_A2: 'R'}

            lattice_str = 'Lattice="' + ' '.join(map(str, np.reshape(self.lattice,9,order='F'))) + '"'

            props_str =  ':'.join(map(':'.join,
                                      zip(props,
                                          [pmap[self.properties.get_type_and_size(k)[0]] for k in props],
                                          [str(ncols(self.properties, k)) for k in props])))

            return props_str, lattice_str+' Properties='+props_str+' '+str(params)


        if properties is None:
            props = at.properties.keys()
        else:
            props = properties

        if 'species' in props:
            i = props.index('species')
            props[0], props[i] = props[i], props[0]

        if 'pos' in props:
            i = props.index('pos')
            props[1], props[i] = props[i], props[1]

        species = getattr(at, props[0].lower())
        if species.shape != (TABLE_STRING_LENGTH, at.n) or species.dtype.kind != 'S':
            raise ValueError('First property must be species like')

        pos = getattr(at, props[1].lower())
        if pos.shape != (3, at.n) or pos.dtype.kind != 'f':
            raise ValueError('Second property must be position like')

        # Set extra params entries for cutoff etc.
        params = at.params.copy()
        params['nneightol'] = at.nneightol
        params['cutoff'] = at.cutoff
        params['pbc'] = at.pbc

        props_str, comment = properties_comment(at, props, params)

        subprops = OrderedDict.frompairs([(k, at.properties[k]) for k in props])
        data = np.zeros(at.n, _props_dtype(parse_properties(props_str)))

        for prop in props:
            value = getattr(at, prop)
            if at.properties.get_type_and_size(prop)[0] in (T_INTEGER_A2, T_REAL_A2):
                for c in range(value.shape[0]):
                    data[prop+str(c)] = value[c+1,:]
            else:
                if value.dtype.kind == 'S':
                    value = a2s(value)
                data[prop] = value

        format = ''.join([_getfmt(data.dtype.fields[name][0]) for name in data.dtype.names])+'\n'

        self.xyz.write('%d\n' % at.n)
        self.xyz.write('%s\n' % comment)
        for i in range(at.n):
            self.xyz.write(format % tuple(data[i]))

        if self.string: return self.xyz.getvalue()

    def close(self):
        if self.opened: self.xyz.close()

AtomsReaders['pupyxyz'] = PuPyXYZReader
AtomsWriters['pupyxyz'] = AtomsWriters['stdout'] = PuPyXYZWriter

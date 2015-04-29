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

import logging, StringIO
from math import pi

from quippy.atoms import make_lattice, get_lattice_params
from quippy.io import atoms_reader, AtomsReaders, AtomsWriters
from quippy.periodictable import ElementName
from quippy.dictionary import (Dictionary, PROPERTY_INT, PROPERTY_REAL,
                               PROPERTY_STR, PROPERTY_LOGICAL,
                               T_NONE, T_INTEGER, T_REAL, T_COMPLEX,
                               T_CHAR, T_LOGICAL, T_INTEGER_A,
                               T_REAL_A, T_COMPLEX_A, T_CHAR_A,
                               T_LOGICAL_A, T_INTEGER_A2, T_REAL_A2)
from quippy.table import TABLE_STRING_LENGTH
from quippy.ordereddict import *
from quippy.farray import *
from quippy import netcdf_file

__all__ = ['NetCDFReader']

def netcdf_dimlen(obj, name):
    """Return length of dimension 'name'. Works for both netCDF4 and pupynere."""
    n = obj.dimensions[name]
    try:
        return len(n)
    except TypeError:
        return n

@atoms_reader(netcdf_file)
@atoms_reader('nc')
def NetCDFReader(source, frame=None, start=0, stop=None, step=1, format=None):

    opened = False
    if isinstance(source, str):
        opened = True
        source = netcdf_file(source)

    from quippy import Atoms

    DEG_TO_RAD = pi/180.0

    remap_names = {'coordinates': 'pos',
                   'velocities': 'velo',
                   'cell_lengths': None,
                   'cell_angles': None,
                   'cell_lattice': None,
                   'cell_rotated': None}

    prop_type_to_value = {PROPERTY_INT: 0,
                          PROPERTY_REAL: 0.0,
                          PROPERTY_STR: "",
                          PROPERTY_LOGICAL: False}

    prop_dim_to_ncols = {('frame','atom','spatial'): 3,
                         ('frame','atom','label'): 1,
                         ('frame', 'atom'): 1}

    if frame is not None:
        start = frame
        stop = frame+1
        step = 1
    else:
        if stop is None:
            stop = source.variables['cell_lengths'].shape[0]

    for frame in range(start, stop, step):
        cl = source.variables['cell_lengths'][frame]
        ca = source.variables['cell_angles'][frame]
        lattice = make_lattice(cl[0],cl[1],cl[2],ca[0]*DEG_TO_RAD,ca[1]*DEG_TO_RAD,ca[2]*DEG_TO_RAD)

        at = Atoms(n=netcdf_dimlen(source, 'atom'), lattice=lattice, properties={})

        for name, var in source.variables.iteritems():
            name = remap_names.get(name, name)

            if name is None:
                continue

            name = str(name) # in case it's a unicode string

            if 'frame' in var.dimensions:
                if 'atom' in var.dimensions:
                    # It's a property
                    value = var[frame]
                    if value.dtype.kind != 'S': value = value.T
                    at.add_property(name, value)
                else:
                    # It's a param
                    if var.dimensions == ('frame','string'):
                        # if it's a single string, join it and strip it
                        at.params[name] = ''.join(var[frame]).strip()
                    else:
                        if name == 'cutoff':
                            at.cutoff = var[frame]
                        elif name == 'cutoff_skin':
                            at.cutoff_skin = var[frame]
                        elif name == 'nneightol':
                            at.nneightol = var[frame]
                        else:
                            at.params[name] = var[frame].T

        if 'cell_rotated' in source.variables:
            cell_rotated = source.variables['cell_rotated'][frame]
            orig_lattice = source.variables['cell_lattice'][frame]
            #if cell_rotated == 1:
            at.set_lattice(orig_lattice, True)

        yield at

    def close(self):
        if self.opened: source.close()

class NetCDFWriter(object):

    def __init__(self, dest, frame=None, append=False, netcdf4=True):
        self.opened = False
        self.frame = frame
        if isinstance(dest, str):
            self.opened = True
            try:
                FORMAT = {True:  'NETCDF4',
                          False: 'NETCDF3_CLASSIC'}
                MODE = {True:  'a',
                        False: 'w'}

                self.dest = netcdf_file(dest, MODE[append], format=FORMAT[netcdf4])
            except ValueError:
                self.dest = netcdf_file(dest, 'w')
        else:
            self.dest = dest

    def write(self, at, **kwargs):
        remap_names_rev = {'pos': 'coordinates',
                           'velocities': 'velo'}


        prop_type_ncols_to_dtype_dim = {(T_INTEGER_A2,3):       ('i', ('frame','atom','spatial')),
                                        (T_REAL_A2,3):     ('d', ('frame','atom','spatial')),
                                        (T_LOGICAL_A,3):  ('d', ('frame','atom','spatial')),
                                        (T_INTEGER_A,1):      ('i', ('frame', 'atom')),
                                        (T_REAL_A, 1):    ('d', ('frame', 'atom')),
                                        (T_LOGICAL_A, 1): ('i', ('frame', 'atom')),
                                        (T_CHAR_A,1):      ('S1', ('frame','atom','label'))}


        param_type_to_dtype_dim = {T_INTEGER:   ('i', ('frame',)),
                                   T_REAL:      ('d', ('frame',)),
                                   T_CHAR:      ('S1', ('frame', 'string')),
                                   T_LOGICAL:   ('i', ('frame',)),
                                   T_INTEGER_A: ('i', ('frame', 'spatial')),
                                   T_REAL_A:    ('d', ('frame', 'spatial')),
                                   T_LOGICAL_A: ('i', ('frame', 'spatial'))}

        # Temporary hack to retain compatibility with old NetCDF property codes
        convert_ptype = {
           T_INTEGER_A2: PROPERTY_INT,
           T_REAL_A2: PROPERTY_REAL,
           T_LOGICAL_A: PROPERTY_LOGICAL,
           T_INTEGER_A: PROPERTY_INT,
           T_REAL_A:   PROPERTY_REAL,
           T_LOGICAL_A: PROPERTY_LOGICAL,
           T_CHAR_A : PROPERTY_STR
           }

        if self.dest.dimensions == {}:
            self.dest.createDimension('frame', None)
            self.dest.createDimension('spatial', 3)
            self.dest.createDimension('atom', at.n)
            self.dest.createDimension('cell_spatial', 3)
            self.dest.createDimension('cell_angular', 3)
            self.dest.createDimension('label', 10)
            self.dest.createDimension('string', 1024)

            self.dest.createVariable('spatial', 'S1', ('spatial',))
            self.dest.variables['spatial'][:] = list('xyz')

            self.dest.createVariable('cell_spatial', 'S1', ('cell_spatial',))
            self.dest.variables['cell_spatial'][:] = list('abc')

            self.dest.createVariable('cell_angular', 'S1', ('cell_angular', 'label'))
            self.dest.variables['cell_angular'][0,:] = list('alpha     ')
            self.dest.variables['cell_angular'][1,:] = list('beta      ')
            self.dest.variables['cell_angular'][2,:] = list('gamma     ')

            self.dest.createVariable('cell_lengths', 'd', ('frame', 'spatial'))
            self.dest.createVariable('cell_angles', 'd', ('frame', 'spatial'))

        if self.frame is None:
            self.frame = netcdf_dimlen(self.dest, 'frame')
            if self.frame is None: self.frame = 0


        assert at.n == netcdf_dimlen(self.dest, 'atom')

        a, b, c, alpha, beta, gamma = get_lattice_params(at.lattice)

        RAD_TO_DEG = 180.0/pi

        self.dest.variables['cell_lengths'][self.frame] = (a, b, c)
        self.dest.variables['cell_angles'][self.frame] = (alpha*RAD_TO_DEG, beta*RAD_TO_DEG, gamma*RAD_TO_DEG)

        for origname in at.properties:
            name = remap_names_rev.get(origname, origname)
            ptype, s1, (s2, s3) = at.properties.get_type_and_size(origname)
            if ptype in (T_INTEGER_A, T_REAL_A, T_LOGICAL_A, T_CHAR_A):
                ncols = 1
            elif ptype in (T_INTEGER_A2, T_REAL_A2):
                ncols = s2
            else:
                raise TypeError('bad property type %d' % ptype)
            dtype, dims = prop_type_ncols_to_dtype_dim[(ptype, ncols)]

            if not name in self.dest.variables:
                self.dest.createVariable(name, dtype, dims)
                self.dest.variables[name].type = convert_ptype[ptype]

            assert self.dest.variables[name].dimensions == dims
            self.dest.variables[name][self.frame] = getattr(at, origname.lower()).T

        params = at.params.copy()
        params['nneightol'] = at.nneightol
        params['cutoff'] = at.cutoff

        for name in params.keys():
            t, s, shape, = params.get_type_and_size(name)
            dtype, dims = param_type_to_dtype_dim[t]

            if not name in self.dest.variables:
                self.dest.createVariable(name, dtype, dims)
                self.dest.variables[name].type = t

            assert self.dest.variables[name].type == t
            assert self.dest.variables[name].dimensions == dims

            self.dest.variables[name][self.frame] = params[name]

        self.frame += 1


    def close(self):
        if self.opened: self.dest.close()


AtomsWriters[netcdf_file] = NetCDFWriter
if not 'nc' in AtomsWriters: AtomsWriters['nc'] = NetCDFWriter

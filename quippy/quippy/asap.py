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


"""
This module contains utility functions for use with the ASAP code, which is
developed by Paul Tangney and Sandro Scandolo.
"""

from __future__ import print_function, absolute_import, division

import time, os, itertools, sys
from .farray import farray, fzeros, frange

from quippy.atoms import Atoms
from quippy.io import AtomsReaders, AtomsWriters, atoms_reader
from quippy.periodictable import ElementMass, atomic_number
from quippy.units import MASSCONVERT, BOHR, HARTREE, RYDBERG, GPA

__all__ = ['PosCelWriter', 'PosCelReader']


class PosCelWriter(object):

    def __init__(self, basename=None, pos='pos.in', cel='cel.in',
                 force='force.in', energy='energy.in', stress='stress.in',
                 step_name='', species_map={'O':1, 'Si':2}, cel_angstrom=False,
                 pos_angstrom=False, rydberg=True):
        if basename == 'stdout':
            pos = cel = energy = stress = force = sys.stdout
        elif basename is not None:
            basename = os.path.splitext(basename)[0]
            pos = '%s.pos' % basename
            cel = '%s.cel' % basename
            energy = '%s.ene' % basename
            stress = '%s.str' % basename
            force = '%s.for' % basename
        self.pos = pos
        self.cel = cel
        self.force = force
        self.energy = energy
        self.stress = stress
        self.step_name = step_name
        self.species_map = species_map
        self.cel_angstrom = cel_angstrom
        self.pos_angstrom = pos_angstrom
        self.rydberg = rydberg
        self.it = 0

    def write(self, at):

        self.doenergy = hasattr(at, 'energy')
        self.doforce  = hasattr(at, 'force')
        self.dostress = hasattr(at, 'virial')

        if isinstance(self.pos, str): self.pos = open(self.pos, 'w')
        if isinstance(self.cel, str): self.cel = open(self.cel, 'w')
        if self.doenergy and isinstance(self.energy, str): self.energy = open(self.energy, 'w')
        if self.doforce  and isinstance(self.force,  str): self.force  = open(self.force,  'w')
        if self.dostress and isinstance(self.stress, str): self.stress = open(self.stress, 'w')

        objind = 1

        self.pos.write('\n')# % (self.step_name, self.it))
        for i in frange(at.n):
            p = at.pos[:,i].copy()
            if not self.pos_angstrom: p /= BOHR
            self.pos.write('%20.10e%20.10e%20.10e%4d%4d X\n' % (p[1], p[2], p[3], 
                           self.species_map[at.species[i].stripstrings()], objind))

        self.cel.write(' %s %d\n' % (self.step_name, self.it))
        for i in (1,2,3):
            L = at.lattice[:,i].copy()
            if not self.cel_angstrom: L /= BOHR
            self.cel.write('%20.10e%20.10e%20.10e\n' % (L[1], L[2], L[3]))

        if self.doenergy:
            e = at.energy
            if self.rydberg:
                e /= RYDBERG
            self.energy.write('%20.10f Ry\n' % e)

        if self.doforce:
            self.force.write('\n')
            for i in frange(at.n):
                f = at.force[i].copy()
                if self.rydberg:
                    f /= (RYDBERG/BOHR)
                self.force.write('%20.10e%20.10e%20.10e%4d%4d X\n' % (f[1], f[2], f[3], self.species_map[at.species[i].stripstrings()], objind))

        if self.dostress:
            self.stress.write('(kbar)\n')
            for v in at.virial:
                self.stress.write('%20.10e%20.10e%20.10e\n' % tuple(v*(10.0*GPA)/at.cell_volume()))

        self.it += 1

    def close(self):
        if self.pos != sys.stdout: self.pos.close()
        if self.cel != sys.stdout: self.cel.close()
        if self.doenergy and self.energy != sys.stdout: self.energy.close()
        if self.doforce  and self.force  != sys.stdout: self.force.close()
        if self.dostress and self.stress != sys.stdout: self.stress.close()


AtomsWriters['pos'] = PosCelWriter

@atoms_reader('pos')
def PosCelReader(basename=None, pos='pos.in', cel='cel.in', force='force.in', energy='energy.in', stress='stress.in',
                 species_map={'O':1, 'Si':2}, cel_angstrom=False, pos_angstrom=False, rydberg=True, format=None):

    if basename is not None:
        basename = os.path.splitext(basename)[0]
        pos = '%s.pos' % basename
        cel = '%s.cel' % basename
        energy = '%s.ene' % basename
        stress = '%s.str' % basename
        force = '%s.for' % basename

    doenergy = os.path.exists(energy)
    doforce  = os.path.exists(force)
    dostress = os.path.exists(stress)

    if isinstance(pos, str): pos = open(pos)
    if isinstance(cel, str): cel = open(cel)
    if doenergy and isinstance(energy, str): energy = open(energy)
    if doforce  and isinstance(force,  str): force  = open(force)
    if dostress and isinstance(stress, str): stress = open(stress)

    pos = iter(pos)
    cel = iter(cel)
    if doenergy: energy = iter(energy)
    if doforce:  force  = iter(force)
    if dostress: stress = iter(stress)

    next(pos)    # throw away blank line at start
    if doforce:
        next(force)

    rev_species_map = dict(list(zip(list(species_map.values()), list(species_map.keys()))))

    while True:

        poslines = list(itertools.takewhile(lambda L: L.strip() != '' and not L.strip().startswith('STEP'), pos))
        if poslines == []:
            break

        cellines = list(itertools.islice(cel, 4))
        #lattice = farray([ [float(x) for x in L.split()] for L in cellines[1:4] ]).T
        lattice = fzeros((3,3))
        for i in (1,2,3):
            lattice[:,i] = [float(x) for x in cellines[i].split()]
        if not cel_angstrom: lattice *= BOHR

        at = Atoms(n=len(poslines), lattice=lattice)
        at.pos[:] = farray([ [float(x) for x in L.split()[0:3] ] for L in poslines ]).T
        if not pos_angstrom: at.pos[:] *= BOHR
        species = [ rev_species_map[int(L.split()[3])] for L in poslines ]
        elements = [ not el.isdigit() and atomic_number(el) or el for el in species ]
        at.set_atoms(elements)

        if doenergy:
            at.params['energy'] = float(energy.next().split()[0])
            if rydberg:
                at.params['energy'] *= RYDBERG

        if dostress:
            stress_lines = list(itertools.islice(stress, 4))
            virial = farray([ [float(x) for x in L.split()] for L in stress_lines[1:4] ])
            virial *= at.cell_volume()/(10.0*GPA)
            at.params['virial'] = virial

        if doforce:
            at.add_property('force', 0.0, n_cols=3)
            force_lines = list(itertools.takewhile(lambda L: L.strip() != '', force))
            if len(force_lines) != at.n:
                raise ValueError("len(force_lines) (%d) != at.n (%d)" % (len(force_lines), at.n))
            at.force[:] = farray([ [float(x) for x in L.split()[0:3] ] for L in force_lines ]).T
            if rydberg:
                at.force[:] *= RYDBERG/BOHR

        yield at

logical = lambda x:  x in (1, True, '1', 'True', 't', 'T', '.t.', '.T')
real = lambda s: float(s.replace('d', 'e'))

param_format_traj = [
   [('mass', real, 'nspecies', False)],
   [('species', str, 'nspecies', False)],
   [('nsp', int, 'nspecies', False)],
   [('tewald', logical, 1, False),('raggio', real, 1, False),('a_ew', real, 1, False),('gcut', real, 1, False),('iesr', int, 3, False),('rcut', real, 4, False)],
   [('z', real, 'nspecies-1', False)],
   [('alphaij', real, 'triangle(nspecies)', False)],
   [('bij', real, 'triangle(nspecies)', False)],
   [('cij', real, 'triangle(nspecies)', False)],
   [('dij', real, 'triangle(nspecies)', False)],
   [('eij', real, 'triangle(nspecies)', False)],
   [('nij', real, 'triangle(nspecies)', False)],
   [('b_tt1', real, 'triangle(nspecies)', False)],
   [('b_tt2', real, 'triangle(nspecies)', False)],
   [('b_tt3', real, 'triangle(nspecies)', False)],
   [('d_ms', real, 'triangle(nspecies)', False)],
   [('gamma_ms', real, 'triangle(nspecies)', False)],
   [('r_ms', real, 'triangle(nspecies)', False)],
   [('pol', real, 'nspecies', False)],
   [('betapol', real, 1, False),('maxipol', int, 1, False),('tolpol', real, 1, False), ('pred_order', int, 1, False), ('tdip_sr', logical, 1, False)],
   [('bpol', real, 'triangle(nspecies)', False)],
   [('cpol', real, 'triangle(nspecies)', False)],
   [('taimsp', logical, 'nspecies', False),('xgmin', real, 1, False),('xgmax', real, 1, False)],
   [('adist',   real, 'triangle(nspecies)', False)],
   [('bdist',   real, 'triangle(nspecies)', False)],
   [('cdist',   real, 'nspecies', False)],
   [('ddist',   real, 'triangle(nspecies)', False)],
   [('sigma1',  real, 'nspecies', False)],
   [('sigma2',  real, 'nspecies', False)],
   [('sigma3',  real, 'nspecies', False)],
   [('sigma4',  real, 'nspecies', False)],
   [('bu1',     real, 'triangle(nspecies)', False)],
   [('alphau1', real, 'triangle(nspecies)', False)],
   [('bu2',     real, 'triangle(nspecies)', False)],
   [('alphau2', real, 'triangle(nspecies)', False)],
   [('c_harm',  real, 'triangle(nspecies)', False)],
   [('rmin', real, 'triangle(nspecies)', False)],
   [('a_s', real, 1, False),('n_s', real, 1, False),('smooth', logical, 1, False)],
   [('tyukawa', logical, 1, False), ('yukalpha', real, 1, False),('yuksmoothlength', real, 1, False)]
   ]


param_format_gen = [
   [('testewald', logical, 1, False), ('time', logical, 1, False), ('tpbc', logical, 1, False), ('tangstrom', logical, 1, False),
    ('trydberg', logical, 1, False), ('tev', logical, 1, False)],
   [('nat', int, 1, False), ('nspecies', int, 1, False), ('npar', int, 1, False)],
   [('ntmin', int, 1, False)],
   [('nts', int ,'ntmin', False)],
   [('filepos', str, 1, False)],
   [('fileforce', str, 1, False)],
   [('filestress', str, 1, False)],
   [('filecel', str, 1, False)],
   [('fileene', str, 1, False)],
   [('testforce', logical, 1, False), ('ngrid', int, 1, False), ('gridmin1', real, 1, False), ('gridmax1', real, 1, False),
    ('gridmin2', real, 1, False), ('gridmax2', real, 1, False), ('isp_tf', int, 1, False)],
   [('tquickpar', logical, 1, False)],
   [('mass', real, 'nspecies', False)],
   [('nsp', int, 'nspecies', False)],
   [('tewald', logical, 1, False),('raggio', real, 1, False),('a_ew', real, 1, False),('gcut', real, 1, False),('iesr', int, 3, False),('rcut', real, 4, False)],
   [('z', real, 'nspecies-1', False), ('tz', logical, 'nspecies-1', False)],
   [('alphaij', real, 'triangle(nspecies)', True)],
   [('bij', real, 'triangle(nspecies)', True)],
   [('cij', real, 'triangle(nspecies)', True)],
   [('dij', real, 'triangle(nspecies)', True)],
   [('eij', real, 'triangle(nspecies)', True)],
   [('nij', real, 'triangle(nspecies)', True)],
   [('b_tt1', real, 'triangle(nspecies)', True)],
   [('b_tt2', real, 'triangle(nspecies)', True)],
   [('b_tt3', real, 'triangle(nspecies)', True)],
   [('d_ms', real, 'triangle(nspecies)', True)],
   [('gamma_ms', real, 'triangle(nspecies)', True)],
   [('r_ms', real, 'triangle(nspecies)', True)],
   [('pol', real, 'nspecies', True)],
   [('betapol', real, 1, False),('maxipol', int, 1, False),('tolpol', real, 1, False), ('tdip_sr', logical, 1, False)],
   [('bpol', real, 'triangle(nspecies)', True)],
   [('cpol', real, 'triangle(nspecies)', True)],
   [('taimsp', logical, 'nspecies', False),('xgmin', real, 1, False),('xgmax', real, 1, False)],
   [('adist',   real, 'triangle(nspecies)', True)],
   [('bdist',   real, 'triangle(nspecies)', True)],
   [('cdist',   real, 'nspecies', True)],
   [('ddist',   real, 'triangle(nspecies)', True)],
   [('sigma1',  real, 'nspecies', True)],
   [('sigma2',  real, 'nspecies', True)],
   [('sigma3',  real, 'nspecies', True)],
   [('sigma4',  real, 'nspecies', True)],
   [('bu1',     real, 'triangle(nspecies)', True)],
   [('alphau1', real, 'triangle(nspecies)', True)],
   [('bu2',     real, 'triangle(nspecies)', True)],
   [('alphau2', real, 'triangle(nspecies)', True)],
   [('c_harm',  real, 'triangle(nspecies)', True)],
   [('rmin', real, 'triangle(nspecies)', False)],
   [('a_s', real, 1, False),('n_s', real, 1, False),('smooth', logical, 1, False)],
   [('tyukawa', logical, 1, False), ('yukalpha', real, 1, False),('yuksmoothlength', real, 1, False)]
   ]

def triangle(n):
    """Return list of length n which sums to nth triangular number"""
    res = []
    for i in range(1,n+1):
        res.append(i)
    res.reverse()
    return res

def inv_trig(t):
    """Given the nth triangular number t_n, return n by inverting t_n = 1/2 n (n+1)"""
    return int(0.5*(((1+8*t)**0.5)-1))

output_converters = {
   (real, 'triangle(nspecies)'): lambda v: '\n'.join(['%E '*n for n in triangle(inv_trig(len(v)))]) % tuple(v),
   int: lambda v: ' '.join(['%d' % x for x in v]),
   real: lambda v: ' '.join(['%f' % x for x in v]),
   logical: lambda v: ' '.join([x and '.t.' or '.f.' for x in v]),
   str: lambda v: ' '.join(v)
   }

type_map = {}
for line in param_format_gen:
    for key, conv, nfields, interleave in line:
        if (key, conv, nfields) in output_converters:
            invconv = output_converters[(key, conv, nfields)]
        elif (conv, nfields) in output_converters:
            invconv = output_converters[(conv, nfields)]
        elif conv in output_converters:
            invconv = output_converters[conv]
        else:
            invconv = lambda s: str(s)
        type_map[key] = (conv, invconv)

def read_traj_gen(format, gen_file, defaults={}):

    lines = list(gen_file)
    lines = [x for x in lines if not (x.strip().startswith('--')
                                      or x.strip().startswith('**')
                                      or x.strip().startswith('%%'))]

    lengths = [ len(x) for x in format ]

    origlines = lines[:]

    params = defaults.copy()
    fixed = {}
    fixed_order = []

    evaldict = defaults.copy()
    evaldict['triangle'] = triangle
    for entries in format:
        n = 0
        gotline = False
        for key, conv, nfields, interleave in entries:
            if isinstance(nfields, str):
                nfields = eval(nfields, evaldict)


            if isinstance(nfields, list) and len(nfields) != 1:
                # multiple lines

                if interleave:
                    mylines = lines[:2*len(nfields):2]
                    ilines = lines[1:2*len(nfields):2]
                    del lines[:2*len(nfields)]
                else:
                    mylines = lines[:len(nfields)]
                    ilines = ['']*len(mylines)
                    del lines[:len(nfields)]

                for line, iline, nf in zip(mylines, ilines, nfields):
                    fields = line.split()
                    values = [conv(x) for x in fields[:nf]]
                    if key in params:
                        params[key] = params[key] + values
                        evaldict[key] = params[key] + values
                    else:
                        if len(values) == 1: values = values[0]
                        params[key] = values
                        evaldict[key] = values

                    if interleave:
                        tfields = iline.split()
                        tvalues = [logical(x) for x in tfields[:nf]]

                        if key in fixed:
                            fixed[key] += tvalues
                        else:
                            fixed_order.append(key)
                            fixed[key] = tvalues

            else:
                # just one line, possibly multiple values per line

                if not gotline:
                    line = lines[0]
                    del lines[0]
                    if interleave:
                        iline = lines[0]
                        del lines[0]
                    gotline = True

                fields = line.split()
                values = [conv(x) for x in fields[n:n+nfields]]
                if len(values) == 1: values = values[0]
                params[key] = values
                evaldict[key] = values

                if interleave:
                    tvalues = [logical(x) for x in iline.split()[n:n+nfields]]

                    if key in fixed:
                        fixed[key] += tvalues
                    else:
                        fixed_order.append(key)
                        fixed[key] = tvalues

                n += nfields

    # special case: tz logical option is on same line as z
    if 'tz' in params:
        fixed['z'] = [params['tz']]
        fixed_order.insert(0,'z')
        del params['tz']

    opt_vars = []
    for var in fixed_order:
        if len(fixed[var]) == 1:
            if fixed[var]: opt_vars.append(var)
        else:
            for i,v in enumerate(fixed[var]):
                if v: opt_vars.append((var, i))

    return params, opt_vars

traj_format_str = """%(mass)s
%(species)s
%(nsp)s
%(tewald)s %(raggio)s %(a_ew)s %(gcut)s %(iesr)s %(rcut)s
%(z)s
-------------Alphaij---------------
%(alphaij)s
-------------Bij--------------------
%(bij)s
---------------Cij------------------
%(cij)s
-----------------Dij----------------
%(dij)s
---------------Eij------------------
%(eij)s
---------------Nij------------------
%(nij)s
---------Tang-Toennies-------------
%(b_tt1)s
---------Tang-Toennies-------------
%(b_tt2)s
---------Tang-Toennies-------------
%(b_tt3)s
---------------D_ms----------------
%(d_ms)s
---------------Gamma_ms------------
%(gamma_ms)s
----------------R_ms---------------
%(r_ms)s
--------------Polarization---------
%(pol)s
%(betapol)s %(maxipol)s %(tolpol)s %(pred_order)s
---------------Bpol----------------
%(bpol)s
---------------Cpol----------------
%(cpol)s
--------Aspherical-Ion-Model-------
%(taimsp)s %(xgmin)s %(xgmax)s
-----------Adist-------------------
%(adist)s
---------------Bdist---------------
%(bdist)s
---------------Cdist---------------
%(cdist)s
--------------Ddist----------------
%(ddist)s
-------------Sigma1----------------
%(sigma1)s
-------------Sigma2----------------
%(sigma2)s
------------Sigma3-----------------
%(sigma3)s
-------------Sigma4----------------
%(sigma4)s
-------------Bu1-------------------
%(bu1)s
--------------Alphau1--------------
%(alphau1)s
---------------Bu2-----------------
%(bu2)s
--------------Alphau2--------------
%(alphau2)s
***********Spring Constant*********
%(c_harm)s
***********Spring Cutoff***********
%(rmin)s
**********Smoothing***************
%(a_s)s %(n_s)s %(smooth)s
*************Yukawa***************
%(yukalpha)s %(yuksmoothlength)s %(tdip_sr)s"""

def param_to_str(format, params):

    out_params = {}

    for key in params:
        try:
            len(params[key])
            out_params[key] = type_map[key][1](params[key])
        except TypeError:
            out_params[key] = type_map[key][1]([params[key]])

    return format % out_params


def param_to_xml(params, encoding='iso-8859-1'):
    from xml.sax.saxutils import XMLGenerator
    from io import StringIO
    output = StringIO()
    xml = XMLGenerator(output, encoding)
    xml.startDocument()

    xml.startElement('TS_params', {'cutoff': ' '.join([str(x) for x in params['rcut']]),
                                     'n_types': str(params['nspecies']),
                                     'betapol': str(params['betapol']),
                                     'maxipol': str(params['maxipol']),
                                     'tolpol': str(params['tolpol']),
                                     'pred_order': str(params['pred_order']),
                                     'yukalpha': str(params['yukalpha']),
                                     'yuksmoothlength': str(params['yuksmoothlength']),
                                     'tewald': params['tewald'] and 'T' or 'F',
                                     'raggio': str(params['raggio']),
                                     'a_ew': str(params['a_ew']),
                                     'gcut': str(params['gcut']),
                                     'iesr': ' '.join([str(x) for x in params.get('iesr', [0,0,0])])
                                     })
    ti_tj_to_index = {}
    n = 0
    for ti in range(params['nspecies']):
        for tj in range(params['nspecies']):
            if tj > ti: continue
            ti_tj_to_index[(ti,tj)] = n
            n += 1

    for ti in range(params['nspecies']):
        zi = atomic_number(params['species'][ti])
        xml.startElement('per_type_data', {'type':str(ti+1),
                                           'atomic_num': str(zi),
                                           'pol': str(params['pol'][ti]),
                                           'z': str(params['z'][ti])})
        xml.endElement('per_type_data')

        for tj in range(params['nspecies']):
            if tj > ti: continue
            idx = ti_tj_to_index[(ti,tj)]
            zj = atomic_number(params['species'][tj])
            xml.startElement('per_pair_data', {'atnum_i':str(zi),
                                               'atnum_j':str(zj),
                                               'D_ms': str(params['d_ms'][idx]),
                                               'gamma_ms': str(params['gamma_ms'][idx]),
                                               'R_ms': str(params['r_ms'][idx]),
                                               'B_pol': str(params['bpol'][idx]),
                                               'C_pol': str(params['cpol'][idx]),
                                               })
            xml.endElement('per_pair_data')

    xml.endElement('TS_params')
    xml.endDocument()
    return output.getvalue()

def update_vars(params, opt_vars, param_str, verbose=False):
    out_params = params.copy()
    if isinstance(param_str, str):
        opt_values = [real(x) for x in param_str.split()]
    else:
        opt_values = param_str

    assert(len(opt_values) == len(opt_vars))

    for var, value in zip(opt_vars, opt_values):
        if isinstance(var, str):
            if var not in out_params: raise ValueError('var %s missing' % var)
            if verbose: print('%s   %f -> %f' % (var, params[var], value))
            out_params[var] = value
        else:
            key, idx = var
            if key not in out_params: raise ValueError('var %s missing' % key)
            if verbose: print('%s[%d] %f -> %f' % (key, idx, params[key][idx], value))
            out_params[key][idx] = value

    return out_params



default_params = {'species': ['O', 'Si'],
                  'pred_order': 2}

type_map['species'] = (str, output_converters[str])
type_map['pred_order'] = (int, output_converters[int])


def save_ref_config(config_list):

    for at in config_list:
        at.params['dft_energy'] = at.params['energy'] / HARTREE
        at.add_property('dft_force', 0.0, n_cols=3)
        at.dft_force[:] = at.force / (HARTREE/BOHR)
        at.params['dft_stress'] = at.virial/at.cell_volume() / (HARTREE/(BOHR**3))


def costfn(config_list, pot, wf=1.0, ws=0.5, we=0.1, bulk_mod=2000.0/294156.6447):

    normweight = (wf*wf + ws*ws + we*we)**0.5
    wf /= normweight
    ws /= normweight
    we /= normweight

    dist_f = 0.0; norm_f = 0.0
    dist_e = 0.0; norm_e = 0.0
    dist_s = 0.0; norm_s = 0.0

    for ati, at in enumerate(config_list):
        pot.calc(at, virial=True, force=True, energy=True)

        at.params['md_energy'] = at.params['energy'] / HARTREE
        at.params['md_stress'] = at.virial/at.cell_volume() / (HARTREE/(BOHR**3))
        at.md_force = at.force / (HARTREE/BOHR)

        for i in frange(at.n):
            for j in (1,2,3):
                dist_f += (at.md_force[j,i] - at.dft_force[j,i])**2.0
                norm_f += at.dft_force[j,i]**2.0

        for i in (1,2,3):
            for j in (1,2,3):
                dist_s += (at.md_stress[i,j] - at.dft_stress[i,j])**2.0
                norm_s += bulk_mod**2.0


    for i, ati in enumerate(config_list):
        for j, atj in enumerate(config_list):
            if j <= i: continue

            diff_md = ati.md_energy - atj.md_energy
            diff_ai = ati.dft_energy - atj.dft_energy
            norm_e += diff_ai**2.0
            dist_e += (diff_md - diff_ai)**2.0

    dist_f = (dist_f / norm_f)**0.5
    dist_s = (dist_s / norm_s)**0.5
    dist_e = (dist_e / norm_e)**0.5

    return config_list, (wf*dist_f + ws*dist_s + we*dist_e, dist_f, dist_s, dist_e)

def read_gen_and_param_files(gen_file, param_file, verbose=True):

    if isinstance(gen_file, str):
        gen_file = open(gen_file, 'r')

    params = default_params.copy()
    gen_params, opt_vars = read_traj_gen(param_format_gen, gen_file)
    params.update(gen_params)

    if verbose: print('opt_vars = ', opt_vars)

    params = update_vars(params, opt_vars, param_file, verbose=verbose)

    # fix masses
    params['mass'] = [ElementMass[x]/MASSCONVERT for x in params['species']]

    # fix charges
    params['z'] = [params['z'], -2.0*params['z']]

    return params


def read_minimisation_progress(file):
    if isinstance(file, str):
        file = open(file)

    params = []
    costs = []
    while True:
        try:
            line = next(file) # skip blank line
            if not line:
                break
            costs.append([float(x) for x in next(file).split()[1:]])
            params.append([float(x) for x in (next(file) + next(file)).split()])
        except StopIteration:
            break

    return farray(costs), farray(params)

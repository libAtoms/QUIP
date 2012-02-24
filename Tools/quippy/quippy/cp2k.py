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

import sys
import os
import numpy as np
from quippy.system import verbosity_push, verbosity_pop, PRINT_SILENT
from quippy.atoms import Atoms, atoms_reader
from quippy.farray import fzeros, s2a, farray, fenumerate
from quippy.table import TABLE_STRING_LENGTH
from quippy.units import MASSCONVERT, HARTREE, BOHR
from quippy.clusters import HYBRID_NO_MARK
from quippy.potential import Potential
from quippy.ordereddict import OrderedDict
from quippy.util import read_text_file
from quippy.cp2k_driver import do_cp2k_calc, read_output, qmmm_qm_abc

__all__ = ['CP2KPotential', 'CP2KInputHeader']

class CP2KPotential(Potential):
    """
    QUIP Potential interface to the CP2K code.

    Calls do_cp2k_calc() in QUIP_FilePot_Drivers/cp2k_driver_module.f95 to do the heavy lifting
    """
    def __init__(self, fortran_indexing=True,
                 fpointer=None, finalise=True,
                 error=None):
        Potential.__init__(self, callback=self.do_cp2k_calc)

    def do_cp2k_calc(self, at):
        args_str = at.params.get('calc_args_str', '')
        cp2k_force, cp2k_energy = do_cp2k_calc(at, infile='', args_str=args_str, n0=3, n1=at.n)
        at.add_property('force', cp2k_force, overwrite=True)
        at.params['energy'] = cp2k_energy


class CP2KInputHeader(OrderedDict):
    """
    Read a cp2k_input.inp.header file and parse into an ordered dictionary.
    """

    def __init__(self, filename=None):
        OrderedDict.__init__(self)
        if filename is not None:
            self.read(filename)

    def read(self, filename):
        filename, lines = read_text_file(filename)
        for line in lines:
            dummy, key, value = line.split()
            self[key] = value

    def write(self, filename=None):
        if filename is None:
            filename = sys.stdout
        if isinstance(filename, basestring):
            file = open(filename, 'w')
        else:
            file = filename
    
        for (key, value) in self.iteritems():
            file.write('@SET %s %s\n' % (key, value))
            
        if isinstance(filename, basestring):            
            file.close()


def cp2k_run_type(cp2k_output=None, cp2k_input_header=None):
    """
    Find the run type of the CP2K calculation.

    Either `cp2k_output` or `cp2k_input_header` should be present.
    When both are givenm we check output is consistent with input,
    raising a ValueError if not.

    Returns one of 'QS', 'MM' or 'QMMM'.
    """

    if cp2k_output is None and cp2k_input_header is None:
        raise ValueError("atlest one of cp2k_output and cp2k_input_header must be present")

    if cp2k_output is not None:
        got_fist = " MODULE FIST:  ATOMIC COORDINATES IN angstrom\n" in cp2k_output
        got_qs   = " MODULE QUICKSTEP:  ATOMIC COORDINATES IN angstrom\n" in cp2k_output
        got_qmmm = " MODULE QM/MM:  ATOMIC COORDINATES IN angstrom\n" in cp2k_output
        got_modules = (got_fist, got_qs, got_qmmm)

    if cp2k_input_header is not None:
        do_fist = int(cp2k_input_header['DO_MM'])
        do_qs = int(cp2k_input_header['DO_DFT'])
        do_qmmm = int(cp2k_input_header['DO_QMMM'])
        got_modules = (do_fist, do_qs, do_qmmm)

    if cp2k_output is not None and cp2k_input_header is not None:
        if (do_fist != got_fist or do_qs != got_qs or do_qmmm != got_qmmm):
            raise ValueError("CP2K output inconsitent with cp2k input header")

    if got_modules == (True, True, True):
        return "QMMM"
    elif got_modules == (True, False, False):
        return "MM"
    elif got_modules == (False, True, False):
        return "QS"
    else:
        raise ValueError("Got inconsistent set of modules got_fist=%s, got_qs=%s, got_qmmm+%s" % got_modules)
    

@atoms_reader('cp2k_output.out')
@atoms_reader('cp2k_output')
def CP2KOutputReader(fh, module=None, type_map=None, kind_map=None):

    # mapping from run type to (default module index, list of available module)
    run_types = {
        'QS':   ['QUICKSTEP'],
        'QMMM': ['FIST', 'QM/MM', 'QUICKSTEP'],
        'MM' :  ['FIST']
        }

    filename, lines = read_text_file(fh)
    run_type = cp2k_run_type(cp2k_output=lines)
    
    if type_map is None:
        type_map = {}
    if kind_map is None:
        kind_map = {}

    try:
        available_modules = run_types[run_type]
    except KeyError:
        raise ValueError('Unknown CP2K run type %s' % run_type)

    if module is None:
        module = available_modules[0]

    try:
        cell_index = available_modules.index(module)
    except IndexError:
        raise ValueError("Don't know how to read module %s from file %s" % (module, filename))

    cell_lines = [i for i,line in enumerate(lines) if line.startswith(" CELL| Vector a")]
    if cell_lines == []:
        raise ValueError("Cannot find cell in file %s" % filename)

    try:
        cell_line = cell_lines[cell_index]
    except IndexError:
        raise ValueError("Cannot find cell with index %d in file %s for module %s" %
                         (cell_index, filename, module))

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


@atoms_reader('cp2k_run_dir')
def CP2KDirectoryReader(run_dir, at_ref=None, proj='quip', calc_qm_charges=None,
                        calc_virial=False, out_i=None, qm_vacuum=6.0, run_suffix='_extended'):
    if at_ref is None:
        filepot_xyz = os.path.join(run_dir, 'filepot.xyz')
        if not os.path.exists(filepot_xyz):
            # try looking up one level 
            filepot_xyz = os.path.join(run_dir, '../filepot.xyz')

        if os.path.exists(filepot_xyz):
            at_ref = Atoms(filepot_xyz)
        else:
            at_ref = Atoms(os.path.join(run_dir, 'cp2k_output.out'),
                           format='cp2k_output')

    at = at_ref.copy()

    cp2k_output_filename, cp2k_output = read_text_file(os.path.join(run_dir, 'cp2k_output.out'))
    cp2k_params = CP2KInputHeader(os.path.join(run_dir, 'cp2k_input.inp.header'))
    at.params.update(cp2k_params)

    run_type = cp2k_run_type(cp2k_output=cp2k_output, cp2k_input_header=cp2k_params)

    try:
        cluster_mark = getattr(at, 'cluster_mark'+run_suffix)
        qm_list_a = ((cluster_mark != HYBRID_NO_MARK).nonzero()[0]).astype(np.int32)
    except AttributeError:
        qm_list_a = fzeros(0, dtype=np.int32)

    if calc_qm_charges is None:
        calc_qm_charges = ''

    try:
        cur_qmmm_qm_abc = [float(cp2k_params['QMMM_ABC_X']),
                           float(cp2k_params['QMMM_ABC_Y']),
                           float(cp2k_params['QMMM_ABC_Z'])]
    except KeyError:
        if 'QM_cell'+run_suffix in at.params:
            cur_qmmm_qm_abc = at.params['QM_cell'+run_suffix]
        else:
            cur_qmmm_qm_abc = qmmm_qm_abc(at, qm_list_a, qm_vacuum)

    verbosity_push(PRINT_SILENT)
    cp2k_energy, cp2k_force = read_output(at, qm_list_a, cur_qmmm_qm_abc, run_dir, proj,
                                          calc_qm_charges, calc_virial, True, 3, at.n, out_i)
    verbosity_pop()

    qm_list = None
    if os.path.exists(os.path.join(run_dir, 'cp2k_input.qmmm_qm_kind')):
        qm_kind_grep_cmd = "grep MM_INDEX %s/cp2k_input.qmmm_qm_kind | awk '{print $2}'" % run_dir
        qm_list = [int(i) for i in os.popen(qm_kind_grep_cmd).read().split()]

    if qm_list is not None:
        if run_type == 'QMMM':
            reordering_index = getattr(at, 'reordering_index', None)

            at.add_property('qm', False, overwrite=True)
            if reordering_index is not None:
                qm_list = reordering_index[qm_list]
            at.qm[qm_list] = True
        elif run_type == 'QS':

            # check for a reverse sort index, and apply to qm_list if found
            rev_sort_index_file = os.path.join(run_dir, 'quip_rev_sort_index')
            if os.path.exists(rev_sort_index_file):
                rev_sort_index = farray(np.loadtxt(rev_sort_index_file))
                rev_sort_index = rev_sort_index.reshape(rev_sort_index.size).astype(int)
                sort_index = rev_sort_index.argsort()
            else:
                sort_index = farray(np.arange(1, at.n+1))
            
            at.add_property('qm_orig_index', 0, overwrite=True)
            for i, qm_at in fenumerate(qm_list):
                at.qm_orig_index[i] = sort_index[qm_at]

            
        

    at.add_property('force', cp2k_force, overwrite=True)
    at.params['energy'] = cp2k_energy
    yield at



def read_cp2k_qm_kind(fh):
    """
    Read a file in cp2k_input.qmmm_qm_kind

    This is the format produced by the QUIP cp2k_driver.
    Returns a list of tuples (qm_kind, atom_indices).
    """
    
    filename, lines = read_text_file(fh)
    qm_kinds = []
    current_qm_kind = None
    current_indices = []
    for line in lines:
        line = line.strip()
        if line.startswith('&QM_KIND'):
            current_qm_kind = line.split()[1]
            current_indices = []
        elif line.startswith('&END QM_KIND'):
            qm_kinds.append((current_qm_kind,current_indices))
        elif line.startswith('MM_INDEX'):
            index = int(line.split()[1])
            current_indices.append(index)
        else:
            raise ValueError('Unexpected line %r' % line)

    return qm_kinds
    

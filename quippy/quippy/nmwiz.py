# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy writer for normal mode (`.nmd`) files
# HQ X   Copyright Max Veit 2018
# HQ X   File format used by ProDy (http://prody.csb.pitt.edu/)
# HQ X   and VMD's NMWiz plugin
# HQ X   (http://www.ks.uiuc.edu/Research/vmd/plugins/nmwiz/)
# HQ X
# HQ X   Part of:
# HQ X   quippy: Python interface to QUIP atomistic simulation library
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
"""nmwiz: Write atoms and normal modes to '.nmd' format.

Contents:
    NMDWriter   Writer for the .nmd format
"""


from __future__ import division, print_function, unicode_literals
from itertools import izip
import sys

from numpy import abs, sqrt
import quippy


class NMDWriter(object):

    """Writer for the .nmd format

    The Atoms object must have normal modes -- both eigenvectors and
    eigenvalues (force constants) of the Hessian matrix -- stored in
    ``atoms.info['hessians_forceconst_X'`` and
    ``atoms.arrays['hessians_vector_X']``, where ``X`` ranges from 1
    to the number of stored normal modes.

    This writer implements the context manager protocol, so you can
    use it like so::

        with NMDWriter(filename) as writer:
            writer.write(atoms)

    and the associated file will be closed automatically when done.

    The ``.nmd`` file format only supports single configurations, so
    this writer doesn't accept trajectories.
    """

    def __init__(self, filename):
        if filename == 'stdout':
            self._file = sys.stdout
        else:
            self._file = open(filename, 'w')

    def __enter__(self):
        # The writer interface requires we open the file in __init__(),
        # so do nothing here.
        # In particular, don't reopen the file if it's already been closed.
        return self

    def write(self, atoms, title='quippy atoms'):
        """Write out the atoms and modes with an optional title"""
        # Get the modes from the Atoms object
        try:
            mode_idx = 1
            eigvecs = []
            mode_idces = []
            while True:
                mode_string = 'hessians_vector_{:d}'.format(mode_idx)
                eigvecs.append(atoms.get_array(mode_string))
                mode_idces.append(mode_idx)
                mode_idx += 1
        except KeyError as kerr:
            if not mode_idces:
                # Oh great, Python 2 doesn't support exception chaining
                #py2kfacepalm
                #raise ValueError("Couldn't find any mode vectors "
                #                 "in Atoms object") from kerr
                raise ValueError(
                    "Couldn't find any mode vectors in Atoms object "
                    "(key '{:s}' missing)".format(mode_string))
        try:
            eigvals = []
            for mode_idx in mode_idces:
                mode_string = 'hessians_forceconst_{:d}'.format(mode_idx)
                eigvals.append(atoms.info[mode_string])
        except KeyError as kerr:
            # This is how easy it would be in Python 3
            #raise ValueError("Couldn't find force constant for mode number "
            #                 "{:d}".format(mode_idx)) from kerr
            raise ValueError("Couldn't find force constant for mode number "
                             "{:d} (key '{:s}')".format(mode_idx, mode_string))
        # Write the actual file
        self._file.write('title ' + title + '\n')
        self._file.write('coordinates ')
        atoms.get_positions().tofile(self._file, ' ', '%.6f')
        self._file.write('\n')
        self._file.write('names ')
        self._file.write(' '.join(atoms.get_chemical_symbols()))
        self._file.write('\n')
        for idx, eigval, eigvec in izip(mode_idces, eigvals, eigvecs):
            self._file.write('mode {:d} {:.6f} '.format(idx, eigval))
            eigvec.tofile(self._file, ' ', '%.6f')
            self._file.write('\n')

    def __exit__(self, exc_type, exc_value, traceback):
        if self._file is not sys.stdout:
            self._file.close()
        # Don't suppress any exceptions
        return False

    def close(self):
        self.__exit__(None, None, None)


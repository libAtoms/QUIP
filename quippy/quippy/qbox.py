"""Utilities for writing Qbox input files

Public classes:
    QboxInputWriter         Write Qbox input files from Atoms objects

Public functions:
    format_filename_seq     Format a filename with an index

Only writing is implemented here, as the Qbox distribution includes a
utility for converting output files to XYZ (without cell info or
energies, though).

This code is copyright (c) 2016 Max Veit and released under the GNU
General Public License, version 2.

See qboxcode.org for more about Qbox.
"""


import collections
import errno
import os
import shutil

import numpy as np
from quippy.periodictable import ElementName
from quippy.units import BOHR, AU_FS
from quippy.io import AtomsWriters


def format_filename_seq(name, idx):
    """Format a filename with an index (placed before any suffix)"""
    path_split = os.path.splitext(name)
    return path_split[0] + '_{:d}'.format(idx) + path_split[1]


class QboxInputWriter(object):

    """Writer that makes Qbox input files from quippy.Atoms objects

    Since each input file can only contain one geometry (Atoms object),
    files are opened and closed on each write() call, rather than
    opening on initialization and closing on the close() call.

    Writes AtomsLists as a sequence of multiple files: The index is
    inserted into the output filename before any suffix.  This mode is
    triggered by calling write() more than once on a single object.  On
    the second write() call, the first file written is moved to
    <filename>_0[.<suffix>] and subsequent files are written using the
    same format.

    NOTE: The files written do not contain any Qbox species definition.
    The species are named by the corresponding entry in
    quippy.periodictable.ElementName.  The final qbox input needs to
    contain species definitions, such as 'species H <species_def.xml>'
    for hydrogen.
    """

    def __init__(self, filename):
        """Create a writer for the given filename."""
        self.out_fname = filename
        self._traj_idx = 0
        self._wrote_single = True
        self._closed = False
        fname_head = os.path.split(filename)[0]
        if fname_head:
            try:
                os.makedirs(fname_head)
            except OSError as ose:
                if ose.errno != errno.EEXIST:
                    raise ose

    def write(self, atoms):
        """Write the Atoms object to a Qbox input file.

        See the class documentation for specifics on writing
        trajectories, filename manipulation, and file formatting.
        """
        if self._closed:
            raise RuntimeError("Attempting to write using a closed "
                               "QboxInputWriter")
        if self._traj_idx == 0:
            filename_cur = self.out_fname
        else:
            filename_cur = format_filename_seq(self.out_fname, self._traj_idx)
        if self._traj_idx > 0 and self._wrote_single:
            shutil.move(self.out_fname, format_filename_seq(self.out_fname, 0))
            self._wrote_single = False
        with open(filename_cur, 'w') as outf:
            cell_bohr = atoms.cell / BOHR
            cell_param_str = ' '.join(str(num) for num in cell_bohr.flat)
            outf.write('set cell ' + cell_param_str + '\n')
            species_count = collections.defaultdict(lambda: 1)
            atom_fmt = 'atom\t{:s}\t{:s}' + ('\t{:.8f}' * 3)
            atom_props = atoms.pos / BOHR
            if atoms.has_property('velo'):
                atom_fmt += ('\t{:.8f}' * 3)
                atom_props = np.concatenate(
                    (atoms.pos / BOHR,
                     atoms.velo * AU_FS / BOHR), axis=0)
            for at_idx in atoms.indices:
                sp_name = ElementName[atoms.Z[at_idx]]
                at_name = sp_name + str(species_count[atoms.Z[at_idx]])
                species_count[atoms.Z[at_idx]] += 1
                atom_line = atom_fmt.format(
                    at_name, sp_name, *atom_props[:, at_idx])
                outf.write(atom_line + '\n')
        self._traj_idx += 1

    def close(self):
        """Close the writer."""
        del self.out_fname
        self._closed = True


AtomsWriters['qbox'] = QboxInputWriter


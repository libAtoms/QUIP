# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2019
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

import unittest
import os
import xml.etree.ElementTree as ET
import json
import subprocess
import string

import numpy as np
import quippytest

from quippy.potential import Potential
from ase.constraints import voigt_6_to_full_3x3_stress
from ase.io import read, write

@unittest.skipIf(os.environ['HAVE_GAP'] != '1', 'GAP support not enabled')
class TestGAP_fit(quippytest.QuippyTestCase):
    alpha_tol = 1e-5

    def setUp(self):
        self.cl_template = string.Template(
            "at_file=train_sub4.xyz gap={soap l_max=8 n_max=8 atom_sigma=0.5 "
            "zeta=4 cutoff=4.0 cutoff_transition_width=1.0 central_weight=1.0 "
            "n_sparse=500 delta=3.0 f0=0.0 covariance_type=dot_product "
            "$SPARSE_METHOD} "
            "default_sigma={0.01 1.0 0.5 0.0} "
            "energy_parameter_name=dft_energy force_parameter_name=dft_force "
            "virial_parameter_name=dft_virial config_type_parameter_name=config_type "
            "sparse_jitter=1.0e-3 e0_offset=2.0 gp_file=gp.xml rnd_seed=1 $COND_NUM_NORM"
        )
        with open('si_gap_fit_test.json') as f:
            self.ref_data = json.load(f)

    def check_gap_fit(self, command_line, ref_data, new_test=False, prefix=''):
        if new_test:
            command_line = command_line.replace('$SPARSE_METHOD',
                                                'sparse_method=cur_points print_sparse_index=sparse_file')
            if os.path.exists('sparse_file'):
                os.unlink('sparse_file') # ensure we don't append to an old file
        else:
            with open('sparse_file', 'w') as fh:
                for sp in ref_data['sparse_points']:
                    fh.write(f'{sp}\n')
            command_line = command_line.replace('$SPARSE_METHOD',
                                                'sparse_method=index_file sparse_file=sparse_file')

        build_dir = os.environ.get('BUILDDIR')
        program = os.path.join(build_dir, 'gap_fit')
        full_command = f'{prefix} {program} {command_line}'
        print(full_command)
        proc = subprocess.run(full_command, shell=True, env=os.environ)
        assert proc.returncode == 0, proc

        tree = ET.parse('gp.xml')
        root = tree.getroot()
        gp_label = root.tag

        # check GP coefficients match expected values
        idx = np.array([int(tag.attrib['i']) for tag in root[1][1][0].findall('sparseX')])
        idx -= 1 # convert from one- to zero-based indexing
        alpha = np.array([float(tag.attrib['alpha']) for tag in root[1][1][0].findall('sparseX')])
        alpha = alpha[idx] # reorder correctly
        if new_test:
            sparse_points = np.loadtxt('sparse_file').astype(int)
            gap_dict = {'sparse_points': sparse_points.tolist(),
                        'alpha': alpha.tolist()}
            with open('dump.json', 'w') as f:
                json.dump(gap_dict, f, indent=4)
        else:
            print('max abs error in alpha =', np.abs((alpha - ref_data['alpha'])).max())
            assert np.abs((alpha - ref_data['alpha'])).max() < self.alpha_tol

    def test_gap_fit_silicon(self):
        command_line = self.cl_template.safe_substitute(COND_NUM_NORM="condition_number_norm=I")
        self.check_gap_fit(command_line, self.ref_data)

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_gap_fit_silicon_scalapack(self):
        command_line = self.cl_template.safe_substitute(SPARSE_METHOD='sparse_method=FILE sparse_file=si_gap_fit_sparseX.inp')
        self.check_gap_fit(command_line, self.ref_data, prefix='mpirun -np 2')

if __name__ == '__main__':
    unittest.main()

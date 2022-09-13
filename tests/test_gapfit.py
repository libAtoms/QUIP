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
import hashlib
from pathlib import Path

import numpy as np
import quippytest
import quippy

def file2hash(filename, chunksize=4096):
    hasher = hashlib.sha256()
    with open(filename, 'rb') as f:
        chunk = f.read(chunksize)
        while chunk != b'':
            hasher.update(chunk)
            chunk = f.read(chunksize)
    return hasher.hexdigest()


@unittest.skipIf(os.environ['HAVE_GAP'] != '1', 'GAP support not enabled')
class TestGAP_fit(quippytest.QuippyTestCase):
    alpha_tol = 1e-5
    log_name = 'gap_fit.log'
    xml_name = 'gp.xml'
    config_name = 'gap_fit.config'
    here = Path('.')
    if 'BUILDDIR' in os.environ:
        prog_path = Path(os.environ.get('BUILDDIR')) / 'gap_fit'
    else:
        prog_path = Path(quippy.__path__[0]) / 'gap_fit'
        
    if not os.path.isfile(prog_path):
        raise unittest.SkipTest(f"gap_fit exectuable does not exist at {prog_path}")
        
    with open('si_gap_fit_test.json') as f:
        ref_data = json.load(f)
    si_sparsex_hash = 'bf3d99356e16bc666cee1f1abc6a2cfc63e98a8f69658bcc5ab84e01d9e3ab2d'

    def setUp(self):
        self.cl_template = string.Template(
            "at_file=train_sub4.xyz gap={soap l_max=8 n_max=8 atom_sigma=0.5 "
            "zeta=4 cutoff=4.0 cutoff_transition_width=1.0 central_weight=1.0 "
            "n_sparse=500 delta=3.0 f0=0.0 covariance_type=dot_product "
            "$SPARSE_METHOD} "
            "default_sigma={0.01 1.0 0.5 0.0} "
            "energy_parameter_name=dft_energy force_parameter_name=dft_force "
            "virial_parameter_name=dft_virial config_type_parameter_name=config_type "
            "sparse_jitter=1.0e-3 e0_offset=2.0 gp_file=gp.xml rnd_seed=1 "
        )
        self.env = os.environ.copy()
        self.env['OMP_NUM_THREADS'] = '1'

    def tearDown(self):
        for fname in [self.log_name, self.config_name]:
            try:
                os.remove(self.here / fname)
            except FileNotFoundError:
                pass
        for path in self.here.glob(self.xml_name + '*'):
            os.remove(path)

    def run_gap_fit(self, command_line, new_test=False, prefix=''):        
        if new_test:
            if os.path.exists('sparse_file'):
                os.unlink('sparse_file') # ensure we don't append to an old file
        else:
            with open('sparse_file', 'w') as fh:
                for sp in self.ref_data['sparse_points']:
                    fh.write(f'{sp}\n')

        full_command = f'{prefix} {self.prog_path} {command_line}'
        with open(self.log_name, 'w') as f:
            self.proc = subprocess.run(full_command, shell=True, env=self.env, stdout=f, stderr=f)

    def check_gap_fit(self, new_test=False):
        assert self.proc.returncode == 0, self.proc

        tree = ET.parse(self.xml_name)
        root = tree.getroot()

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
            print('max abs error in alpha =', np.abs((alpha - self.ref_data['alpha'])).max())
            assert np.abs((alpha - self.ref_data['alpha'])).max() < self.alpha_tol

    def check_latest_sparsex_file_hash(self):
        files = self.here.glob(self.xml_name + '.*')
        flast = max(files, key=os.path.getctime)
        hash = file2hash(flast)
        self.assertEqual(hash, self.si_sparsex_hash)

    def test_gap_fit_silicon_sparsify_only(self):
        with open(self.config_name, 'w') as f:
            config = self.cl_template.safe_substitute(SPARSE_METHOD='sparse_method=index_file sparse_file=sparse_file')
            config += ' sparsify_only_no_fit=T'
            print(config, file=f)
        command_line = f'config_file={self.config_name}'
        self.run_gap_fit(command_line)
        with open(self.log_name) as f:
            self.assertEqual(f.read().count('Number of partial derivatives of descriptors: 0'), 1)

    # new test: 'sparse_method=cur_points print_sparse_index=sparse_file'
    def test_gap_fit_silicon(self):
        self.env['OMP_NUM_THREADS'] = '2'
        command_line = self.cl_template.safe_substitute(SPARSE_METHOD='sparse_method=index_file sparse_file=sparse_file')
        command_line += ' condition_number_norm=I'
        self.run_gap_fit(command_line)
        self.check_gap_fit()

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_gap_fit_silicon_scalapack(self):
        command_line = self.cl_template.safe_substitute(SPARSE_METHOD='sparse_method=FILE sparse_file=si_gap_fit_sparseX.inp')
        command_line += ' mpi_blocksize=101'  # not a divisor of n_sparse
        self.run_gap_fit(command_line, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit()
        self.check_latest_sparsex_file_hash()

if __name__ == '__main__':
    unittest.main()

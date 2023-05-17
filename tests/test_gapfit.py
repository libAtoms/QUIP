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

# To update a reference, add `self.make_ref_file(config, ref_file)` after
# `self.run_gap_fit(config)`. To update a sparseX.inp file, prevent their
# deletion in tearDown() and convert them via bin/gap_prepare_sparsex_input.py.
# New xyz need reordered non-MPI frames per n_sparse due to unstable sorting.

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
    index_name = 'sparse_index'
    config_name = 'gap_fit.config'
    here = Path('.')
    if 'BUILDDIR' in os.environ:
        prog_path = Path(os.environ.get('BUILDDIR')) / 'gap_fit'
    else:
        prog_path = Path(quippy.__path__[0]) / 'gap_fit'

    if not os.path.isfile(prog_path):
        raise unittest.SkipTest(f"gap_fit exectuable does not exist at {prog_path}")

    ref_files = {
        ('Si', 'soap', 'cur_points'): 'Si.soap.cur_points.json',
        ('Si', 'distance_2b', 'uniform'): 'Si.distance_2b.uniform.json',
    }

    config_default_mapping = {'XML_NAME': xml_name, 'EXTRA': ''}
    config_template = string.Template(
        "at_file=$XYZ_FILE config_type_parameter_name=config_type default_sigma={0.01 1.0 0.5 0.0} "
        " do_copy_atoms_file=F e0_offset=2.0 energy_parameter_name=dft_energy force_parameter_name=dft_force"
        " gap={$GAP} gp_file=$XML_NAME rnd_seed=1 sparse_jitter=1.0e-3 virial_parameter_name=dft_virial $EXTRA")

    sparse_method_index = f'sparse_method=index_file sparse_file={index_name}'
    gap_default_mapping = {'INDEX_NAME': index_name, 'SPARSE_METHOD': sparse_method_index}
    gap_distance_2b_template = string.Template(
        "{distance_2b covariance_type=ard_se cutoff=6.0 delta=1.0 n_sparse=20 print_sparse_index=$INDEX_NAME"
        " theta_uniform=0.1 $SPARSE_METHOD}")
    gap_soap_template = string.Template(
        "{soap atom_sigma=0.5 central_weight=1.0 covariance_type=dot_product cutoff=4.0 cutoff_transition_width=1.0"
        " delta=3.0 l_max=8 n_max=8 n_sparse=100 print_sparse_index=$INDEX_NAME zeta=4 $SPARSE_METHOD}")

    def setUp(self):
        self.env = os.environ.copy()
        self.env['OMP_NUM_THREADS'] = '1'

    def tearDown(self):
        for fname in [self.log_name, self.config_name, self.index_name]:
            try:
                os.remove(self.here / fname)
            except FileNotFoundError:
                pass
        for path in self.here.glob(self.xml_name + '*'):
            os.remove(path)
        for path in self.here.glob('*.xyz.idx'):
            os.remove(path)

    def check_gap_fit(self, ref_file=None):
        """check GP coefficients match expected values"""
        assert self.proc.returncode == 0, self.proc

        alpha = self.get_alpha_from_xml(self.xml_name)
        ref_data = self.read_ref_data(ref_file)
        self.assertLess(np.abs((alpha - ref_data['alpha'])).max(), self.alpha_tol)
        self.check_latest_sparsex_file_hash(ref_data['sparse_hash'])

    def check_latest_sparsex_file_hash(self, ref):
        files = self.here.glob(self.xml_name + '.*')
        flast = max(files, key=os.path.getctime)
        hash = file2hash(flast)
        self.assertEqual(hash, ref)

    def get_alpha_from_xml(self, xml_file):
        root = ET.parse(xml_file).getroot()

        idx = np.array([int(tag.attrib['i']) - 1 for tag in root[1][1][0].findall('sparseX')])
        alpha = np.array([float(tag.attrib['alpha']) for tag in root[1][1][0].findall('sparseX')])
        alpha = alpha[idx]  # reorder correctly
        return alpha

    def get_config(self, xyz_file, gap, *, extra=''):
        config = self.config_template.substitute(self.config_default_mapping, XYZ_FILE=xyz_file, GAP=gap, EXTRA=extra)
        return config

    def get_gap(self, template, sparse_method):
        gap = template.substitute(self.gap_default_mapping, SPARSE_METHOD=sparse_method)
        return gap

    def make_index_file_from_ref(self, ref_file, index_file):
        ref_data = self.read_ref_data(ref_file)
        with open(index_file, 'w') as f:
            for sp in ref_data['indices']:
                f.write(f'{sp}\n')

    def make_ref_file(self, config, ref_file='ref.json'):
        alpha = self.get_alpha_from_xml(self.xml_name)
        indices = np.loadtxt(self.index_name).astype(int)

        root = ET.parse(self.xml_name).getroot()
        gp_coord = root.find('GAP_params/gpSparse/gpCoordinates')
        sparse_hash = file2hash(gp_coord.attrib['sparseX_filename'])

        gap_dict = {
            'alpha': alpha.tolist(),
            'config': config,
            'indices': indices.tolist(),
            'sparse_hash': sparse_hash,
        }
        with open(ref_file, 'w') as f:
            json.dump(gap_dict, f, indent=4)

    def read_ref_data(self, ref_file):
        with open(ref_file) as f:
            return json.load(f)

    def run_gap_fit(self, config, prefix=''):
        command = f'{prefix} {self.prog_path} {config}'
        with open(self.log_name, 'w') as f:
            self.proc = subprocess.run(command, shell=True, env=self.env, stdout=f, stderr=f)

    def test_si_soap_cur_points(self):
        self.env['OMP_NUM_THREADS'] = '2'
        ref_file = self.ref_files[('Si', 'soap', 'cur_points')]
        gap = self.get_gap(self.gap_soap_template, 'sparse_method=cur_points')
        config = self.get_config('Si.np1.sp100.xyz', gap, extra='condition_number_norm=I')
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_si_distance_2b_uniform(self):
        self.env['OMP_NUM_THREADS'] = '2'
        ref_file = self.ref_files[('Si', 'distance_2b', 'uniform')]
        gap = self.get_gap(self.gap_distance_2b_template, 'sparse_method=uniform')
        config = self.get_config('Si.np1.sp20.xyz', gap)
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_si_distance_2b_index(self):
        self.env['OMP_NUM_THREADS'] = '2'
        ref_file = self.ref_files[('Si', 'distance_2b', 'uniform')]
        gap = self.get_gap(self.gap_distance_2b_template, f'sparse_method=index_file sparse_file={self.index_name}')
        config = self.get_config('Si.np1.sp20.xyz', gap)
        self.make_index_file_from_ref(ref_file, self.index_name)
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_si_config_file_sparsify_only(self):
        gap = self.get_gap(self.gap_soap_template, sparse_method='sparse_method=cur_points')
        config = self.get_config('Si.np1.sp100.xyz', gap, extra='sparsify_only_no_fit=T')
        with open(self.config_name, 'w') as f:
            print(config, file=f)
        config = f'config_file={self.config_name}'
        self.run_gap_fit(config)
        with open(self.log_name) as f:
            self.assertEqual(f.read().count('Number of partial derivatives of descriptors: 0'), 1)

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_si_scalapack_soap_file(self):
        gap = self.get_gap(self.gap_soap_template, 'sparse_method=file sparse_file=Si.cur_points.sparseX.inp')
        config = self.get_config('Si.np2.xyz', gap, extra='mpi_blocksize_rows=101')
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit(self.ref_files[('Si', 'soap', 'cur_points')])

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_si_scalapack_soap_cur_points(self):
        gap = self.get_gap(self.gap_soap_template, 'sparse_method=cur_points')
        config = self.get_config('Si.np2.xyz', gap)
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit(self.ref_files[('Si', 'soap', 'cur_points')])

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_si_scalapack_distance_2b_uniform(self):
        gap = self.get_gap(self.gap_distance_2b_template, 'sparse_method=uniform')
        config = self.get_config('Si.np2.xyz', gap)
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit(self.ref_files[('Si', 'distance_2b', 'uniform')])

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_si_scalapack_big_blocksize(self):
        gap = self.get_gap(self.gap_soap_template, 'sparse_method=cur_points')
        config = self.get_config('Si.np2.xyz', gap, extra='mpi_blocksize_cols=37838')  # too large for 32bit integer
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("too large for 32bit work array in ScaLAPACK" in f.read())


if __name__ == '__main__':
    unittest.main()

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
# `self.run_gap_fit(config)`. New xyz need reordered non-MPI frames to compare
# with MPI ones.

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


class NumPyJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        else:
            return super().default(self, obj)


@unittest.skipIf(os.environ['HAVE_GAP'] != '1', 'GAP support not enabled')
class TestGAP_fit(quippytest.QuippyTestCase):
    alpha_tol = 1e-5
    sparsex_tol = 1e-8
    log_name = 'gap_fit.log'
    xml_name = 'gp.xml'
    index_name = 'sparse_index'
    index_name_inp = 'sparse_index.inp'
    index_name_out = 'sparse_index.out'
    config_name = 'gap_fit.config'
    sparsex_name_inp = 'sparsex.inp'
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
        ('Si', 'two_descriptors'): 'Si.two_descriptors.json',
        ('SiC', 'distance_2b', 'uniform'): 'SiC.distance_2b.uniform.json',
    }

    config_default_mapping = {'XML_NAME': xml_name, 'EXTRA': ''}
    config_template = string.Template(
        "at_file=$XYZ_FILE config_type_parameter_name=config_type default_sigma={0.01 1.0 0.5 0.0} "
        " do_copy_atoms_file=F e0_offset=2.0 energy_parameter_name=dft_energy force_parameter_name=dft_force"
        " gap={$GAP} gp_file=$XML_NAME rnd_seed=1 sparse_jitter=1.0e-3 virial_parameter_name=dft_virial $EXTRA")

    sparse_method_index = f'sparse_method=index_file sparse_file={index_name_inp}'
    gap_default_mapping = {'INDEX_NAME_OUT': index_name_out, 'SPARSE_METHOD': sparse_method_index}
    gap_distance_2b_template = string.Template(
        "{distance_2b covariance_type=ard_se cutoff=6.0 delta=1.0 n_sparse=20 print_sparse_index=$INDEX_NAME_OUT"
        " theta_uniform=0.1 $SPARSE_METHOD}")
    gap_soap_template = string.Template(
        "{soap atom_sigma=0.5 central_weight=1.0 covariance_type=dot_product cutoff=4.0 cutoff_transition_width=1.0"
        " delta=3.0 l_max=8 n_max=8 n_sparse=100 print_sparse_index=$INDEX_NAME_OUT zeta=4 $SPARSE_METHOD}")

    def setUp(self):
        self.env = os.environ.copy()
        self.env['OMP_NUM_THREADS'] = '1'

    def tearDown(self):
        for fname in [self.config_name, self.index_name_inp, self.index_name_out,
                      self.log_name, self.sparsex_name_inp]:
            try:
                os.remove(self.here / fname)
            except FileNotFoundError:
                pass
        for path in self.here.glob(self.xml_name + '*'):
            os.remove(path)
        for path in self.here.glob('*.xyz.idx'):
            os.remove(path)

    def check_gap_committee(self, err_tol=0.2):
        np.random.seed(1)
        ncomm = 200 # num_committee members to generate
        xml_committee = quippy.gap_tools.get_xml_committee(self.xml_name, ncomm)

        comm_weights = np.array([comm.weights for comm in xml_committee])
        true_weights = xml_committee[0].mean_weights

        weight_stds = np.std(comm_weights, axis=0)
        weight_mean = np.average(comm_weights, axis=0)

        mean_errs = np.abs(weight_mean - true_weights) / weight_stds
        print("max mean err = ", np.max(mean_errs))
        assert np.all(mean_errs < err_tol)
        np.random.seed(None)

    def check_gap_fit(self, ref_file=None):
        """check GP coefficients match expected values"""
        assert self.proc.returncode == 0, self.proc

        ref_data = self.read_ref_data(ref_file)
        self.check_index_file(ref_data)

        data = self.get_data_from_xml(self.xml_name)
        for coord, ref in zip(data['coords'], ref_data['coords']):
            self.assertLess(np.abs(coord['sparsex'] - ref['sparsex']).max(), self.sparsex_tol)
            self.assertLess(np.abs(coord['alpha'] - ref['alpha']).max(), self.alpha_tol)

    def check_index_file(self, ref):
        file = self.here / self.index_name_out
        if not file.exists():
            return
        index = self.get_index_from_file(file)
        self.assertEqual(index, ref['index'])

    def check_latest_sparsex_file_hash(self, ref):
        files = self.here.glob(self.xml_name + '.*')
        flast = max(files, key=os.path.getctime)
        hash = file2hash(flast)
        self.assertEqual(hash, ref['sparse_hash'])

    def get_data_from_xml(self, xml_file):
        root = ET.parse(xml_file).getroot()
        coords = root.findall('GAP_params/gpSparse/gpCoordinates')
        data = {'coords': []}
        for coord in coords:
            alpha = np.array([float(tag.attrib['alpha']) for tag in coord.findall('sparseX')])
            cutoff = np.array([float(tag.attrib['sparseCutoff']) for tag in coord.findall('sparseX')])
            dimensions = int(coord.attrib['dimensions'])
            sparsex_file = Path(xml_file).parent / coord.attrib['sparseX_filename']
            sparsex = np.loadtxt(sparsex_file)
            data['coords'].append({
                'dimensions': dimensions,
                'alpha': alpha,
                'cutoff': cutoff,
                'sparsex': sparsex,
            })
        return data

    def get_config(self, xyz_file, gap, *, extra=''):
        config = self.config_template.substitute(self.config_default_mapping, XYZ_FILE=xyz_file, GAP=gap, EXTRA=extra)
        return config

    def get_gap(self, template, sparse_method):
        gap = template.substitute(self.gap_default_mapping, SPARSE_METHOD=sparse_method)
        return gap

    def get_index_from_file(self, index_file):
        with open(index_file) as f:
            return [[int(i) for i in line.split()] for line in f]

    def make_first_index_file_from_ref(self, ref_file, index_file):
        ref_data = self.read_ref_data(ref_file)
        with open(index_file, 'w') as f:
            print(*ref_data['index'][0], sep='\n', file=f)

    def make_first_sparsex_input_file_from_ref(self, ref, sparsex_file):
        coord = ref['coords'][0]
        with open(sparsex_file, 'w') as f:
            size = coord['dimensions']
            vectors = (coord['sparsex'][i:i+size] for i in range(0, len(coord['sparsex']), size))
            for cutoff, vector in zip(coord['cutoff'], vectors):
                print(cutoff, file=f)
                print(*vector, sep='\n', file=f)

    def make_ref_file(self, config, ref_file='ref.json'):
        data = self.get_data_from_xml(self.xml_name)
        data['config'] = config
        data['index'] = self.get_index_from_file(self.index_name_out)

        with open(ref_file, 'w') as f:
            json.dump(data, f, indent=4, cls=NumPyJSONEncoder)

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
        config = self.get_config('Si.np1.xyz', gap, extra='condition_number_norm=I')
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_si_distance_2b_uniform(self):
        self.env['OMP_NUM_THREADS'] = '2'
        ref_file = self.ref_files[('Si', 'distance_2b', 'uniform')]
        gap = self.get_gap(self.gap_distance_2b_template, 'sparse_method=uniform')
        config = self.get_config('Si.np1.xyz', gap)
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_si_two_descriptors(self):
        self.env['OMP_NUM_THREADS'] = '2'
        ref_file = self.ref_files[('Si', 'two_descriptors')]
        gap1 = self.get_gap(self.gap_distance_2b_template, 'sparse_method=uniform')
        gap2 = self.get_gap(self.gap_soap_template, 'sparse_method=cur_points')
        gap = ":".join([gap1, gap2])
        config = self.get_config('Si.np1.xyz', gap)
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_sic_distance_2b_uniform(self):
        self.env['OMP_NUM_THREADS'] = '2'
        ref_file = self.ref_files[('SiC', 'distance_2b', 'uniform')]
        gap = self.get_gap(self.gap_distance_2b_template, 'sparse_method=uniform')
        config = self.get_config('SiC.np1.xyz', gap)
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_si_distance_2b_index(self):
        self.env['OMP_NUM_THREADS'] = '2'
        ref_file = self.ref_files[('Si', 'distance_2b', 'uniform')]
        gap = self.get_gap(self.gap_distance_2b_template, f'sparse_method=index_file sparse_file={self.index_name_inp}')
        config = self.get_config('Si.np1.xyz', gap)
        self.make_first_index_file_from_ref(ref_file, self.index_name_inp)
        self.run_gap_fit(config)
        self.check_gap_fit(ref_file)

    def test_si_config_file_sparsify_only(self):
        gap = self.get_gap(self.gap_soap_template, sparse_method='sparse_method=cur_points')
        config = self.get_config('Si.np1.xyz', gap, extra='sparsify_only_no_fit=T')
        with open(self.config_name, 'w') as f:
            print(config, file=f)
        config = f'config_file={self.config_name}'
        self.run_gap_fit(config)
        with open(self.log_name) as f:
            self.assertEqual(f.read().count('Number of partial derivatives of descriptors: 0'), 1)

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_si_scalapack_soap_file(self):
        ref_file = self.ref_files[('Si', 'soap', 'cur_points')]
        ref_data = self.read_ref_data(ref_file)
        self.make_first_sparsex_input_file_from_ref(ref_data, self.sparsex_name_inp)
        gap = self.get_gap(self.gap_soap_template, f'sparse_method=file sparse_file={self.sparsex_name_inp}')
        config = self.get_config('Si.np2.xyz', gap, extra='mpi_blocksize_rows=101')
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit(ref_file)

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
        config = self.get_config('Si.np2.xyz', gap, extra="export_covariance=T")
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit(self.ref_files[('Si', 'distance_2b', 'uniform')])
        self.check_gap_committee()

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_si_scalapack_two_descriptors(self):
        ref_file = self.ref_files[('Si', 'two_descriptors')]
        gap1 = self.get_gap(self.gap_distance_2b_template, 'sparse_method=uniform')
        gap2 = self.get_gap(self.gap_soap_template, 'sparse_method=cur_points')
        gap = ":".join([gap1, gap2])
        config = self.get_config('Si.np2.xyz', gap)
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit(ref_file)

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_sic_scalapack_distance_2b_uniform(self):
        ref_file = self.ref_files[('SiC', 'distance_2b', 'uniform')]
        gap = self.get_gap(self.gap_distance_2b_template, 'sparse_method=uniform')
        config = self.get_config('SiC.np2.xyz', gap)
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("Using ScaLAPACK to solve QR" in f.read())
        self.check_gap_fit(ref_file)

    @unittest.skipIf(os.environ.get('HAVE_SCALAPACK') != '1', 'ScaLAPACK support not enabled')
    def test_si_scalapack_big_blocksize(self):
        gap = self.get_gap(self.gap_soap_template, 'sparse_method=cur_points')
        config = self.get_config('Si.np2.xyz', gap, extra='mpi_blocksize_cols=37838')  # too large for 32bit integer
        self.run_gap_fit(config, prefix='mpirun -np 2')
        with open(self.log_name) as f:
            self.assertTrue("too large for 32bit work array in ScaLAPACK" in f.read())


if __name__ == '__main__':
    unittest.main()

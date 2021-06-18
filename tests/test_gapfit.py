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

import numpy as np
import quippytest

from quippy.potential import Potential
from ase.constraints import voigt_6_to_full_3x3_stress
from ase.io import read, write

@unittest.skipIf(os.environ['HAVE_GAP'] != '1', 'GAP support not enabled')
class TestGAP_fit(quippytest.QuippyTestCase):
    
    def check_gap_fit(self, command_line, ref_data, new_test=False):
        if new_test:
            command_line = command_line.replace('SPARSE_METHOD', 
                                                'sparse_method=cur_points print_sparse_index=sparse_file')
            if os.path.exists('sparse_file'):
                os.unlink('sparse_file') # ensure we don't append to an old file
        else:
            with open('sparse_file', 'w') as fh:
                for sp in ref_data['sparse_points']:
                    fh.write(f'{sp}\n')
            command_line = command_line.replace('SPARSE_METHOD',
                                                'sparse_method=index_file sparse_file=sparse_file')
            
        print(command_line)
        stat = os.system('gap_fit '+command_line)
        assert stat == 0
        
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
            assert np.abs((alpha - ref_data['alpha'])).max() < 1e-2
                
    def test_gap_fit_silicon(self):
        train_filename = 'train_sub4.xyz'
        command_line = (f"at_file={train_filename} gap={{soap l_max=8 n_max=8 atom_sigma=0.5 "
                        "zeta=4 cutoff=4.0 cutoff_transition_width=1.0 central_weight=1.0 "
                        "n_sparse=500 delta=3.0 f0=0.0 covariance_type=dot_product "
                        "SPARSE_METHOD } "
                        "default_sigma={0.001 0.1 0.05 0.0} "
                        "energy_parameter_name=dft_energy force_parameter_name=dft_force "
                        "virial_parameter_name=dft_virial config_type_parameter_name=config_type "
                        "sparse_jitter=1.0e-8 e0_offset=2.0 gp_file=gp.xml rnd_seed=1")
        with open('si_gap_fit_test.json') as f:
            ref_data = json.load(f)
        self.check_gap_fit(command_line, ref_data, new_test=False)

if __name__ == '__main__':
    unittest.main()

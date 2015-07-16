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

# Tests of LOTF used to embed one SW potential within another

from quippy import *
import unittest, quippy
from quippytest import *

if hasattr(quippy, 'Potential') and hasattr(quippy, 'Potential') and quippy.have_lotf:
    
    class TestLOTF(QuippyTestCase):
        def setUp(self):
            xml ="""<quip_params>
<SW_params n_types="2" label="PRB_31_plus_H">
<comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984), extended for other elements </comment>
<per_type_data type="1" atomic_num="1" />
<per_type_data type="2" atomic_num="14" />
<per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
      p="0" q="0" a="5.0" sigma="1.0" eps="0.0" />
<per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827"
      p="4" q="0" a="1.25" sigma="2.537884" eps="2.1672" />
<per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
      p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

<!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
<per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
      lambda="21.0" gamma="1.20" eps="2.1675" />
<per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.1675" />
<per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.1675" />

<per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
      lambda="21.0" gamma="1.20" eps="2.1675" />
<per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.1675" />
<per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.1675" />
</SW_params>

<SW_params n_types="2" label="eps_2.3">

<per_type_data type="1" atomic_num="1" />
<per_type_data type="2" atomic_num="14" />
<per_pair_data atnum_i="1" atnum_j="1" AA="0.0" BB="0.0"
      p="0" q="0" a="5.0" sigma="1.0" eps="0.0" />
<per_pair_data atnum_i="1" atnum_j="14" AA="8.581214" BB="0.0327827"
      p="4" q="0" a="1.25" sigma="2.537884" eps="2.3" />
<per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
      p="4" q="0" a="1.80" sigma="2.0951" eps="2.3" />

<!-- triplet terms: atnum_c is the center atom, neighbours j and k -->
<per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="1"
      lambda="21.0" gamma="1.20" eps="2.3" />
<per_triplet_data atnum_c="1"  atnum_j="1"  atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.3" />
<per_triplet_data atnum_c="1"  atnum_j="14" atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.3" />

<per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="1"
      lambda="21.0" gamma="1.20" eps="2.3" />
<per_triplet_data atnum_c="14" atnum_j="1"  atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.3" />
<per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.3" />
</SW_params>
</quip_params>"""

            system_reseed_rng(1984068303)
            self.pot1 = Potential('IP SW label="PRB_31_plus_H"', param_str=xml)
            self.pot2 = Potential('IP SW label="eps_2.3"', param_str=xml)
         
            self.dia = diamond(5.44, 14)
            self.at = supercell(self.dia, 4, 4, 4)
            randomise(self.at.pos, 0.1)
            self.at.set_cutoff(self.pot1.cutoff()+2.0)
            self.at.calc_connect()
            
            self.at.add_property('hybrid', 0)
            self.embedlist = Table(4,0,0,0)
            self.embedlist.append((1,0,0,0))
            self.at.bfs_grow_list(self.embedlist, 1, nneighb_only=True)
            self.at.hybrid[self.embedlist.int[1,:]] = 1

            self.lotf = Potential('ForceMixing method=lotf_adj_pot_svd fit_hops=3 buffer_hops=3 '+
                                      'randomise_buffer=T qm_args_str={carve_cluster=T single_cluster=T cluster_calc_connect=T}', 
                                      self.pot1, self.pot2, bulk_scale=self.dia)
            self.forcemix = Potential('ForceMixing method=force_mixing buffer_hops=3 '+
                                          'qm_args_str={carve_cluster=T single_cluster=T cluster_calc_connect=T}', 
                                          self.pot1, self.pot2, bulk_scale=self.dia)

        def test_no_predictor_corrector(self):
            ds = DynamicalSystem(self.at)
            f = fzeros((3,self.at.n))
            f_hyb = fzeros((3,self.at.n))
            verbosity_push(PRINT_SILENT)
            self.lotf.calc(ds.atoms, force=f)
            self.forcemix.calc(ds.atoms, force=f_hyb)
            ds.advance_verlet(1.0, f_hyb)
            verbosity_pop()
            self.assertArrayAlmostEqual([rms_diff(f_hyb, f), abs(f_hyb -f).max()],
                                        [0.00032825, 0.00179225])

        def do_predictor_corrector(self, n_extrap):
            f = fzeros((3,self.at.n))
            f_hyb = fzeros((3,self.at.n))
            extrap_force_err = []
            interp_force_err = []
            self.lotf.calc(self.ds.atoms, force=f) # bootstrap

            # Extrapolation
            self.ds_saved.save_state(self.ds)
            for j in frange(n_extrap):
                if j == 1:
                    self.lotf.calc(self.ds.atoms, force=f, lotf_do_qm=False, lotf_do_init=True, lotf_do_map=True)
                else:
                    self.lotf.calc(self.ds.atoms, force=f, lotf_do_qm=False, lotf_do_init=False)
                self.forcemix.calc(self.ds.atoms, force=f_hyb)
                self.ds.advance_verlet(1.0, f)
                self.ds.print_status('E', instantaneous=True)
                extrap_force_err.append([rms_diff(f_hyb, f), abs(f_hyb -f).max()])
            
            # Force computation
            self.lotf.calc(self.ds.atoms, force=f, lotf_do_qm=True, lotf_do_init=False, lotf_do_fit=True)

            self.ds.restore_state(self.ds_saved)

            # Interpolation
            for j in frange(n_extrap):
                self.lotf.calc(self.ds.atoms, force=f, lotf_do_qm=False, lotf_do_init=False, lotf_do_interp=True, lotf_interp=float(j-1)/float(n_extrap))
                self.forcemix.calc(self.ds.atoms, force=f_hyb)
                self.ds.advance_verlet(1.0, f)
                self.ds.print_status('I', instantaneous=True)
                interp_force_err.append([rms_diff(f_hyb, f), abs(f_hyb -f).max()])

            return extrap_force_err, interp_force_err

        def test_extrap_10_steps(self):
            self.ds = DynamicalSystem(self.at)
            self.ds_saved = DynamicalSystem(self.at)

            verbosity_push(PRINT_SILENT)
            extrap_force_err, interp_force_err = self.do_predictor_corrector(10)
            verbosity_pop()

            self.assertArrayAlmostEqual(extrap_force_err, [[ 0.00032825,  0.00179226],
                                                           [ 0.00032781,  0.00229827],
                                                           [ 0.00034314,  0.00421474],
                                                           [ 0.00039979,  0.00706289],
                                                           [ 0.00051401,  0.0108113 ],
                                                           [ 0.00068497,  0.01541874],
                                                           [ 0.00090493,  0.02083476],
                                                           [ 0.00116668,  0.02700031],
                                                           [ 0.00146453,  0.03384861],
                                                           [ 0.00179375,  0.04130589]])
            self.assertArrayAlmostEqual(interp_force_err, [[ 0.00032825,  0.00179226],
                                                           [ 0.0003675 ,  0.00490657],
                                                           [ 0.00043828,  0.00749941],
                                                           [ 0.0004971 ,  0.00924506],
                                                           [ 0.00052894,  0.01017118],
                                                           [ 0.00052948,  0.01031334],
                                                           [ 0.00049851,  0.00971477],
                                                           [ 0.00043801,  0.00842611],
                                                           [ 0.00035217,  0.00650499],
                                                           [ 0.00025062,  0.00401573]])
            
        def test_extrap_2_cycles_5_steps(self):
            self.ds = DynamicalSystem(self.at)
            self.ds_saved = DynamicalSystem(self.at)

            verbosity_push(PRINT_SILENT)
            extrap_force_err_1, interp_force_err_1 = self.do_predictor_corrector(5)
            self.ds.atoms.calc_connect()
            extrap_force_err_2, interp_force_err_2 = self.do_predictor_corrector(5)
            verbosity_pop()

            self.assertArrayAlmostEqual(extrap_force_err_1, [[ 0.00032825,  0.00179226],
                                                             [ 0.00032781,  0.00229827],
                                                             [ 0.00034314,  0.00421474],
                                                             [ 0.00039979,  0.00706289],
                                                             [ 0.00051401,  0.0108113 ]])
            self.assertArrayAlmostEqual(interp_force_err_1, [[ 0.00032825,  0.00179226],
                                                             [ 0.00033528,  0.00310028],
                                                             [ 0.00034133,  0.0038855 ],
                                                             [ 0.00033204,  0.00382098],
                                                             [ 0.00030716,  0.00293322]])
            self.assertArrayAlmostEqual(extrap_force_err_2, [[ 0.00028088,  0.00151502],
                                                             [ 0.00035584,  0.00652683],
                                                             [ 0.00056809,  0.01265381],
                                                             [ 0.00084781,  0.01945528],
                                                             [ 0.00116842,  0.02685806]])
            self.assertArrayAlmostEqual(interp_force_err_2, [[ 0.00028088,  0.00151502],
                                                             [ 0.00026973,  0.00239851],
                                                             [ 0.00025789,  0.00287206],
                                                             [ 0.00023529,  0.0027123 ],
                                                             [ 0.00020142,  0.0019848 ]])



if __name__ == '__main__':
   unittest.main()



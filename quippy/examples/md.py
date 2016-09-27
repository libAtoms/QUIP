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
Simple example of molecular dynamics of 64 atoms of bulk silicon at 1000 K
James Kermode 2009
"""

from quippy import *
from atomeye import view

# Set up atomic configuation
s = supercell(diamond(5.44, 14), 3, 3, 3)

# Initialise potential from XML string
pot = Potential('IP SW', param_str="""<SW_params n_types="1">
<comment> Stillinger and Weber, Phys. Rev. B  31 p 5262 (1984)</comment>
<per_type_data type="1" atomic_num="14" />

<per_pair_data atnum_i="14" atnum_j="14" AA="7.049556277" BB="0.6022245584"
      p="4" q="0" a="1.80" sigma="2.0951" eps="2.1675" />

<per_triplet_data atnum_c="14" atnum_j="14" atnum_k="14"
      lambda="21.0" gamma="1.20" eps="2.1675" />
</SW_params>
""")

s.set_cutoff(pot.cutoff()+2.0)
s.calc_connect()

# Set up dynamical system at 1000K
ds = DynamicalSystem(s)
ds.rescale_velo(1000.0)
ds.zero_momentum()

outf = CInOutput('si-1000.xyz', OUTPUT)
traj = AtomsList(ds.run(pot, dt=1.0, n_steps=100, save_interval=10, trajectory=outf))
outf.close()
view(traj)
raw_input('Press ENTER to terminate')


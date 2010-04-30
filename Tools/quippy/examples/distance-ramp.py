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

from quippy import *

verbosity_set_minimum(VERBOSE)

d = diamond(5.44, 14)

at = supercell(d, 3, 3, 3)
at.calc_connect()

at.add_property('hybrid_mark', HYBRID_NO_MARK)

embed = Table()
embed.append((1,0,0,0))
embed.append(at.bfs_step(embed))

at.hybrid_mark[embed.int[1,:]] = HYBRID_ACTIVE_MARK

create_hybrid_weights(at, trans_width=3, buffer_width=0, weight_interpolation='distance_ramp')

print at.weight_region1[at.hybrid_mark == HYBRID_ACTIVE_MARK]
print at.weight_region1[at.hybrid_mark == HYBRID_TRANS_MARK]

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

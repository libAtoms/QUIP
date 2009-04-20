from numpy import *
from quippy import *

s = supercell(diamond(5.44, 14), 2, 2, 2)
s.set_cutoff(5.0)
s.calc_connect()

pot = Potential('IP SW',"/home/jk2/Code/QUIP/QUIP_Core/parameters/ip.parms.SW.xml")

ds = DynamicalSystem(s)
ds.rescale_velo(300.0)
ds.zero_momentum()

al = ds.run(pot, dt=1.0, n_steps=10, save_interval=1)

al.show()

raw_input()



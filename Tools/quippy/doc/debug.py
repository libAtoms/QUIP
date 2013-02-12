from quippy.structures import diamond, supercell
from atomeye import AtomEyeViewer

d = diamond(5.43, 14)
viewer = AtomEyeViewer(d)
viewer.capture('si8.png')

at = supercell(d, 2, 2, 2)
viewer.show(at)

viewer.change_bgcolor((0, 0, 0))
viewer.resize(400,300)
viewer.capture('si2x2x2.png')

raw_input()

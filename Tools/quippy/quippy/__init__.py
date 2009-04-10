"""pylibatoms package. Contains python bindings to the libAtoms F95 code 
<http://www.libatoms.org> and optionally also to QUIP. """

import libatoms
from libatoms import *

try:
    import quip
    from quip import *
except ImportError:
    pass

try:
    import atomeye
except ImportError:
    pass


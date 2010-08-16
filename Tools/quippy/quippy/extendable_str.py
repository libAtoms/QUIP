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

from quippy import FortranExtendable_str


class Extendable_str(FortranExtendable_str):

    def __init__(self, s=None):
        # Allow an Extendable_str to be created from a Python string or another Extendable_str
        if isinstance(s, Extendable_str):
            FortranExtendable_str.__init__(self, s)
        else:
            FortranExtendable_str.__init__(self)
            if s is not None:
                self.concat(str(s))

    @property
    def string(self):
        return ''.join(self.s[1,1:self.len])

    @string.setter
    def string(self, s):
        self.initialise()
        self += s

    def __iadd__(self, s):
        self.concat(s)

    def __str__(self):
        return self.string
    
    def __repr__(self):
        return 'Extendable_str("%s")' % str(self)

    def __len__(self):
        return self.len

    def __getitem__(self, i):
        return self.s[1,i]
    

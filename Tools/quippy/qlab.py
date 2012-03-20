from quippy import *
from atomeye import *
import numpy as np
from quippy.cinoutput import CInOutputReader

class AtomsViewer(AtomsReader, AtomEyeViewer):

    def __init__(self, source, **kwargs):
        AtomsReader.__init__(self, source, **kwargs)
        AtomEyeViewer.__init__(self, self)

    def show(self, **showargs):
        AtomEyeViewer.show(self, self, **showargs)

    def clip(self):
        indices = self.get_visible()
        if isinstance(self.reader, CInOutputReader):
            self.reader.source.indices = indices
        mask = fzeros(len(self.current_atoms), dtype=np.bool)
        mask[:] = True
        mask[indices] = False
        self._unclipped_atoms = self.current_atoms
        self.current_atoms.remove_atoms(mask=mask)
        self.show()

    def unclip(self):
        if hasattr(self, '_unclipped_atoms'):
            self.current_atoms = self._unclipped_atoms
        if isinstance(self.reader, CInOutputReader):
            self.reader.source.indices = None
        self.show()


def at(source):
    return Atoms(source)

def ar(source):
    return AtomsReader(source)

def al(source):
    return AtomsList(source)

def av(source):
    return AtomsViewer(source)

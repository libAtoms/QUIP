from quippy import *
from atomeye import *

class AtomsViewer(AtomsReader, AtomEyeViewer):

    def __init__(self, source, **kwargs):
        AtomsReader.__init__(self, source, **kwargs)
        AtomEyeViewer.__init__(self, self)

    def show(self, **showargs):
        AtomEyeViewer.show(self, self, **showargs)

    def clip(self):
        self.reader.indices = self.get_visible()


def at(source):
    return Atoms(source)

def ar(source):
    return AtomsReader(source)

def al(source):
    return AtomsList(source)

def av(source):
    return AtomsViewer(source)

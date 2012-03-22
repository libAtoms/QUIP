from numpy import *
from pylab import *
from quippy import *
from atomeye import *

from quippy.cinoutput import CInOutputReader

class AtomsViewer(AtomsReader, AtomEyeViewer):

    def __init__(self, source, cache=True, **kwargs):
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, **kwargs)
        AtomEyeViewer.__init__(self, self)

    def show(self, **showargs):
        AtomEyeViewer.show(self, self, **showargs)

    def clip_visible(self):
        self.redraw()
        self.wait()
        indices = self.get_visible()
        print 'Clipping view to include %d atoms' % len(indices)

        clipped_atoms = self.current_atoms.select(list=indices)
        self.atoms[self.frame] = clipped_atoms
        self.redraw()

        if isinstance(self.reader, CInOutputReader):
            self.reader.source.indices = indices

def at(source, **kwargs):
    return Atoms(source, **kwargs)

def ar(source, **kwargs):
    return AtomsReader(source, **kwargs)

def al(source, **kwargs):
    return AtomsList(source, **kwargs)

def av(source, **kwargs):
    return AtomsViewer(source, **kwargs)

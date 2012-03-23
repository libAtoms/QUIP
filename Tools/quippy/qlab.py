from numpy import *
from pylab import *
from quippy import *
from atomeye import *

from quippy.cinoutput import CInOutputReader

class AtomEyeViewerMixin(AtomEyeViewer):
    """Contains quippy-specific extensions to AtomEyeViewer"""
    
    def __init__(self):
        AtomEyeViewer.__init__(self, self,
                               enter_hook=self.enter_hook,
                               exit_hook=self.exit_hook,
                               click_hook=self.click_hook,
                               redraw_hook=self.redraw_hook)
        self.selection_mark = None

    def show(self, *args, **kwargs):
        AtomEyeViewer.show(self, self, *args, **kwargs)
        
    def enter_hook(self, atoms):
        def hook(at):
            self.redraw()
        atoms.redraw_hook = hook
        atoms.add_hook(atoms.redraw_hook)

    def exit_hook(self, atoms):
        atoms.remove_hook(atoms.redraw_hook)

    def click_hook(self, atoms, idx):
        if self.selection_mark is not None:
            self.selection_mark[idx] = self.selection_value
            print idx,
            sys.stdout.flush()
            self.redraw()
        elif self.verbose:
            print
            atoms.print_atom(idx)
            sys.stdout.flush()

    def redraw_hook(self, atoms):
        for (key, value) in atoms.params.iteritems():
            print '%-20s = %r' % (key, value)
    
    def clip_visible(self, orig_index=True):
        self.redraw()
        self.wait()
        indices = self.get_visible()
        print 'Clipping view to include %d atoms' % len(indices)

        mask = fzeros(self.current_atoms.n, dtype=bool)
        mask[:] = True
        mask[indices] = False
        self.current_atoms.remove_atoms(mask=mask)

        if orig_index:
            self.current_atoms.add_property('orig_index',
                                            logical_not(mask).nonzero()[0],
                                            overwrite=True)
        self.redraw()

        if hasattr(self, 'reader') and isinstance(self.reader, CInOutputReader):
            self.reader.source.indices = indices

    def select_atoms(self, reset=True, markname='selection_mark', value=True):
        if reset:
            self.current_atoms.add_property(markname, False, overwrite=True)
        self.selection_mark = getattr(self.current_atoms, markname)
        self.selection_value = value
        self.show(markname)
        saved_verbose = self.verbose
        self.verbose = False
        print 'Click to select atoms. Press ENTER to finish.'
        print 'indices = [',
        raw_input()
        print ']'
        indices = list(self.selection_mark.nonzero()[0])
        self.selection_mark = None
        self.verbose = saved_verbose
        return indices

        

class AtomsViewer(Atoms, AtomEyeViewerMixin):
    def __init__(self, source, **kwargs):
        Atoms.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self)

class TrajectoryViewer(AtomsReader, AtomEyeViewerMixin):
    def __init__(self, source, cache=True, **kwargs):
        if cache:
            total_mem, free_mem = mem_info()
            kwargs['cache_mem_limit'] = 0.5*free_mem
        AtomsReader.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self)

class TrajectoryEditor(AtomsList, AtomEyeViewerMixin):
    def __init__(self, source, **kwargs):
        AtomsList.__init__(self, source, **kwargs)
        AtomEyeViewerMixin.__init__(self)
            

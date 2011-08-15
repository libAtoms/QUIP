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

from quippy import AtomsWriters
import os
import atomeye

class AtomEyeWriter(object):
    "Write atoms to image file (png/eps/jpg) using AtomEye"

    def __init__(self, image, width=None, height=None, aspect=0.75, shift=True,
                 commands=None, script=None, nowindow=True, number_frames=False):
        self.image = image
        self.frame = None
        self.width = width
        self.height = height
        self.aspect = aspect
        self.shift = shift
        self.commands = commands
        self.script = script
        self.nowindow = nowindow
        self.first_config = True

    def write(self, at, centre=None, frame=None, *showargs, **showkwargs):
        if self.first_config:
            # If it's the first time, create an AtomEye viewer
            self.view = atomeye.AtomEyeView(nowindow=self.nowindow, echo=True, block=True)

        self.view.show(at, *showargs, **showkwargs)

        if self.first_config:
            self.first_config = False

            if self.commands is not None:
                for command in self.commands:
                    self.view.run_command(command)

            if self.script is not None:
                self.view.run_script(self.script)

            if self.width is not None or self.height is not None:
                if self.width  is None: self.width = int(self.height/self.aspect)
                if self.height is None: self.height = int(self.width*self.aspect)
                self.view.resize(self.width, self.height)

            if self.shift:
                self.view.xtal_origin_goto([.5, .5, .5])

        if centre is not None:
            if isinstance(centre, int):
                self.view.run_command('set n->anchor %d' % centre)
            else:
                self.view.run_command('set n->anchor -1')
                self.view.run_command('set n->hook %f %f %f' % tuple(centre + np.diag(at.lattice)/2.))
            self.view.look_at_the_anchor()

        if frame is not None:
            self.view.capture('%s%05d%s' % (os.path.splitext(self.image)[0], frame, os.path.splitext(self.image)[1]))
        else:
            self.view.capture(self.image)

    def close(self):
        pass

AtomsWriters['eps'] = AtomsWriters['png'] = AtomsWriters['jpg'] = AtomEyeWriter

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
from farray import *
import sys

class PovrayWriter(object):

    def __init__(self, f, center=(0.0,0.0,0.0), camerapos=(0.0,0.0,-100),
               lights=(0.0,0.0,-100), colour=None, rotate=None, radius=1.0,
               skip_hydrogens=False):

        if type(f) == type(''):
            if f == 'stdout':
                self.f = sys.stdout
                self.opened = False
            else:
                self.opened = True
                self.f = open(f,'w')
        else:
            self.opened = False
            self.f = f

        self.center = center
        self.camerapos = camerapos
        if len(lights) == 3:
            lights = (lights,)
        self.lights = lights
        self.colour = colour
        self.rotate = rotate
        self.radius = radius
        self.skip_hydrogens = skip_hydrogens

    def write(self, at):
        self.f.write('''#include "colors.inc"

camera
{
  right x*400/600
  location <%f,%f,%f>
  look_at <0,0,0>
}

background
{
  colour White
}
''' % tuple(self.camerapos))


        for light in self.lights:
            self.f.write('''
light_source
{
  <%f,%f,%f>
  colour White
}
''' % tuple(light))

        center = farray(self.center)
        self.f.write('union {\n')
        for i in frange(at.n):
            if self.skip_hydrogens and str(at.species[i]) == 'H': continue
            p = at.pos[i] - center
            if self.colour is not None:
                c = self.colour[i]
            else:
                c = farray((1.,1.,1.))

            self.f.write('sphere { <%f,%f,%f>, %r pigment {color <%f,%f,%f> } finish { phong .8 } }\n' %
                    (p[1],-p[2],p[3],self.radius,c[1],c[2],c[3]))

        if self.rotate is not None:
            self.f.write('rotate <%f,%f,%f>\n' % tuple(self.rotate))
        self.f.write('}')

    def close(self):
        if self.opened:  self.f.close()


AtomsWriters['pov'] = PovrayWriter

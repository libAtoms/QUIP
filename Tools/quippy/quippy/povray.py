from quippy import atoms_writer
from farray import *

@atoms_writer('pov')
def write_pov(f, center=(0.0,0.0,0.0), camerapos=(0.0,0.0,-100),
              lights=(0.0,0.0,-100), colour=None, rotate=None, radius=1.0,
              skip_hydrogens=False):

   at = yield None

   if type(f) == type(''):
      f = open(f,'w')

   f.write('''#include "colors.inc"

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
''' % tuple(camerapos))

   if len(lights) == 3:
      lights = (lights,)

   for light in lights:
      f.write('''
light_source
{
  <%f,%f,%f>
  colour White
}
''' % tuple(light))

   center = farray(center)
   f.write('union {\n')
   for i in frange(at.n):
      if skip_hydrogens and str(at.species[i]) == 'H': continue
      p = at.pos[i] - center
      if colour is not None:
         c = colour[i]
      else:
         c = farray((1.,1.,1.))
   
      f.write('sphere { <%f,%f,%f>, %r pigment {color <%f,%f,%f> } finish { phong .8 } }\n' %
              (p[1],-p[2],p[3],radius,c[1],c[2],c[3]))

   if rotate is not None:
      f.write('rotate <%f,%f,%f>\n' % tuple(rotate))
   f.write('}')

   yield None

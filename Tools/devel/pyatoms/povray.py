from numpy import *
from atoms import *


def write_pov(at, f, center=(0,0,0.0,0.0), camerapos=(0.0,0.0,-100),
              lights=(0.0,0.0,-100), colour=None, rotate=None, radius=1.0):
   opened = False
   if type(f) == type(''):
      f = open(f,'w')
      opened = True

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

   f.write('union {\n')
   for i in range(at.n):
      if at.species[i].startswith('H'): continue
      p = at.pos[i] - center
      if colour is not None:
         c = colour[i]
      else:
         c = array((1.,1.,1.))
   
      f.write('sphere { <%f,%f,%f>, %r pigment {color <%f,%f,%f> } finish { phong .8 } }\n' %
              (p[0],-p[1],p[2],radius,c[0],c[1],c[2]))

   if rotate is not None:
      f.write('rotate <%f,%f,%f>\n' % tuple(rotate))
   f.write('}')


   if opened:
      f.close()

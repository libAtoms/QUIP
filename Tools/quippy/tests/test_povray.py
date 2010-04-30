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

from quippy import *
import unittest
from quippytest import *

class TestPovray(QuippyTestCase):

   def setUp(self):
      from StringIO import StringIO
      self.s = StringIO()
      self.at = diamond(5.44, 14)

      self.pov_ref = """#include "colors.inc"

camera
{
  right x*400/600
  location <0.000000,0.000000,-100.000000>
  look_at <0,0,0>
}

background
{
  colour White
}

light_source
{
  <0.000000,0.000000,-100.000000>
  colour White
}
union {
sphere { <0.000000,-0.000000,0.000000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
sphere { <1.360000,-1.360000,1.360000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
sphere { <2.720000,-2.720000,0.000000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
sphere { <4.080000,-4.080000,1.360000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
sphere { <2.720000,-0.000000,2.720000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
sphere { <4.080000,-1.360000,4.080000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
sphere { <0.000000,-2.720000,2.720000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
sphere { <1.360000,-4.080000,4.080000>, 1.0 pigment {color <1.000000,1.000000,1.000000> } finish { phong .8 } }
}"""

   def testpov(self):
      self.at.write(self.s, format='pov')
      self.assertEqual(self.pov_ref.replace('-0.00','0.00'), self.s.getvalue().replace('-0.00','0.00'))


if __name__ == '__main__':
   unittest.main()





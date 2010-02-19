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





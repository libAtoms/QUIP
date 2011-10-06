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

from quippy.farray import *
import numpy
from quippytest import *
import unittest

class TestFortranArray(QuippyTestCase):

   def setUp(self):
      self.z0 =fzeros(0)
      self.z1 = fzeros(1)
      self.z2 = fzeros((2,2))
      self.f = farray((1,2,3))
      self.f2 = farray(((1,2,3),(4,5,6))).T
      self.na =  numpy.array(self.f2)
      self.s = farray(('s1','s2','s3'))

   def testfrange(self):
      self.assertEqual(list(frange(3)), [1,2,3])
      self.assertEqual(list(frange(3,6,2)), [3, 5])

   def testfenumerate(self):
      self.assertEqual(list(fenumerate(frange(3))), [(1,1), (2,2), (3,3)])

   def testfzeros0(self):
      self.assertEqual(self.z0.shape, (0,))

   def testfzeros1(self):
      self.assertEqual(self.z1.shape, (1,))
      self.assertEqual(self.z1[1], 0)

   def testfzeros2(self):
      self.assertEqual(self.z2.shape, (2,2))
      self.assertAlmostEqual(self.z2[1,1], 0)
      self.assertAlmostEqual(self.z2[1,2], 0)
      self.assertAlmostEqual(self.z2[2,1], 0)
      self.assertAlmostEqual(self.z2[2,2], 0)

   def testfindentity(self):
      i2 = fidentity(2)
      self.assertEqual(i2.shape, (2,2))
      self.assertAlmostEqual(i2[1,1], 1)
      self.assertAlmostEqual(i2[1,2], 0)
      self.assertAlmostEqual(i2[2,1], 0)
      self.assertAlmostEqual(i2[2,2], 1)
      
   def testfvar(self):
      fvar('z')
      self.assertEqual(z.shape, ())
      self.assertEqual(z, 0.0)

      fvar('abc')
      self.assertEqual(a.shape, ())
      self.assertEqual(a, 0.0)
      self.assertEqual(b.shape, ())
      self.assertEqual(b, 0.0)
      self.assertEqual(c.shape, ())
      self.assertEqual(c, 0.0)
      
   def testshape(self):
      self.assertEqual(self.f.shape, (3,))

   def testlen(self):
      self.assertEqual(len(self.f), 3)

   def testindexzero(self):
      self.assertRaises(IndexError, self.f.__getitem__, 0)

   def testindex(self):
      self.assert_(isinstance(self.f[1], numpy.int))
      self.assertEqual(self.f[1], 1)
      self.assertEqual(self.f[2], 2)
      self.assertEqual(self.f[3], 3)

   def testsetzero(self):
      self.assertRaises(IndexError, self.f.__setitem__, 0, 1)

   def testsetitem(self):
      self.f[1] = -1
      self.assertEqual(self.f[1], -1)
      self.assertEqual(self.f[2], 2)
      self.assertEqual(self.f[3], 3)

   def testnegindex(self):
      self.assert_(isinstance(self.f[-1], numpy.int))
      self.assertEqual(self.f[-1], 3)
      self.assertEqual(self.f[-2], 2)
      self.assertEqual(self.f[-3], 1)

   def testslice(self):
      sl = self.f[1:]
      self.assert_(isinstance(sl, FortranArray))
      self.assertEqual(sl.shape, (3,))
      self.assertEqual(sl[1], 1)
      self.assertEqual(sl[2], 2)
      self.assertEqual(sl[3], 3)

   def testsetslice(self):
      self.f[:] = 0
      self.assertEqual(self.f[1], 0)
      self.assertEqual(self.f[2], 0)
      self.assertEqual(self.f[3], 0)

   def testlist(self):
      self.assertEqual(list(self.f), [1,2,3])

   def testslice2(self):
      sl = self.f[2:]
      self.assert_(isinstance(sl, FortranArray))
      self.assertEqual(sl.shape, (2,))
      self.assertEqual(sl[1], 2)
      self.assertEqual(sl[2], 3)

   def testnegslice(self):
      sl = self.f[:-1]
      self.assert_(isinstance(sl, FortranArray))
      self.assertEqual(sl.shape, (2,))
      self.assertEqual(sl[1], 1)
      self.assertEqual(sl[2], 2)

   def testnegslice2(self):
      sl = self.f[2:-1]
      self.assert_(isinstance(sl, FortranArray))
      self.assertEqual(sl.shape, (1,))
      self.assertEqual(sl[1], 2)

# Negative slices don't work correctly -- this test fails
# To workaround this, would somehow need to remove __setslice__ and __getslice__ from MRO
##    def testnegslice3(self):
##       sl = self.f[-2:]
##       self.assertEqual(list(sl), [2,3])

   def testrows(self):
      self.assertEqual(farray(0).rows.next(), 0)
      self.assertEqual(list(self.f.rows), [1,2,3])

   def testcols(self):
      self.assertEqual(farray(0).cols.next(), 0)
      self.assertEqual(list(self.f.cols), [1,2,3])

   def testextsetlist(self):
      self.f[[1,2]] = 0
      self.assertEqual(list(self.f), [0,0,3])

   def testextsetfarray(self):
      self.f[farray([1])] = 0
      self.assertEqual(list(self.f), [0,2,3])

   def testextsetlistbool(self):
      self.f[[True,False,False]] = 0
      self.assertEqual(list(self.f), [0,2,3])

   def testextsetlistbool(self):
      self.f[farray([True,False,False])] = 0
      self.assertEqual(list(self.f), [0,2,3])

   def test2dshape(self):
      self.assertEqual(self.f2.shape, (3,2))

   def test2dindexzero(self):
      self.assertRaises(IndexError, self.f2.__getitem__, 0)
      self.assertRaises(IndexError, self.f2.__getitem__, (0,0))
   
   def test2dindex1(self):
      self.assert_(isinstance(self.f2[1], FortranArray))
      self.assertEqual(list(self.f2[1]), [1,2,3])
      self.assertEqual(list(self.f2[1]), list(self.f2[:,1]))
      self.assertEqual(list(self.f2[1]), list(self.f2[...,1]))

   def test2dindex2(self):
      self.assert_(isinstance(self.f2[1,1], numpy.int))
      self.assertEqual(self.f2[1,1], 1)
      self.assertEqual(self.f2[2,1], 2)
      self.assertEqual(self.f2[3,1], 3)
      self.assertEqual(self.f2[1,2], 4)
      self.assertEqual(self.f2[2,2], 5)
      self.assertEqual(self.f2[3,2], 6)

   def test2dsetitem(self):
      self.f2[1,1] = 0
      self.assertEqual(list(self.f2[1]), [0,2,3])
      self.assertEqual(list(self.f2[2]), [4,5,6])

   def test2dsetslice(self):
      self.f2[:,1] = 0
      self.assertEqual(list(self.f2[1]), [0,0,0])
      self.assertEqual(list(self.f2[2]), [4,5,6])

   def test2dsetslice2(self):
      self.f2[1,:] = 0
      self.assertEqual(list(self.f2[1]), [0,2,3])
      self.assertEqual(list(self.f2[2]), [0,5,6])

   def test2dsetslice3(self):
      self.f2[:,:] = 0
      self.assertEqual(list(self.f2[1]), [0,0,0])
      self.assertEqual(list(self.f2[2]), [0,0,0])

   def test2dsetslice4(self):
      self.f2[:] = 0
      self.assertEqual(list(self.f2[1]), [0,0,0])
      self.assertEqual(list(self.f2[2]), [0,0,0])

   def test2dsetslice5(self):
      self.f2[1::2,1] = 0
      self.assertEqual(list(self.f2[1]), [0,2,0])
      self.assertEqual(list(self.f2[2]), [4,5,6])

   def test2dsetslice6(self):
      self.f2[1:] = 0
      self.assertEqual(list(self.f2[1]), [0,0,0])
      self.assertEqual(list(self.f2[2]), [0,0,0])

   def test2dsetellipsis(self):
      self.f2[...] = 0
      self.assertEqual(list(self.f2[1]), [0,0,0])
      self.assertEqual(list(self.f2[2]), [0,0,0])

   def test2dsetellipsis2(self):
      self.f2[...,1] = 0
      self.assertEqual(list(self.f2[1]), [0,0,0])
      self.assertEqual(list(self.f2[2]), [4,5,6])

   def test2dsetellipsis3(self):
      self.f2[1,...] = 0
      self.assertEqual(list(self.f2[1]), [0,2,3])
      self.assertEqual(list(self.f2[2]), [0,5,6])

   def test2dsetcol(self):
      self.f2[1] = 0
      self.assertEqual(list(self.f2[1]), [0,0,0])
      self.assertEqual(list(self.f2[2]), [4,5,6])         

   def test2dcols(self):
      col1 = list(self.f2.cols.next())
      self.assertEqual(col1, [1,2,3])
      self.assertEqual(col1, list(self.f2[1]))

   def test2drows(self):
      row1 = list(self.f2.rows.next())
      self.assertEqual(row1, [1,4])
      self.assertEqual(row1, list(self.f2[1,:]))

   def test2diter(self):
      col1 = iter(self.f2).next()
      self.assertEqual(list(col1), [1,2,3])

   def test2dfromseq(self):
      fa = FortranArray(((1,2,3),(4,5,6),(7,8,9))).T
      self.assert_(isinstance(fa, FortranArray))
      self.assertEqual(list(fa[1]), [1,2,3])
      self.assertEqual(list(fa[2]), [4,5,6])
      self.assertEqual(list(fa[3]), [7,8,9])
      self.assertEqual(list(fa[:,1]), [1,2,3])
      self.assertEqual(list(fa[:,2]), [4,5,6])
      self.assertEqual(list(fa[:,3]), [7,8,9])
      self.assertEqual(list(fa[1,:]), [1,4,7])
      self.assertEqual(list(fa[2,:]), [2,5,8])
      self.assertEqual(list(fa[3,:]), [3,6,9])

   def test2dextindexintlist(self):
      self.assert_((self.f2[...,[1]] == self.na[...,[0]]).all())
      self.assert_((self.f2[...,[1,2]] == self.na[...,[0,1]]).all())
      self.assert_((self.f2[...,[2,1],:] == self.na[...,[1,0],:]).all())
      self.assert_((self.f2[...,[2,1],:] == self.na[...,[1,0],:]).all())

   def test2dextindexintarray(self):
      self.assert_((self.f2[farray([1])] == self.na[...,numpy.array([0])]).all())
      self.assert_((self.f2[farray([1,2])] == self.na[...,numpy.array([0,1])]).all())
      self.assert_((self.f2[farray([2,1]),:] == self.na[...,numpy.array([1,0]),:]).all())

   def test2dextindexsetlist(self):
      self.f2[farray([1,2]),:] = 0
      self.assertEqual(list(self.f2[1]), [0,0,3])
      self.assertEqual(list(self.f2[2]), [0,0,6])

   def test2dextindexsetarray(self):
      self.f2[[1,2],:] = 0
      self.assertEqual(list(self.f2[1]), [0,0,3])
      self.assertEqual(list(self.f2[2]), [0,0,6])

   def test2dextindexintfail(self):
      self.assertRaises(ValueError, self.f2.__getitem__, farray(0))
      self.assertRaises(ValueError, self.f2.__getitem__, [[0]])
      self.assertRaises(ValueError, self.f2.__getitem__, farray('str'))
      self.assertRaises(ValueError, self.f2.__getitem__, [{}])
      
   def test2dextindexbool(self):
      bf = self.f2[self.f2 > 1]
      bn = self.na[self.na > 1]
      self.assert_((bf == bn).all())

   def test2dextindexbool2(self):
      self.f2[self.f2 > 1] = 0
      self.na[self.na > 1] = 0
      self.assert_((self.f2 == self.na).all())
      
   def testnonzero(self):
      self.assertEqual(list((self.f > 1).nonzero()[0]), [2,3])
      self.assert_((self.f2[(self.f2 > 1).nonzero()] == self.f2[self.f2 > 1]).all())

   def testcount(self):
      self.assertEqual((self.f2 > 1).sum(), 5)
      self.assertEqual((self.f2 == 1).sum(), 1)

   def testinteq(self):
      la = self.f2 == 1
      self.assert_(isinstance(la, FortranArray))
      self.assertEqual(list(la[1]), [True, False, False])
      self.assertEqual(list(la[2]), [False, False, False])

   def testintne(self):
      la2 = self.f2 != 1
      self.assert_(isinstance(la2, FortranArray))
      self.assertEqual(list(la2[1]), [False, True, True])
      self.assertEqual(list(la2[2]), [True, True, True])

   def testintge(self):
      la3 = self.f2 >= 1
      self.assert_(isinstance(la3, FortranArray))
      self.assertEqual(list(la3[1]), [True, True, True])
      self.assertEqual(list(la3[2]), [True, True, True])

   def testintle(self):
      la4 = self.f2 <= 1
      self.assert_(isinstance(la4, FortranArray))
      self.assertEqual(list(la4[1]), [True, False, False])
      self.assertEqual(list(la4[2]), [False, False, False])

   def testcmpstring(self):
      la = self.s == 's1'
      self.assert_(isinstance(la, FortranArray))
      self.assertEqual(list(la), [True, False, False])

      la2 = self.s != 's1'
      self.assert_(isinstance(la2, FortranArray))
      self.assertEqual(list(la2), [False, True, True])

   def testargminmax(self):
      self.assertEqual(self.f2.argmin(), 1)
      self.assertEqual(list(self.f2.argmin(1)), [1,1])
      self.assertEqual(list(self.f2.argmin(2)), [1,1,1])

      self.assertEqual(self.f2.argmax(), 6)
      self.assertEqual(list(self.f2.argmax(1)), [3,3])
      self.assertEqual(list(self.f2.argmax(2)), [2,2,2])

   def testargsort(self):
      self.assertEqual(list(self.f2.argsort()), [1,3,5,2,4,6])
      self.assertEqual(list(self.f2.argsort(1)[1]), [1,2,3])
      self.assertEqual(list(self.f2.argsort(2)[1]), [1,1,1])

   def testtake(self):
      self.assertEqual(list(self.f2.take((1,2,3))), [1,4,2])
      self.assertEqual(list(self.f2.take((1,2),axis=1)[1]), [1,2])

   def testput(self):
      self.f2.put((1, 4), (0, 0))
      self.assertEqual(list(self.f2[1]), [0, 2, 3])
      self.assertEqual(list(self.f2[2]), [4, 0, 6])

   def testreprstr(self):
      self.assertEqual(repr(self.f2),"""FortranArray([[1, 4],
              [2, 5],
              [3, 6]])""")
      self.assertEqual(str(self.f2),"""[[1 4]
 [2 5]
 [3 6]]""")

      self.assertEqual(str(self.s), "['s1' 's2' 's3']")
      self.assertEqual(str(farray('s')), 's')

      s2 = farray((('s','1'), ('s','2'))).T
      self.assertEqual([str(x) for x in s2], ["['s' '1']", "['s' '2']"])
      
   def testnorm(self):
      self.assertAlmostEqual(farray(1.0).norm(), 1.0)
      self.assertRaises(ValueError, self.f2.T.norm2)
      self.assertAlmostEqual(self.f2[1].norm2(), 14.0)
      n2 = self.f2.norm2()
      self.assertAlmostEqual(n2[1], 14.0)
      self.assertAlmostEqual(n2[2], 77.0)
      self.assert_((numpy.sqrt(n2) - self.f2.norm() < 1e-6).all())
      self.assertRaises(ValueError, fzeros((3,3,3)).norm2)

   def testall(self):
      self.assert_((self.f2 >= 1).all())
      self.assert_(not (self.f2 == 1).all())
      self.assertEqual(list((self.f2 >= 1).all(axis=1)), [True, True])
      self.assertEqual(list((self.f2 > 2).all(axis=2)), [False, False, True])      

   def testany(self):
      self.assert_((self.f2 >= 1).any())
      self.assert_((self.f2 == 1).any())
      self.assertEqual(list((self.f2 >= 1).any(axis=1)), [True, True])
      self.assertEqual(list((self.f2 > 2).any(axis=2)), [True, True, True])      
      
      
if __name__ == '__main__':
   unittest.main()

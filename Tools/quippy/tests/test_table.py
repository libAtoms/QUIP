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

import unittest, numpy
from quippytest import *

from quippy import Table, farray

class TestTable(QuippyTestCase):

   def testempty(self):
      t = Table()
      self.assertEqual(t.n, 0)

   def testsingleint(self):
      t = Table()
      t.append(1)
      t.append(2)
      t.append(3)
      self.assertEqual(list(t.int), [1,2,3])
      t.append([4,5,6])
      self.assertEqual(list(t.int), [1,2,3,4,5,6])
      self.assertArrayAlmostEqual(t.int_part(1,t.n), farray([1,2,3,4,5,6]))

   def testsinglereal(self):
      t = Table()
      t.append(1.0)
      t.append(2.0)
      t.append(3.0)
      self.assertArrayAlmostEqual(t.real[1,:], farray([1.0, 2.0, 3.0]))
      t.append(realpart=[4.0, 5.0, 6.0])
      self.assertArrayAlmostEqual(t.real[1,:], farray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))
      self.assertArrayAlmostEqual(t.real_part(column=1,n0=t.n), farray([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]))

   def testintreal(self):
      t = Table()
      t.append(1, 1.0)
      t.append(2, 2.0)
      t.append(3, 3.0)
      t.append(intpart=[4,5,6], realpart=[4.0, 5.0, 6.0])
      self.assertArrayAlmostEqual(t.int[1,:],  farray([1,2,3,4,5,6]))
      self.assertArrayAlmostEqual(t.int_part(column=1,n0=t.n),  farray([1,2,3,4,5,6]))
      self.assertArrayAlmostEqual(t.real[1,:], farray([1,2,3,4,5,6]))
      self.assertArrayAlmostEqual(t.real_part(column=1,n0=t.n), farray([1,2,3,4,5,6]))

   def testsinglestr(self):
      t = Table()
      t.append('str1')
      t.append('str2')
      t.append('str3')
      t.append('str4')
      self.assertEqual(list(t.str[:,1,:].stripstrings()), ['str1', 'str2', 'str3', 'str4'])

   def testsinglelogical(self):
      t = Table()
      t.append(False)
      t.append(False)
      t.append(True)
      t.append(False)
      t.append(False)
      t.append(False)
      self.assertEqual(list(t.logical[1,:]), [False, False, True, False, False, False])
      self.assertEqual(list(t.logical_part(1, t.n)), [False, False, True, False, False, False])

   def testmultintmultlog(self):
      t = Table()
      t.append(realpart=[1.0,2.0,3.0],logicalpart=[False,True,False])
      t.append(realpart_2d=numpy.reshape([4.0,5.0,6.0, 7.0,8.0,9.0],[3,2],order='F'),
               logicalpart_2d=numpy.reshape([False,False,False,True,True,True],[3,2],order='F'))
      self.assertArrayAlmostEqual(t.real, farray([[1.0,2.0,3.0],[4.0,5.0,6.0],[7.0,8.0,9.0]]).T)
      self.assertArrayAlmostEqual(t.logical, farray([[False,True,False],[False,False,False],[True,True,True]]).T)

   def testtableappend(self):
      t1 = Table(2,1,0,0)
      t1.append([1,2], 3.0)
      t1.append([4,5], 6.0)
      t1.append([7,8], 9.0)
      t2 = Table()
      t2.append(t1)
      t2.append(t1)
      self.assertEqual(list(t2.int[1,:]), [1,4,7,1,4,7])
      self.assertEqual(list(t2.int[2,:]), [2,5,8,2,5,8])
      self.assertArrayAlmostEqual(t2.real[1,:], farray([3.0,6.0,9.0,3.0,6.0,9.0]))

   def testtableappendstr(self):
      t1 = Table()
      t1.append('hello')
      t1.append('world')
      t2 = Table()
      t2.append(t1)
      t2.append(t1)
      self.assertEqual(list(t2.str[:,1,:].stripstrings()), ['hello', 'world', 'hello', 'world'])

   def testappendcolumn(self):
      t1 = Table(1,0,0,0)
      t1.append([1,2,3])
      t1.append_column([4,5,6],n_cols=2)
      self.assertEqual(list(t1.int[1,:]), [1,2,3])
      self.assertEqual(list(t1.int[2,:]), [4,5,6])
      self.assertEqual(list(t1.int[3,:]), [4,5,6])

   def testsubtable(self):
      t1 = Table(1,0,0,0) 
      t1.append([1,2,3])
      t2 = t1.subtable([1,3])
      self.assertEqual(list(t2.int), [1,3])

   def testremovecolumns(self):
      # Table.removecolumns assigns to 'this', so it's not a proper method
      # hence it doesn't work with quippy
      pass
      #t1 = Table(1,0,0,0)
      #t1.append([1,2,3])
      #t1.append_column([4,5,6])
      #t1.remove_columns(int_col_min=1)
      #self.assertEqual(list(t1.int), [1,2,3])


   def testinsert(self):
      t = Table(1,0,0,0)
      t.append([1,2,3])
      t.insert(3, intpart=[4])
      self.assertEqual(list(t.int), [1,2,4,3])

   def testfind(self):
      t = Table(1,0,0,0)
      t.append([1,2,3,4,5,6,7,8])
      self.assertEqual(t.find(5), 5)

   def testsort(self):
      t = Table(1,0,0,0)
      t.append([4,3,7,1,6,2,8,5])
      t.sort()
      self.assertEqual(list(t.int), [1,2,3,4,5,6,7,8])

   def testsearch(self):
      t = Table(1,0,0,0)
      t.append([1,2,3,4,5,6,7,8])
      self.assertEqual(t.search(5), 5)

   def testdelete(self):
      t = Table(1,0,0,0)
      t.append([2,4,6,8,10])
      t.delete(1, keep_order=False)
      self.assertEqual(list(t.int), [10,4,6,8])
      t.delete([4], keep_order=True)
      self.assertEqual(list(t.int), [10,6,8])

   def testrecord_delete_multiple(self):
      t = Table(1,0,0,0)
      t.append([2,4,6,8,10])
      t.record_delete_multiple([1,2])
      self.assertEqual(list(t.int), [10,8,6])

   def testwipe(self): 
      t = Table(1,0,0,0)
      t.append([2,4,6,8,10])
      t.wipe()
      self.assertEqual(t.n, 0)
      self.assertEqual(list(t.int), [])
      self.assertEqual(t.intsize, 1)
      
      
if __name__ == '__main__':
   unittest.main()

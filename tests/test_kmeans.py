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
import unittest, quippy, numpy
from quippytest import *

class TestKMeans(QuippyTestCase):

    def setUp(self):
        pass

    def test_initial_means(self):
        # deterministic test where initial means are chosen to yield three clusters
        k = 3
        n = 7
        m = 1
        
        data = fzeros((n,m))
        data[:,1] = [1.0,1.1,1.2,3.0,5.0,5.1,5.2]
        means = fzeros((k,m))
        means[1,:] = 1.1
        means[2,:] = 3
        means[3,:] = 5.1
        assign = fzeros(n,dtype=numpy.int32)
        
        kmeans(data,k,means,assign)
        self.assertEqual(list(assign), [1,1,1,2,3,3,3])

    def test_zero_means(self):
        # deterministic test where initial means are set to zero
        # results in two clusters.
        k = 3
        n = 7
        m = 1
        
        data = fzeros((n,m))
        data[:,1] = [1.0,1.1,1.2,3.0,5.0,5.1,5.2]
        means = fzeros((k,m))
        assign = fzeros(n,dtype=numpy.int32)
        err = farray(0.0)
        
        kmeans(data,k,means,assign,err)
        self.assertArrayAlmostEqual(means.T, [[4.575, 1.1, 0.0]])
        self.assertEqual(list(assign), [2,2,2,1,1,1,1])
        self.assertAlmostEqual(err, 0.96589272460937403)


    def test_random_partition(self):
        # probabilistic test based on random partitions
        # most probable outcome is three clusters - we do 1000 trials
        # and check it comes up at least 75% of the time
        
        k = 3
        n = 7
        m = 1
        
        data = fzeros((n,m))
        data[:,1] = [1.0,1.1,1.2,3.0,5.0,5.1,5.2]
        means = fzeros((k,m))
        assign = fzeros(n,dtype=numpy.int32)

        count = 0
        for trial in range(1000):
            kmeans(data,k,means,assign, initialisation='random_partition')
            order = list(numpy.argsort(means[:,1]))
            reordered_assign =  [order.index(assign[i])+1 for i in frange(n)]
            if reordered_assign == [1,1,1,2,3,3,3]: count += 1

        # Check that most probable outcome is three clusters
        self.assert_(count/1000.0 > 0.75)

    def test_random_means(self):
        # probabilistic test based on random initial means
        # most probable outcome is two clusters - we do 1000 trials
        # and check it comes up at least 90% of the time
        
        k = 2
        n = 6
        m = 1
        
        data = fzeros((n,m))
        data[:,1] = [1.0,1.1,1.2,5.0,5.1,5.2]
        means = fzeros((k,m))
        assign = fzeros(n,dtype=numpy.int32)

        count = 0
        for trial in range(1000):
            kmeans(data,k,means,assign, initialisation='random_means')
            order = list(numpy.argsort(means[:,1]))
            reordered_assign =  [order.index(assign[i])+1 for i in frange(n)]
            if reordered_assign == [1, 1, 1, 2, 2, 2]: count += 1

        # Check that most probable outcome is two clusters
        self.assert_(count/1000.0 > 0.90)

        

if __name__ == '__main__':
    unittest.main()

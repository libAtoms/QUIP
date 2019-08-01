# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2019
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

import unittest
import logging
import numpy as np
import quippy


# def string_to_array(s):
#     return loadtxt(StringIO(s)).T


class QuippyTestCase(unittest.TestCase):

    def assertDictionariesEqual(self, d1, d2, skip_keys=[], ignore_case=True):

        def lower_if_ignore_case(k):
            if ignore_case:
                return k.lower()
            else:
                return k

        d1 = dict([(lower_if_ignore_case(k), v) for (k, v) in d1.iteritems() if k not in skip_keys])
        d2 = dict([(lower_if_ignore_case(k), v) for (k, v) in d2.iteritems() if k not in skip_keys])

        if sorted(d1.keys()) != sorted(d2.keys()):
            self.fail('Dictionaries differ: d1.keys() (%r) != d2.keys() (%r)' % (d1.keys(), d2.keys()))
        for key in d1:
            v1, v2 = d1[key], d2[key]
            if not np.array_equal(v1, v2):
                self.fail('Dictionaries differ: key=%s value1=%r value2=%r' % (key, v1, v2))

    def assertEqual(self, first, second, msg=None):
        if first == second:
            return
        # Repeat comparison with debug-level logging
        import logging
        level = logging.root.level
        logging.root.setLevel(logging.DEBUG)
        first == second
        logging.root.setLevel(level)
        self.fail('%s != %s' % (first, second))

    def assertArrayAlmostEqual(self, first, second, tol=1e-7):
        first = np.array(first)
        second = np.array(second)
        self.assertEqual(first.shape, second.shape)

        if np.isnan(first).any():
            self.fail('Not a number (NaN) found in first array')
        if np.isnan(second).any():
            self.fail('Not a number (NaN) found in second array')

        absdiff = abs(first - second)
        if np.max(absdiff) > tol:
            print('First array: \n', first)
            print('\n \n Second array: \n', second)
            print('\n \n Abs Difference: \n', absdiff)
            self.fail('Maximum abs difference between array elements is %e at location %r' % (np.max(absdiff),
                                                                                              np.argmax(absdiff)))


def skip(f):
    def g(self):
        logging.warning('skipping test %s' % f.__name__)

    return g


def profile(f):
    import cProfile, pstats, functools

    @functools.wraps(f)
    def g(self):
        cProfile.runctx('f(self)', globals(), locals(), f.__name__ + '.profile')
        p = pstats.Stats(f.__name__ + '.profile')
        p.strip_dirs().sort_stats('cumulative').print_stats()
        return f(self)

    return g

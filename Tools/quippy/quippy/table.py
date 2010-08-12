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

from quippy import FortranTable

class Table(FortranTable):

   __doc__ = FortranTable.__doc__

   def __repr__(self):
      return ('Table(n=%d,intsize=%d,realsize=%d,strsize=%d,logicalsize=%d)' %
              (self.n, self.intsize, self.realsize, self.strsize, self.logicalsize))

   def __str__(self):
      return repr(self)

   def copy(self):
      t = Table(self.intsize, self.realsize, self.strsize, self.logicalsize, self.n)
      t.append(blank_rows=self.n)
      if self.intsize != 0: t.int[...] = self.int[...]
      if self.realsize != 0: t.real[...] = self.real[...]
      if self.strsize != 0: t.str[...] = self.str[...]
      if self.logicalsize != 0: t.logical[...] = self.logical[...]
      return t

   def __eq__(self, other):
      return self.equal(other)

   def __ne__(self, other):
      return not self.equal(other)

   def equal(self, other):
      tol = 1e-8
      
      for t1, t2 in zip((self.n,  self.intsize,  self.realsize,  self.strsize,  self.logicalsize),
                        (other.n, other.intsize, other.realsize, other.strsize, other.logicalsize)):
         if t1 != t2: return False

      for n, a1, a2, in zip((self.intsize, self.realsize, self.strsize, self.logicalsize),
                            (self.int,  self.real,  self.str,  self.logical),
                            (other.int, other.real, other.str, other.logical)):

         if n == 0: continue
         try:
            if abs(a1 - a2).max() > tol: return False
         except TypeError:
            if (a1 != a2).any(): return False

      return True


   def _get_array_shape(self, name):
       if name in ('int','real','logical'):
           return (slice(None),slice(1,self.n))
       elif name == 'str':
           return (slice(None),slice(None),slice(1,self.n))
       else:
          return None



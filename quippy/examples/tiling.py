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

import operator
from math import sqrt

def primeListSofE(n):
   sieve = [True for j in xrange(2,n+1)]
   for j in xrange(2,int(sqrt(n))+1):
      i = j-2
      if sieve[i]:
         for k in range(j*j,n+1,j):
            sieve[k-2] = False
 
   return [j for j in xrange(2,n+1) if sieve[j-2]]
 
def factor(n) :
    """Returns all prime factors of n, using trial division by prime
    numbers only. Returns a list of (possibly repeating) prime factors
    """
    ret =[]
    nn = n
    maxFactor = int(n**0.5)
    primes = primeListSofE(maxFactor)
    for p in primes :
       while nn % p == 0 :
          nn //= p
          ret.append(p)
       if nn == 1 :
          break
    if nn != 1 :
        ret.append(nn)

    return ret


def decomp(nproc, na,nb,nc):

   prime_factors = factor(nproc)

   n_dim = sum([n != 1 for n in (na, nb, nc)])

   # Combine into n_dim similar sized factors
   while len(prime_factors) > n_dim:
      prime_factors = [prime_factors[0]*prime_factors[1]] + prime_factors[2:]

   while len(prime_factors) < n_dim:
      prime_factors.append(1)

   return prime_factors

   p1, p2, p3 = sorted(prime_factors, reverse=True)

   n1, n2, n3 = sorted([na,nb,nc], reverse=True)

   return (p1, p2, p3), (n1/p1, n2/p2, n3/p3)


def tile(cell,proc):
   cell_per_proc = cell // proc
   rem = cell % proc

   n = 1
   for i in frange(1, cell-rem, cell_per_proc):
      start = i
      if n == proc:
         end = cell
      else:
         end = i+cell_per_proc-1
      n += 1
      print start, end


def tile2d(cellx, celly, procx, procy):
   res = fzeros((cellx, celly), dtype=int)

   for i in range(1, cellx+1):
      for j in range(1, celly+1):
         px = floor(float(procx)*(i-1)/cellx)
         py = floor(float(procy)*(j-1)/celly)
         res[i,j] = px*procy + py
   

   return res

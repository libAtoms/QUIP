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

"""Utility functions which will be imported into top-level quippy namespace"""

def args_str(D):
   """Construct args string from file, string or mapping object"""
   from dictmixin import PuPyDictionary
   return str(PuPyDictionary(D))

def parse_slice(S):
   """Parse string containing slice in form [start]:[stop]:[range] and return slice instance."""

   class SliceParser(object):
      def __getitem__(self, idx):
         return idx

   return eval('SliceParser()[%s]' % S)

def parse_comma_colon_list(L):
   """Parse a comma or colon seperated string into a list, converting each entry to lower-case."""
   if ':' in L:
      L = L.split(':')
   elif ',' in L:
      L = L.split(',')
   else:
      L = [L]

   return [k.lower() for k in L]

def loadstring(s): 
   import StringIO
   from numpy import loadtxt
   return loadtxt(StringIO.StringIO(s.replace('[','').replace(']','')))

try:
   from pylab import plot

   from farray import convert_farray_to_ndarray
   plot = convert_farray_to_ndarray(plot)

except ImportError:
   pass

        
    

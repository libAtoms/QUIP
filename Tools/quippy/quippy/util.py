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


try:
   from pylab import plot

   from farray import convert_farray_to_ndarray
   plot = convert_farray_to_ndarray(plot)

except ImportError:
   pass

        
    

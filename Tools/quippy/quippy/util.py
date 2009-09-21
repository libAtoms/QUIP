"""Utility functions which will be imported into top-level quippy namespace"""

def args_str(D):
   """Construct args string from file, string or mapping object"""
   from dictmixin import PuPyDictionary
   return str(PuPyDictionary(D))



try:
   from pylab import plot

   from farray import convert_farray_to_ndarray
   plot = convert_farray_to_ndarray(plot)

except ImportError:
   pass

        
    

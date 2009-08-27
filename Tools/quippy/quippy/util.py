"""Utility functions which will be imported into top-level quippy namespace"""

def args_str(D):
   """Construct args string from file, string or mapping object"""
   from dictmixin import PuPyDictionary
   return str(PuPyDictionary(D))



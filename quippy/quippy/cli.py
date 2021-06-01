import os
import sys
import subprocess as sp
import quippy

def gap_fit():
    path = quippy.__path__[0]
    command = os.path.join(path, 'gap_fit')
    sp.call([command] + sys.argv)

def quip():
    path = quippy.__path__[0]
    command = os.path.join(path, 'quip')
    sp.call([command] + sys.argv)    

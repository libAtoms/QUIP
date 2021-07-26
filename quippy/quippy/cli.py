import os
import sys
import subprocess as sp
import quippy
import argparse

def gap_fit():
    path = quippy.__path__[0]
    command = os.path.join(path, 'gap_fit')
    sp.call([command] + sys.argv[1:])

def quip():
    path = quippy.__path__[0]
    command = os.path.join(path, 'quip')
    sp.call([command] + sys.argv[1:]) 

def quip_config():
    parser = argparse.ArgumentParser(description='Configuration tool for QUIP')
    parser.add_argument('--libs', action='store_true', help="Arguments to link to libquip")
    args = parser.parse_args()
    
    if args.libs:        
        libdir = quippy.__path__[0]
        print(f'-L{libdir} -lquip')
    
#!/usr/bin/env python

from __future__ import print_function, unicode_literals

import sys
import os.path
import glob

if len(sys.argv[1:]) == 0:
    dirs = [os.getcwd()]
else:
    dirs = sys.argv[1:]

for directory in dirs:
    for notebook in glob.glob(os.path.join(directory, '*.ipynb')):
        cmd = 'ipython nbconvert --to rst {0}'.format(notebook)
        print(cmd)
        os.system(cmd)


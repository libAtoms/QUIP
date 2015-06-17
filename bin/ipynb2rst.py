#!/usr/bin/env python

import sys
import os
import glob

if len(sys.argv[1:]) == 0:
    dirs = [os.getcwd()]
else:
    dirs = sys.argv[1:]

for dir in dirs:
    for notebook in glob.glob(os.path.join(dir, '*.ipynb')):
        cmd = 'ipython nbconvert --to rst {0}'.format(notebook)
        print cmd
        os.system(cmd)


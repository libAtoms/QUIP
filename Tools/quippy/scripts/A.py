#!/usr/bin/env python

from qlab import *
import os
import sys

for filename in sys.argv[1:]:
    name = os.path.splitext(os.path.basename(filename))[0]
    name = name.replace('-','_').replace('.','_').replace('*','').replace('?','')
    print 'Creating AtomsViewer %s for file %s' % (name, filename)
    setattr(sys.modules[__name__], name, AtomsViewer(filename))

del name, filename

from IPython import embed
embed()


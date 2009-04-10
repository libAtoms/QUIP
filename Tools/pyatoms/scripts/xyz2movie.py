#!/usr/bin/env python

from pyatoms import *
import sys, os, itertools

viewfile = sys.argv[1]
xyzfile = sys.argv[2]
moviefile = sys.argv[3]

cmap = {0: "0.85 0.85 0.25\n",
        1: "1.00 0.00 0.00\n"}

script = open(viewfile).readlines()
script.append('shift_xtal 0 0\n') 

fr = itertools.islice(frame_reader(xyzfile),20)

for i, at in enumerate(fr):
   print 'Frame %d: %d atoms' % (i, at.n)

   colors = [cmap[e] for e in at.embed]

   clr_file = open('frame%05d.clr' % i,'w')
   clr_file.writelines(colors)
   clr_file.close()
   
   at.write_cfg('frame%05d.cfg' % i)
   
   script.append('load_config frame%05d.cfg\n' % i)
   script.append('load_atom_color frame%05d.clr\n' % i)
   script.append('capture_jpg frame%05d.jpg\n' % i)

script.append('quit\n')

script_file = open('atomeye.runscript','w')
script_file.writelines(script)
script_file.close()

os.system('A -nofep -nowindow -f=atomeye.runscript frame00000.cfg')

os.system('ffmpeg -i frame%05d.jpg -r 25 -b 10M '+moviefile)

# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HP X
# HP X   pyatoms: atomistic simulations tools
# HP X
# HP X   Copyright James Kermode 2010
# HP X
# HP X   These portions of the source code are released under the GNU General
# HP X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HP X
# HP X   If you would like to license the source code under different terms,
# HP X   please contact James Kermode, james.kermode@gmail.com
# HP X
# HP X   When using this software, please cite the following reference:
# HP X
# HP X   http://www.jrkermode.co.uk/PyAtoms
# HP X
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

import sys, re, os, os.path, string

from numpy import *
from atoms import *

# Cut out cracktip
def cracktip(atoms, r):
   centre = array([atoms.params['CrackPos'], 0.0, 0.0])
   return atoms.filter(norm(atoms.pos - centre) < r)


# Convert single xyz file to a bunch of .cfg files
def xyz2ecfg(infile, digits=4):
   inf = open(infile)
   basename = os.path.basename(sys.argv[1])

   if basename[-3:] == 'xyz':
      basename = basename[:-4]

   a = Atoms()
   i = 0
   while a.read(inf):
      thisfile = basename+(('%0'+str(digits)+'d') % i)+'.cfg'
      outf = open(thisfile, 'w')
      a.write_cfg(outf)
   
      sys.stdout.write('Atoms %d (%d atoms)     \r' % (i, a.n))
      sys.stdout.flush()
      i = i + 1
   
   print '\n'


# Find offsets in stream f where xyz frames start
def find_frames(f):
   frame_offsets = []

   # Regexp for start of frame matches line containing just n_atoms
   start = re.compile('^([\d]+)\n', re.MULTILINE)

   while 1:
      t = f.tell()
      line = f.readline()
      if not line:
         break

      # Is this line the start of a frame?
      m = start.search(line)

      if not m is None:
         frame_offsets.append(t)
        
         # Quickly skip whole frame without doing any regexp matching
         for j in range(int(m.group(1))+1):
            f.readline()

   return frame_offsets


def resample(infile, outfile, rate):
   inf = open(infile)
   outf = open(outfile, 'w')
   fs = find_frames(inf)
   a = Atoms()

   frame = 0
   for pos in fs[::rate]:
      inf.seek(pos)
      a.read(inf)
      frame = frame + 1
      a.write(outf)
      sys.stdout.write('Read frame %d (%d atoms)        \r'\
                       % (frame, a.n))
      sys.stdout.flush()
   print
   

# Read frames from xyz file
def read_frames(infile, start=None, end=None, step=None):
   inf = open(infile)
   frame_offsets = find_frames(inf)

   frames = range(len(frame_offsets))
   frames = frames[start:end:step]

   atoms_list = []

   for frame in frames:
      inf.seek(frame_offsets[frame])
      atoms_list.append(Atoms())
      if not atoms_list[-1].read(inf):
         raise IOError('error reading frame %d from %s' % (frame, infile))
                       
      sys.stdout.write('Read frame %d (%d atoms) - got %d frames        \r'\
                       % (frame, atoms_list[-1].n,len(atoms_list)))
      sys.stdout.flush()

   print
   inf.close()

   return atoms_list

# s is string list of atom indexes in fortran format (arrays start at 1)
def highlight(atoms, s):
   L = map(lambda x: int(x)-1, s.split())
   atoms.colour[:,:] = 0.0
   for a in L:
      if a < atoms.n:
         atoms.colour[a,:] = array([1,0,0])



## gnuplot = Gnuplot.Gnuplot()

## def plot_param(atomlist, x, ys):

##    data = {}

##    if x == 'Frame':
##       data[x] = range(len(atomlist))
##    else:
##       data[x] = [at.params[x] for at in atomlist]

##    for y, axis, with in ys:
##       data[y] = [at.params[y] for at in atomlist]

##    max_x = float(max(data[x]))
##    min_x = float(min(data[x]))
##    Nmax = len(data[x])

##    gnuplot('set xrange[%f:%f]' % (min_x, max_x))
##    gnuplot("set xlabel '%s'" % x)

##    # Pick first series with axes=x1y1 to label y axis and
##    # first with axes=x1y2 to label y2 axis
##    y_series = None
##    y2_series = None
##    for y, axes, with in ys:
##       if axes == 'x1y1' and not y_series:
##          y_series = y
##       if axes == 'x1y2' and not y2_series:
##          y2_series = y

##    gnuplot("set ylabel '%s'" % y_series)
##    if y2_series:
##       gnuplot("set y2label '%s'" % y2_series)
##       gnuplot("set y2tics border")
##       gnuplot("set ytics nomirror")

##    gnuplot('set terminal x11')

##    plotitems = []
##    for series, saxes, swith in ys:
##       plotitems.append(Gnuplot.Data(data[x][:Nmax], \
##                                     data[series][:Nmax], axes=saxes, \
##                                     with=swith, title=series))
##       gnuplot.plotcmd = 'plot'
##       gnuplot._clear_queue()
##       gnuplot._add_to_queue(plotitems)
##       gnuplot.refresh()


## def plot_param_list(paramslist, x, ys):

##    data = {}

##    if x == 'Frame':
##       data[x] = range(len(paramslist))
##    else:
##       data[x] = [params[x] for params in paramslist]

##    for y, axis, with in ys:
##       data[y] = [params[y] for params in paramslist]

##    max_x = float(max(data[x]))
##    min_x = float(min(data[x]))
##    Nmax = len(data[x])

##    gnuplot('set xrange[%f:%f]' % (min_x, max_x))
##    gnuplot("set xlabel '%s'" % x)

##    # Pick first series with axes=x1y1 to label y axis and
##    # first with axes=x1y2 to label y2 axis
##    y_series = None
##    y2_series = None
##    for y, axes, with in ys:
##       if axes == 'x1y1' and not y_series:
##          y_series = y
##       if axes == 'x1y2' and not y2_series:
##          y2_series = y

##    gnuplot("set ylabel '%s'" % y_series)
##    if y2_series:
##       gnuplot("set y2label '%s'" % y2_series)
##       gnuplot("set y2tics border")
##       gnuplot("set ytics nomirror")

##    gnuplot('set terminal x11')

##    plotitems = []
##    for series, saxes, swith in ys:
##       plotitems.append(Gnuplot.Data(data[x][:Nmax], \
##                                     data[series][:Nmax], axes=saxes, \
##                                     with=swith, title=series))
##       gnuplot.plotcmd = 'plot'
##       gnuplot._clear_queue()
##       gnuplot._add_to_queue(plotitems)
##       gnuplot.refresh()


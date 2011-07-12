#!/usr/bin/env python
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HQ X
# HQ X   quippy: Python interface to QUIP atomistic simulation library
# HQ X
# HQ X   Copyright James Kermode 2010
# HQ X
# HQ X   These portions of the source code are released under the GNU General
# HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HQ X
# HQ X   If you would like to license the source code under different terms,
# HQ X   please contact James Kermode, james.kermode@gmail.com
# HQ X
# HQ X   When using this software, please cite the following reference:
# HQ X
# HQ X   http://www.jrkermode.co.uk/quippy
# HQ X
# HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from quippy import *
from numpy import *
import sys, optparse, shutil, os

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from quippy.progbar import ProgressBar
import atomeye

def skip_bad_frames(seq):
   it = iter(seq)
   i = -1
   while True:
      i = i + 1
      print 'trying to load frame ', i
      try:
         at = it.next()
         yield at
      except RuntimeError:
         print 'skipping frame ', i
         continue

def skip_bad_times(seq, fill_value=9969209968386869046778552952102584320.0):
   """Skip frames with no Time param or with Time == fill_value."""
   
   for at in seq:
      if not 'Time' in at.params or abs(at.params['Time'] - fill_value) <= abs(numpy.finfo(double).eps*fill_value): continue 
      yield at


def skip_time_duplicates(seq, initial_time=-1.0):
   """Skip non-monotonically increasing frames in the Atoms sequence `seq`."""
   
   cur_time = initial_time
   for at in seq:
      if at.params['Time'] > cur_time:
         yield at
         cur_time = at.params['Time']

p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-r', '--range', action='store', help='Range of frames to include in output movie')
p.add_option('-o', '--outfile', action='store', help='Movie output file')
p.add_option('-n', '--nowindow', action='store_true', help='Disable AtomEye viewer window', default=False)
p.add_option('-E', '--encoder', action='store', help='Movie encoder command', default='ffmpeg -i %s -r 25 -b 30M %s')
p.add_option('-P', '--player', action='store', help='Movie player command', default='mplayer %s')
p.add_option('-l', '--load-view', action='store', help='Load view from AtomEye command script')
p.add_option('-s', '--save-view', action='store', help='Save view to AtomEye command script')
p.add_option('-g', '--graph', action='store_true', help='Enable graph plotting')
p.add_option('-x', '--xdata', action='store', help='Data for x axis of graph')
p.add_option('-y', '--ydata', action='store', help='Data for y axis of graph')
p.add_option('--y2data', action='store', help='Data for y2 axis of graph')
p.add_option('-t', '--text', action='store', help="""Annotate images with text. Argument should be a Python
expression which evaluates to a string, e.g. '"G = %f" % at.G'""")
p.add_option('-S', '--skip-bad-frames', action="store_true", help="Skip bad frames which cause an exception when loading")
p.add_option('-d', '--skip-time-duplicates', action="store_true", help="Skip duplicate frames (those with matching 'Time' params)")
p.add_option('-b', '--skip-bad-times', action='store_true', help="Skip frames with no 'Time' param or with where Time equals NetCDF fill value")
p.add_option('-m', '--merge', action='store', help="""Merge two input files. An auxilary input file name should be given.""")
p.add_option('-M', '--merge-properties', action='store', help="""List of properties to overwrite from MERGE file. Default is all properties.""")
p.add_option('-k', '--keep-images', action='store_true', help="""Do not remove frame images after movie has been made.""")
p.add_option('-p', '--property', action='store', help="""Property to use to colour atoms (default none)""")
p.add_option('-a', '--arrows', action='store', help="""Property to use to draw arrows (default none)""")
p.add_option('-e', '--exec_code', action='store', help="""Python code to execute on each frame before writing it to output file. Atoms object is
available as `at`.""")
p.add_option('-u', '--update', action='store_true', help="""Update a previous movie. Requires a previous run with -k, and implies -k on this run.""")
p.add_option('-W', '--width', action='store', help="""Width of output movie, in pixels.""", type='int')
p.add_option('-H', '--height', action='store', help="""Height of output movie, in pixels.""", type='int')
p.add_option('-A', '--aspect', action='store', help="""Aspect ratio. Used if only one of --width or --height is given. Default 0.75.""", default=0.75, type='float')
p.add_option('-R', '--rcut', action='append', help="""Following three arguments should be SYM1 SYM2 INCREMENT, to increment cutoff distance for SYM1-SYM2 bonds.""", nargs=3)

opt, args = p.parse_args()

if len(args) == 0:
   p.error('No input files specified')

if opt.range is not None:
   try:
      opt.range = parse_slice(opt.range)
   except:
      p.error('Cannot parse slice "%s" - should be in format [start]:[stop][:step]')
else:
   # Default is all frames
   opt.range = slice(0, None, None)

if opt.outfile is None:
   opt.outfile = os.path.splitext(args[0])[0] + '.mp4'

sources = [AtomsReader(f) for f in args]

# Try to count number of frames
nframes = sum([len(s) for s in sources])
nframes = len(range(*opt.range.indices(nframes)))

ndigit = 5
basename, ext = os.path.splitext(opt.outfile)
out_fmt = '%s%%0%dd.jpg' % (basename, ndigit)

# Look for existing JPEGs, so we can continue a previous run
if opt.update:
   opt.keep_images = True
   
   existing_frames = [ os.path.exists(out_fmt % i) for i in range(nframes) ]
   restart_frame = existing_frames.index(False)
   
   print 'Restarting from frame %d' % restart_frame
   if opt.range.step is None:
      opt.range = slice(opt.range.start + restart_frame, opt.range.stop, opt.range.step)
   else:
      opt.range = slice(opt.range.start + restart_frame/opt.range.step, opt.range.stop, opt.range.step)

   nframes = sum([len(s) for s in sources])
   nframes = len(range(*opt.range.indices(nframes)))
   img_offset = restart_frame
else:
   img_offset = 0


# Build a chain of iterators over all input files, skipping frames as appropriate
atomseq = itertools.islice(itertools.chain.from_iterable(sources),
                           opt.range.start, opt.range.stop, opt.range.step)

if opt.skip_bad_frames:
   atomseq = skip_bad_frames(atomseq)

if opt.skip_bad_times:
   atomseq = skip_bad_times(atomseq)

if opt.skip_time_duplicates:
   atomseq = skip_time_duplicates(atomseq)

if opt.merge is not None:
   merge_config = Atoms(opt.merge)

   if opt.merge_properties is not None:
      opt.merge_properties = parse_comma_colon_list(opt.merge_properties)
   else:
      opt.merge_properties = merge_config.properties.keys()

a0 = Atoms(args[0], frame=opt.range.start)

if opt.merge is not None:
   for k in opt.merge_properties:
      a0.add_property(k, getattr(merge_config, k))

view = atomeye.AtomEyeView(nowindow=opt.nowindow)
view.show(a0, property=opt.property, arrows=opt.arrows)

if opt.load_view is not None:
   view.run_script(opt.load_view)
   view.redraw()

if opt.rcut is not None:
   for (sym1, sym2, rcut) in opt.rcut:
      view.rcut_patch(sym1, sym2, float(rcut))

if opt.width is not None or opt.height is not None:
   print 'width', opt.width, 'height', opt.height, 'aspect', opt.aspect
   if opt.width  is None: opt.width = int(opt.height/opt.aspect)
   if opt.height is None: opt.height = int(opt.width*opt.aspect)
   view.resize(opt.width, opt.height)

if not opt.nowindow:
   raw_input('Arrange AtomEye view then press enter...')

if opt.save_view is not None:
   view.save(opt.save_view)

postprocess = None

if opt.graph is not None:

   def add_plot(at, i, filename):

      l1.set_data(xdata[:i], ydata[:i])
      if opt.y2data is not None:
         l2.set_data(xdata[:i], y2data[:i])

      basename, ext = os.path.splitext(filename)

      try:
         fig.savefig(basename+'.1.png')
         os.system('montage -geometry +0+0 -tile 1x2 %s %s %s' % (filename, basename+'.1.png', basename+'.2.jpg'))
         shutil.move(basename+'.2.jpg', filename)
      finally:
         os.remove(basename+'.1.png')
   
   postprocess = add_plot
   if opt.xdata == 'frame':
      xdata = farray(range(nframes))
   else:
      xdata = hstack([getattr(s,opt.xdata) for s in sources])[opt.range]
   ydata  = hstack([getattr(s,opt.ydata) for s in sources])[opt.range]
   if opt.y2data is not None:
      y2data = hstack([getattr(s,opt.y2data) for s in sources])[opt.range]

   fig = Figure(figsize=(8,2))
   canvas = FigureCanvas(fig)

   ax1 = fig.add_axes([0.1,0.2,0.8,0.65])
   ax1.set_xlim(xdata.min(), xdata.max())
   ax1.set_ylim(ydata.min(), ydata.max())
   l1, = ax1.plot([], [], 'b-', scalex=False, scaley=False, label='Temp')
   ax1.set_xlabel(opt.xdata)
   ax1.set_ylabel(opt.ydata)

   if opt.y2data is not None:
      ax2 = ax1.twinx()
      ax2.set_ylim(y2data.min(), y2data.max())
      ax2.yaxis.tick_right()
      l2, = ax2.plot([], [], 'r-', scalex=False, scaley=False, label='G')
      ax2.set_ylabel(opt.y2data)

if opt.text is not None:

   def add_text(at, i, filename):
      s = eval(opt.text)
      os.system('mogrify -gravity SouthWest -annotate 0x0+0+0 "%s" -font Helvetica -pointsize 32 %s' % (s, filename))

   postprocess = add_text

progress = nframes is not None
imgs = []
if progress: pb = ProgressBar(0,nframes,80,showValue=True)
try:
   if progress: print 'Rendering frames...'
   for i, at in enumerate(atomseq):
      filename = out_fmt % (i + img_offset)

      if opt.merge:
         for k in opt.merge_properties:
            at.add_property(k, getattr(merge_config, k))

      if opt.exec_code is not None:
         exec(opt.exec_code)
     
      view.show(at, property=opt.property, arrows=opt.arrows)
      view.capture(filename)
      view.wait()
      if postprocess is not None:
          postprocess(at, i, filename)
      imgs.append(filename)
      if progress: pb(i+1)
   if progress: print

   if opt.encoder is not None:
      if progress: print 'Encoding movie'
      os.system(opt.encoder % (out_fmt, opt.outfile))

   if opt.player is not None:
      if progress: print 'Playing movie'
      os.system(opt.player % opt.outfile)

finally:
   if not opt.keep_images:
      view.wait()
      for img in imgs:
         if os.path.exists(img): os.remove(img)
         if os.path.exists(img+'.cmap.eps'): os.remove(img+'.cmap.eps')



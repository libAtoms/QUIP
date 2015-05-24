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
import sys
import optparse
import shutil
import os
import itertools
import select
import time
import atomeye
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from quippy.progbar import ProgressBar


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

def open_sources(args, frames, indices):
   if len(args) == 1:
      atomseq = AtomsReader(args[0], start=frames.start,
                            stop=frames.stop, step=frames.step, indices=indices)
      nframes = len(atomseq)
      sources = [atomseq]
   else:
      sources = [AtomsReader(f, indices=indices) for f in args]

      # Try to count total number of frames
      nframes = sum([len(s) for s in sources])
      nframes = len(range(*frames.indices(nframes)))

      # Build a chain of iterators over all input files, skipping frames as appropriate
      atomseq = itertools.islice(itertools.chain.from_iterable(sources),
                                 frames.start, frames.stop, frames.step)

   return nframes, atomseq, sources


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


def add_text(at, i, filename):
   s = eval(opt.text)
   os.system('mogrify -gravity SouthWest -annotate 0x0+0+0 "%s" -font Helvetica -pointsize 32 %s' % (s, filename))


def block_until_more_frames(filenames, current_length, poll_interval=30.0):
   idx_files = []
   for filename in filenames:
      if not filename.endswith('.xyz'):
         raise IOError('block_until_more_frames() only works with .xyz files')

      idx_file = filename+'.idx'
      if not os.path.exists(idx_file):
         raise IOError('Index file %s not found' % idx_file)
      idx_files.append(idx_file)

   while True:
      length = 0
      for idx_file in idx_files:
         try:
            f = open(idx_file, 'r')
            length += int(f.readline())
         except IOError:
            time.sleep(poll_interval)
            continue
         finally:
            f.close()
      if length > current_length:
         return
      time.sleep(poll_interval)

   

p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-r', '--range', action='store', help='Range of frames to include in output movie')
p.add_option('-o', '--outfile', action='store', help='Movie output file')
p.add_option('-n', '--nowindow', action='store_true', help='Disable AtomEye viewer window', default=False)
p.add_option('-N', '--nomovie', action='store_true', help='Do not make a movie, just render frames', default=False)
p.add_option('-E', '--encoder', action='store', help='Movie encoder command', default='ffmpeg -i %s -r 25 -b 30M %s')
p.add_option('-P', '--player', action='store', help='Movie player command', default='open %s')
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
p.add_option('-w', '--watch', action='store_true', help="""When completed, keep watching for new frames to be written to input stream. Implies -u, -k and -N.""")
p.add_option('-L', '--latest', action='store', help="""Copy latest snapshot image to this file""")
p.add_option('-F', '--fix-indices', action='store_true', help="""Fix the set of atoms loaded from each frame to those visible in viewer window in first frame""")

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

orig_range = opt.range

if opt.outfile is None:
   opt.outfile = os.path.splitext(args[0])[0] + '.mp4'

if opt.watch:
   if sys.platform != 'darwin':
      p.error('-w/--watch mode only implemented on Mac OS X so far')
   opt.update = True
   opt.keep_images = True
   opt.nomovie = True

if opt.nomovie:
   opt.encoder = None
   opt.player = None

view = None
indices = None

# outer loop for opt.watch
while True:

   nframes, atomseq, sources = open_sources(args, orig_range, indices)
   print 'Found %d frames' % nframes
   orig_nframes = nframes
   ndigit = 5
   basename, ext = os.path.splitext(opt.outfile)
   out_fmt = '%s%%0%dd.jpg' % (basename, ndigit)

   if view is None:
      if opt.range.start != 0:
         a0 = Atoms(args[0], frame=opt.range.start)
      else:
         a0 = AtomsReader(args[0])[0]

      if opt.merge is not None:
         for k in opt.merge_properties:
            a0.add_property(k, getattr(merge_config, k))

      view = atomeye.AtomEyeViewer(nowindow=opt.nowindow, verbose=False)
      view.show(a0, property=opt.property, arrows=opt.arrows)
      if opt.load_view is not None:
         view.run_script(opt.load_view)
         view.show()

      rcut_patches = []
      if opt.rcut is not None:
         rcut_patches.extend(opt.rcut)
      for (sym1, sym2, rcut) in rcut_patches:
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

      view.wait()

      if opt.fix_indices:
         indices = view.get_visible()
         print 'Restricing atom indices to the %d visible atoms' % len(indices)
         nframes, atomseq, sources = open_sources(args, orig_range, indices)

   # Look for existing JPEGs, so we can continue a previous run
   if opt.update:
      opt.keep_images = True

      existing_frames = [ os.path.exists(out_fmt % i) for i in range(nframes) ]
      if False in existing_frames:
         restart_frame = existing_frames.index(False)
      else:
         p.error("no new frames after %d for -u/--update mode to process" % nframes)

      print 'Restarting from frame %d' % restart_frame
      if opt.range.step is None:
         opt.range = slice(orig_range.start + restart_frame, orig_range.stop, orig_range.step)
      else:
         opt.range = slice(orig_range.start + restart_frame/orig_range.step, orig_range.stop, orig_range.step)

      # Re-open the sources starting at restart_frame
      nframes, atomseq, sources = open_sources(args, opt.range, indices)
      img_offset = restart_frame
   else:
      img_offset = 0

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

   postprocess = None

   if opt.graph is not None:
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

      if opt.latest is not None:
         shutil.copyfile(imgs[-1], opt.latest)

   finally:
      if not opt.keep_images:
         view.wait()
         for img in imgs:
            if os.path.exists(img): os.remove(img)
            if os.path.exists(img+'.cmap.eps'): os.remove(img+'.cmap.eps')

   if opt.watch:
      print "Waiting for new frames..."
      block_until_more_frames(args, orig_nframes)
   else:
      break
   




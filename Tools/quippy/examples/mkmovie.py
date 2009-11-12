from quippy import *
import sys, optparse, shutil

from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-r', '--range', action='store', help='Range of frames to include in output movie')
p.add_option('-o', '--outfile', action='store', help='Movie output file')
p.add_option('-n', '--nowindow', action='store_true', help='Disable AtomEye viewer window', default=False)
p.add_option('-e', '--encoder', action='store', help='Movie encoder command', default='ffmpeg -i %s -r 25 -b 30M %s')
p.add_option('-p', '--player', action='store', help='Movie player command', default='mplayer %s')
p.add_option('-v', '--viewfile', action='store', help='AtomEye command script for initial view')

opt, args = p.parse_args()

if len(args) == 0:
   p.error('No input files specified')

if opt.range is not None:
   opt.range = slice(*[{True: lambda n: None, False: int}[x == ''](x) for x in (opt.range.split(':') + ['', '', ''])[:3]])
else:
   opt.range = slice(1, None, None)

if opt.outfile is None:
   opt.outfile = os.path.splitext(args[0])[0] + '.mp4'

sources = [AtomsList(f, store=False) for f in args]

# Try to count number of frames
nframes = sum([len(s) for s in sources])

nframes = len(range(*opt.range.indices(nframes)))

# Build a chain of iterators over all input files, skipping frames as appropriate
atomseq = itertools.islice(itertools.chain.from_iterable(sources),
                           opt.range.start+1, opt.range.stop, opt.range.step)

a0 = Atoms(args[0])
view = atomeye.show(a0, nowindow=opt.nowindow)

if opt.viewfile is not None:
   view.load_script(opt.viewfile)
   view.shift_xtal(0, 0)
   view.redraw()
   view.wait()

if not opt.nowindow:
   raw_input('Arrange AtomEye view then press enter...')

fillvalue = 9.969209968386869e+36

xdata = hstack([s.Time for s in sources])[opt.range]
ydata  = hstack([s.Temp for s in sources])[opt.range]
y2data = hstack([s.G for s in sources])[opt.range]

fig = Figure(figsize=(8,2))
canvas = FigureCanvas(fig)

ax1 = fig.add_axes([0.1,0.2,0.8,0.65])
ax1.set_xlim(xdata.min(), xdata.max())
ax1.set_ylim(ydata.min(), ydata.max())
l1, = ax1.plot([], [], 'b-', scalex=False, scaley=False, label='Temp')
ax1.set_xlabel('Time / fs ')
ax1.set_ylabel('Temp / K')

ax2 = ax1.twinx()
ax2.set_ylim(y2data.min(), y2data.max())
ax2.yaxis.tick_right()
l2, = ax2.plot([], [], 'r-', scalex=False, scaley=False, label='G')
ax2.set_ylabel('G / Jm$^{-2}$')


def add_plot(at, i, filename):

   l1.set_data(xdata[:i], ydata[:i])
   l2.set_data(xdata[:i], y2data[:i])

   basename, ext = os.path.splitext(filename)

   try:
      fig.savefig(basename+'.1.png')
      os.system('montage -geometry +0+0 -tile 1x2 %s %s %s' % (filename, basename+'.1.png', basename+'.2.jpg'))
      shutil.move(basename+'.2.jpg', filename)
   finally:
      os.remove(basename+'.1.png')

   

def add_text(at, i, filename):
   s = 'G = %.3f' % at.G
   os.system('mogrify -gravity SouthWest -annotate 0x0+0+0 "%s" -font Helvetica -pointsize 32 %s' % (s, filename))

view.make_movie(atomseq, opt.outfile, postprocess=add_plot,
                movieencoder=opt.encoder, movieplayer=opt.player, nframes=nframes, cleanup=True)


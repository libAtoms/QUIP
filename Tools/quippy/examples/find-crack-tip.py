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

from pylab import *
from quippy import *

# Use percolation algorithm to find crack tip

## def percolation_step(v, nstep):
##    nx, ny, nz = v.shape

##    # copy of v with zeros around edges to avoid IndexError
##    vb = zeros((nx+2,ny+2,nz+2),int)
##    vb[1:nx+1,1:ny+1,1:nz+1] = v.copy()

##    v[v >= 2] += 1

##    for i in frange(nx):
##       for j in frange(ny):
##          for k in frange(nz):

##             # Only consider cells which could be infected
##             if vb[i,j,k] != 1: continue

##             try:
##                if (vb[i+1,j,k] == 2 or vb[i-1,j,k] == 2 or
##                    vb[i,j+1,k] == 2 or vb[i,j-1,k] == 2 or
##                    vb[i,j,k+1] == 2 or vb[i,j,k-1] == 2):
##                   v[i,j,k] = 2
##                   raise StopIteration()

##             except StopIteration:
##                pass

##    print 'alive', (v == 2).sum()
##    return (v == 2).sum() != 0

cmap_data = {'red':   ((0, 0., 0.),(0.365079, 1.000000, 1.000000),(1.0, 1.0, 1.0)),
             'green': ((0, 0., 0.),(0.365079, 0.000000, 0.000000),(0.746032, 1.000000, 1.000000),(1.0, 1.0, 1.0)),
             'blue':  ((0, 0., 0.),(0.746032, 0.000000, 0.000000),(1.0, 1.0, 1.0))}

my_cmap = get_cmap('jet_r')
my_cmap.set_under('k', 1.0)

def percolate(v, fig=None, base=None):

   if fig is not None:
      fig.set_size_inches(6, 6)
      fig.clf()
      ax = fig.add_axes([0,0,1,1],frameon=False)

   nsteps = 0
   while percolation_step(v):
      nsteps = nsteps +1

      if fig is not None:
         z = v.shape[2]/2+1
         vis = flipud(v[z].T)
         if base is not None:
            ax.cla()
            ax.imshow(base, extent=(0,vis.shape[1], 0, vis.shape[0]))
         ax.pcolor(vis, cmap=my_cmap, vmin=1, edgecolor='k', alpha=0.4)
         #ax.set_xlim(12, 75)
         #ax.set_ylim(32, 52)
         ax.set_xticks([])
         ax.set_yticks([])
         draw()
         raw_input()
         #fig.savefig('img%05d.png' % nsteps)

   return nsteps
   
def forest(nx,ny,nz,p):
   v = fzeros((nx,ny,nz))
   for i in frange(nx):
      for j in frange(ny):
         for k in frange(nz):
            if random.random() > 1-p: v[i,j,k] = 1
   # seed first row
   v[1,v[1,:,1] == 1,1] = 2
   return v




def percolate_crack(at, tol, fig, base, start_pos=[0.0, 0.0, 0.0]):

   at = at.copy()
   at.set_cutoff(tol)
   at.calc_connect()

   cells = fzeros((at.connect.cellsna, at.connect.cellsnb, at.connect.cellsnc),int)

   for i in frange(at.connect.cellsna):
      for j in frange(at.connect.cellsnb):
         for k in frange(at.connect.cellsnc):
            cells[i,j,k] = at.connect.cell_n(i,j,k) == 0


   i, j, k = at.connect.cell_of_pos(at.g, start_pos)

   # stop from going backwards
   cells[1:i,:,:] = 0

   # seed percolation at left edge
   cells[i,j,k] = 2

   percolate(cells, fig=fig, base=base)

   # Last cell visited in percolation walk
   crack_cell = (cells == 3).nonzero()

   crack_t = (farray([float(x) for x in crack_cell])/
              farray([at.connect.cellsna, at.connect.cellsnb, at.connect.cellsnc])) - [0.5, 0.5, 0.5]

   crack_pos = dot(at.lattice, crack_t)

   return crack_pos, cells


def img(at):
   view = atomeye.show(at)
   view.capture('tmp.png')
   view.wait()
   im = imread('tmp.png')

   ima = zeros((im.shape[0], im.shape[1], 4))
   ima[:,:,0:3] = im

   # make white transparent
   ima[:,:,3] = 0
   ima[:,:,3] = 1-(im == 1).all(axis=2)

   return ima


               
def is_through_crack(cp, cells):
   # see if there is a connected route through entire system
   vert_slice = cells[max(1,cells.shape[0]/2),:,max(1,cells.shape[2]/2)]
   occupied_cells, = (vert_slice != 1).nonzero()
   top_edge, bottom_edge = occupied_cells[2], occupied_cells[-2]

   print top_edge, bottom_edge

   if cp.crack_double_ended:
      left_edge = 1
      right_edge = -1
   else:
      horz_slice = cells[:,max(1,cells.shape[1]/2),max(1,cells.shape[2]/2)]
      occupied_cells, = (horz_slice != 1).nonzero()
      left_edge, right_edge = occupied_cells[2], occupied_cells[-2]
      
   print left_edge, right_edge
   return any(cells[left_edge,top_edge:bottom_edge,:] > 1) and any(cells[right_edge,top_edge:bottom_edge,:] > 1)


def apply_min_filter(cells, start_pos, min_dist):
   from scipy.ndimage import minimum_filter

   ncells = cells.copy()
##    ncells[ncells == 0] = 2*ncells.max()
   ncells = farray(minimum_filter(ncells, size=min_dist))

   minima = [tuple(x) for x in transpose((logical_and(ncells == cells, ncells != ncells.max())).nonzero())]

   class Duplicate: pass

   while True:
      try:
         for min1,min2 in combinations(minima,2):
            if (farray(min1) - farray(min2)).norm() < min_dist:
               print 'duplicate', min1, min2
               raise Duplicate
         else:
            break
      except Duplicate:
         d1 = (farray(min1) - farray(start_pos)).norm()
         d2 = (farray(min2) - farray(start_pos)).norm()
         print d1, d2
         if d1 < d2:
            print 'removing', min2
            minima.remove(min2)
         else:
            print 'removing', min1
            minima.remove(min1)

   return minima


      
verbosity_push(VERBOSE)

a = Atoms('quartz_movie_2.nc',frame=0)
cp = CrackParams()
io = InOutput('quartz.xml')
cp.read_xml(io)
io.close()

crack_tips = crack_find_tips(a, cp)

print crack_tips

for pos in crack_tips.real:
   a.add_atoms(pos, 29)

a.show()

verbosity_pop()

def count_cells(at,size):
   at.set_cutoff(size)
   at.calc_connect()
   cells = fzeros((at.connect.cellsna, at.connect.cellsnb, at.connect.cellsnc),int)

   for i in frange(at.connect.cellsna):
      for j in frange(at.connect.cellsnb):
         for k in frange(at.connect.cellsnc):
            cells[i,j,k] = at.connect.cell_n(i, j, k)

   return cells


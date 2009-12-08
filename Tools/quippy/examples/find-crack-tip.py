from pylab import *
from quippy import *

# Use percolation algorithm to find crack tip

def percolation_step(v, nstep):
   nx, ny, nz = v.shape

   # copy of v with zeros around edges to avoid IndexError
   vb = zeros((nx+2,ny+2,nz+2),int)
   vb[1:nx+1,1:ny+1,1:nz+1] = v.copy()

   v[v >= 2] += 1

   for i in frange(nx):
      for j in frange(ny):
         for k in frange(nz):

            # Only consider cells which could be infected
            if vb[i,j,k] != 1: continue

            try:
               if (vb[i+1,j,k] == 2 or vb[i-1,j,k] == 2 or
                   vb[i,j+1,k] == 2 or vb[i,j-1,k] == 2 or
                   vb[i,j,k+1] == 2 or vb[i,j,k-1] == 2):
                  v[i,j,k] = 2
                  raise StopIteration()

            except StopIteration:
               pass

   print 'alive', (v == 2).count()
   return (v == 2).count() != 0

def percolate(v):
   nsteps = 0
   imshow(v[1], interpolation='nearest')
   draw()
   while percolation_step(v, nsteps):
      nsteps = nsteps +1
      imshow(v[1], interpolation='nearest')
      draw()
   imshow(v[1], interpolation='nearest')
   draw()
   return nsteps
   
def forest(nx,ny,nz):
   v = farray(numpy.random.randint(2,size=nx*ny*nz).reshape((nx,ny,nz),order='F'))
   # seed first row
   v[1,v[1,:,1] == 1,1] = 2
   return v


def percolate_crack(at, tol):

   at = at.copy()
   at.set_cutoff(tol)
   at.calc_connect()

   cells = fzeros((at.connect.cellsna, at.connect.cellsnb, at.connect.cellsnc),int)

   for i in frange(at.connect.cellsna):
      for j in frange(at.connect.cellsnb):
         for k in frange(at.connect.cellsnc):
            cells[i,j,k] = at.connect.cell_n(i,j,k) == 0


   start_pos = [-at.params['OrigWidth']/2.0, 0.0, 0.0]
   i, j, k = at.connect.cell_of_pos(at.g, start_pos)

   # stop from going backwards
   cells[1:i,:,:] = 0

   # seed percolation at left edge
   cells[i,j,k] = 2

   percolate(cells)

   # Last cell visited in percolation walk
   crack_cell = (cells == 3).nonzero()

   crack_t = (farray([float(x) for x in crack_cell])/
              farray([at.connect.cellsna, at.connect.cellsnb, at.connect.cellsnc])) - [0.5, 0.5, 0.5]

   crack_pos = dot(at.lattice, crack_t)

   return crack_pos

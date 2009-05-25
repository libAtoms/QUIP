from quippy import *
from math import ceil
import time

def alpha_quartz(a=4.9134,c=5.4052, x1=0.4699, x2=0.4141, y2=0.2681, z2=0.7854):
   """Primitive 9-atom orthorhombic alpha quartz cell"""

   a1 = farray((0.5*a, -0.5*sqrt(3.0)*a, 0.0))
   a2 = farray((0.5*a,  0.5*sqrt(3.0)*a, 0.0))
   a3 = farray((0.0,    0.0,             c))

   lattice = fzeros((3,3))
   lattice[:,1] = a1
   lattice[:,2] = a2
   lattice[:,3] = a3
   
   at = Atoms(9,lattice)

   at.set_atoms((14,14,14,8,8,8,8,8,8))

   at.pos[:,1] =  x1*a1 + 2.0/3.0*a3
   at.pos[:,2] =  x1*a2 + 1.0/3.0*a3
   at.pos[:,3] = -x1*a1 - x1*a2
   at.pos[:,4] =  x2*a1 + y2*a2 + z2*a3
   at.pos[:,5] = -y2*a1 + (x2-y2)*a2  + (2.0/3.0 + z2)*a3
   at.pos[:,6] = (y2-x2)*a1 - x2*a2   + (1.0/3.0 + z2)*a3
   at.pos[:,7] = y2*a1 + x2*a2 - z2*a3
   at.pos[:,8] = -x2*a1 + (y2-x2)*a2 + (2.0/3.0 - z2)*a3
   at.pos[:,9] = (x2 - y2)*a1 - y2*a2 + (1.0/3.0 - z2)*a3
   
   return at

def alpha_quartz_cubic(*args, **kwargs):
   """Non-primitive 18-atom cubic quartz cell."""
   a0 = alpha_quartz(*args, **kwargs)
   at = supercell(a0, 4, 4, 1)
   at.map_into_cell()

   lattice = fzeros((3,3))
   lattice[1,1] = a0.lattice[1,1]*2.0
   lattice[2,2] = a0.lattice[2,2]*2.0
   lattice[3,3] = a0.lattice[3,3]

   g = linalg.inv(lattice)
   t = dot(g, at.pos)
   cubic = at.select(logical_and(t >= -0.5, t < 0.5).all(axis=1))
   cubic.set_lattice(lattice)
   return cubic

def beta_quartz():
   pass

def bracket(func, x1, x2, max_iter=50, factor=1.6, **kwargs):
   f1 = func(x1, **kwargs)
   f2 = func(x2, **kwargs)
   for j in range(max_iter):
      if f1*f2 < 0.0:
         return (x1, x2)
      if abs(f1) < abs(f2):
         x1 += factor*(x1 - x2)
         f1 = func(x1, **kwargs)
      else:
         x2 += factor*(x2 - x1)
         f2 = func(x2, **kwargs)
   raise ValueError('Maximum number of iterations exceeded.')

def bisect(func, x1, x2, err=1e-5, max_iter=50, **kwargs):
   f = func(x1, **kwargs)
   fmid = func(x2, **kwargs)
   if f*fmid >= 0.0:
      raise ValueError("Root not bracketed")

   if f < 0.0:
      dx = x2 - x1
      rtb = x1
   else:
      dx = x1 - x2
      rtb = x2

   for j in range(max_iter):
      xmid = rtb + dx
      dx *= 0.5
      fmid = func(xmid, **kwargs)
      if fmid < 0:
         rtb = xmid
      if abs(dx) < err or fmid == 0.0:
         return rtb
   raise ValueError('Maximum number of iterations exceeded.')   

def newton_raphson(func, dfunc, x1, x2, err=1e-5, max_iter=20, **kwargs):
   x = 0.5*(x1 + x2)
   for j in range(max_iter):
      f = func(x, **kwargs)
      df = dfunc(x, **kwargs)
      dx = f/df
      x = x - dx
      print j, x
      if (x1 - x)*(x - x2) < 0.0:
         raise ValueError('Jumped out of brackets')
      if abs(dx) < err:
         return x
   raise ValueError('Maximum number of iterations exceeded.')

def rcut_func(r, alpha, eps):
   return exp(-alpha*r)/r - eps

def rcut_dfunc(r, alpha, eps):
   return -alpha*exp(-alpha*r)/r - exp(-alpha*r)/r**2

def rcut(alpha, eps=1.0/80.0, min=1.0, max=100.0):
   min, max = bracket(rcut_func, min, max, alpha=alpha, eps=eps)
   return bisect(rcut_func, min, max, alpha=alpha, eps=eps)


def force_test(at, p, dx=1e-4):
   analytic_f = fzeros((3,at.n))
   p.calc(at, f=analytic_f)
   num_f = fzeros((3,at.n))
   ep, em = farray(0.0), farray(0.0)

   for i in frange(at.n):
      for j in (1,2,3):
         ap = at.copy()
         ap.pos[j,i] += dx
         p.calc(ap, e=ep)
         print 'e+', j,i,ep
         ap.pos[j,i] -= 2.0*dx
         p.calc(ap, e=em)
         print 'e-', j,i,em
         num_f[j,i] = -(ep - em)/(2*dx)

   return analytic_f, num_f, analytic_f - num_f


times = {}

alpha = 0.0

xml = """<ASAP_params cutoff="4.0" n_types="2">
<per_type_data type="1" atomic_num="8" />
<per_type_data type="2" atomic_num="14" />
<params>
15.9994 28.086
O Si
48 24
.f. 0.0 1.d-9  0.0 1 1 1 %f 21.0 18.0 0.0 0.0 raggio,A_ew,gcut,iesr,rcut
 -1.38257 2.76514
 -------------Alphaij---------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 -------------Bij--------------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 ---------------Cij------------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 -----------------Dij----------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 ---------------Eij------------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 ---------------Nij------------------
 0.0000000E+00   8.0000000E+00
 0.0000000E+00
 ---------Tang-Toennies-------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 ---------Tang-Toennies-------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 ---------Tang-Toennies-------------
 0.0000000E+00   0.0000000E+00
 0.0000000E+00
 ---------------D_ms----------------
 2.4748d-4   1.9033d-3        
-2.0846d-3    
 ---------------Gamma_ms------------
 12.07092 11.15230
 10.45517
 ----------------R_ms---------------
 7.17005 4.63710
 5.75038
 --------------Polarization---------
 8.89378       0.0d0
 0.70,  60,  1.0000000001E-7 2
 ---------------Bpol----------------
 0.0000000E+00   2.02989      
 0.00000000000
 ---------------Cpol----------------
 0.0000000E+00  -1.50435       
 0.00000000    
 --------Aspherical-Ion-Model-------
 F,  F,  7.,  8.
 -----------Adist-------------------
 0.0000000E+00   2.0170894E+00
 2.4232942E+00
 ---------------Bdist---------------
 0.0000000E+00   7.6306646E+01
 1.5861246E+03
 ---------------Cdist---------------
 0.0000000E+00  -1.2069760E-02
 --------------Ddist----------------
 0.0000000E+00  -4.0995369E-02
 2.2483367E-02
 -------------Sigma1----------------
 0.0000000E+00  -1.4410513E+07
 -------------Sigma2----------------
 0.0000000E+00  -5.1477595E+00
 ------------Sigma3-----------------
 0.0000000E+00   1.1143606E+08
 -------------Sigma4----------------
 0.0000000E+00   7.2089861E+00
 -------------Bu1-------------------
 4.1063828E+12   1.8240403E+02
-2.7852429E+04
 --------------Alphau1--------------
 7.2970202E+00   2.2221123E+00
 2.9876383E+00
 ---------------Bu2-----------------
-3.4880044E+13  -2.0079894E+03
 4.0014253E+03
 --------------Alphau2--------------
 7.8085212E+00   3.7185181E+00
 2.4488279E+00
***********Spring Constant*********
 1.0 1.0 
 1.0
***********Spring Cutoff***********
 3.60 2.30
   3.8
**********Smoothing***************
0.0 0.0 .f.
*************Yukawa***************
%f 10.0 .t.
</params>
</ASAP_params>""" 


reps = (1, 2, 3, 4)
alphas = (0.0, 0.01, 0.05, 0.1)

times = {}
fvar('e')
at = alpha_quartz_cubic()
for rep in reps:
   aa = supercell(at, rep, rep, rep)
   for alpha in alphas:
      p = Potential('IP ASAP', xml % (rcut(alpha), alpha))

      t1 = time.time()
      p.calc(aa, e=e)
      t2 = time.time()

      times[(rep,alpha)] = (aa.n, t2-t1)


#for alpha in alphas:
#   x = []
#   y = []
#  for rep in reps:
#      x.append(at.n*rep**3)
#      y.append(times[(rep,alpha)])
#   plot(x,y)
#legend(alphas, 2)

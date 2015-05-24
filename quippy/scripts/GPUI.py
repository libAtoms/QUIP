#!/usr/bin/env python

from numpy import *
# from scipy.special import erf
import sys

if (len(sys.argv) != 10):
   sys.stderr.write(("Usage: %s collective_val_column dE_dcollective_column noise_variance_column\n"% sys.argv[0])+
                    "       len_scale variance_prior min_collective max_collective n_collective min_bracket\n" )
   sys.exit(1)

def p_erf(xi):
   # save the sign of x
   # sign = 1 if x >= 0 else -1
   if (rank(xi) > 0):
      xr=zeros((len(xi)))
      for i in range(len(xi)):
	 x = xi[i]
	 if x >= 0:
	    sign = 1
	 else:
	    sign = -1
	 x = abs(x)

	 # constants
	 a1 =  0.254829592
	 a2 = -0.284496736
	 a3 =  1.421413741
	 a4 = -1.453152027
	 a5 =  1.061405429
	 p  =  0.3275911

	 # A&S formula 7.1.26
	 t = 1.0/(1.0 + p*x)
	 y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
	 xr[i] = sign*y
      # return sign*y # erf(-x) = -erf(x)
      return xr
   else:
      x = xi
      if x >= 0:
	 sign = 1
      else:
	 sign = -1
      x = abs(x)

      # constants
      a1 =  0.254829592
      a2 = -0.284496736
      a3 =  1.421413741
      a4 = -1.453152027
      a5 =  1.061405429
      p  =  0.3275911

      # A&S formula 7.1.26
      t = 1.0/(1.0 + p*x)
      y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*math.exp(-x*x)
      return sign*y


# read command line (should use getopt?)
collective_val_col=int(sys.argv[1])
dE_dcollective_col=int(sys.argv[2])
noise_variance_col=int(sys.argv[3])
len_scale=float(sys.argv[4]) # length scale
variance_prior=float(sys.argv[5]) #prior function variance; large values give less informative priors
min_collective = float(sys.argv[6])
max_collective = float(sys.argv[7])
n_collective = int(sys.argv[8])
min_bracket = float(sys.argv[9])

# read inputs into arrays
collective_list = []
dE_dcollective_list = []
noise_variance_list = []
decimation=1
for i, line in enumerate(sys.stdin):
   if (line[0] == '#'):
      continue
   if ((i % decimation) == 0):
      sys.stderr.write(".")
      f = line.split()
      collective_list.append(float(f[collective_val_col]))
      dE_dcollective_list.append(float(f[dE_dcollective_col]))
      noise_variance_list.append(float(f[noise_variance_col]))
      if ((i % (80*decimation)) == 0):
	 sys.stderr.write("\n")
sys.stderr.write("\n")

collective_vec = array(collective_list)
dE_dcollective_vec = array(dE_dcollective_list)
noise_variance_vec = array(noise_variance_list)


# map 1.0 (minimum) to 0, and transition to new slope around 1.25??
def coord_transform(x):
     #NB return (erf((x-1.0)*4.0)+0.5*x)
     return (p_erf((x-1.0)*4.0)+0.5*x)
def coord_transform_d(x):
     return ((2.0*4.0/sqrt(3.1415965358979))*exp(-(4.0**2)*(x-1.0)**2) + 0.5)
# def coord_transform(x):
#    return x
# def coord_transform_d(x):
#    return 1.0

class GP_Gaussian_kernel:

   def copy(self):
      return GP_Gaussian_kernel()

   # non-periodic kernels
   def ffkernel(self, x1,x2,fvar,l2):
       k = fvar*exp(-0.5*(x2-x1)**2/l2)
       return k
   def dfkernel(self, x1,x2,fvar,l2):
       k = (x2-x1)/l2*self.ffkernel(x1,x2,fvar,l2)
       return k
   def fdkernel(self, x1,x2,fvar,l2):
       k = -(x2-x1)/l2*self.ffkernel(x1,x2,fvar,l2)
       return k
   def ddkernel(self, x1,x2,fvar,l2): 
       k = (1-(x1-x2)**2/l2)/l2*self.ffkernel(x1,x2,fvar,l2)
       return k

class GP:

   import numpy as np

   def __init__(self, kernel=None):
      if (kernel is not None):
	 self.kernel = kernel.copy()
      else:
	 self.kernel = None

   def teach(self, kernel=None, func_r = None, func_val = None, func_noise = None, 
                   deriv_r = None, deriv_val = None, deriv_noise = None, 
		   len_scale=0.0, fvar=0.0):

      if (func_r is None):
	 func_r = GP.np.array([])
      if (func_val is None):
	 func_val = GP.np.array([])
      if (func_noise is None):
	 func_noise = GP.np.array([])
      if (deriv_r is None):
	 deriv_r = GP.np.array([])
      if (deriv_val is None):
	 deriv_val = GP.np.array([])
      if (deriv_noise is None):
	 deriv_noise = GP.np.array([])

      if (kernel is None):
	 if (self.kernel is None):
	    sys.stderr.write("GP.teach got no kernel, and no kernel is set")
	    sys.exit(1)
      else:
	 self.kernel = kernel.copy()

      if ((len(func_r) != len(func_val)) or (len(func_r) != len(func_noise))):
	 sys.stderr.write("mismatching func_r, func_val, func_noise, sizes in teach\n")
	 sys.exit(1)
      if ((len(deriv_r) != len(deriv_val)) or (len(deriv_r) != len(deriv_noise))):
	 sys.stderr.write("mismatching deriv_r, deriv_val, deriv_noise, sizes in teach\n")
	 sys.exit(1)
      if ((len_scale <= 0.0) or (fvar <= 0.0)):
	 sys.stderr.write("invalid len_scale or fvar\n")
	 sys.exit(1)

      self.l2 = len_scale**2
      self.fvar = fvar

      self.n_f = len(func_r)
      self.n_d = len(deriv_r)
      self.n_teach = self.n_f + self.n_d
      self.func_r = func_r.copy()
      self.deriv_r = deriv_r.copy()
      self.Cmat = GP.np.empty((self.n_teach, self.n_teach))
      # f vs. f, d
      for i_f in range(self.n_f):
	 self.Cmat[i_f,0:self.n_f]                 = kernel.ffkernel(func_r[i_f], func_r[0:self.n_f], self.fvar, self.l2)
	 self.Cmat[i_f,self.n_f:self.n_f+self.n_d] = kernel.fdkernel(func_r[i_f], deriv_r[0:self.n_d], self.fvar, self.l2)
	 self.Cmat[i_f,i_f] += func_noise[i_f]
      # d vs. f, d
      for i_d in range(self.n_d):
	 self.Cmat[self.n_f+i_d,0:self.n_f]                 = kernel.dfkernel(deriv_r[i_d], func_r[0:self.n_f], self.fvar, self.l2)
	 self.Cmat[self.n_f+i_d,self.n_f:self.n_f+self.n_d] = kernel.ddkernel(deriv_r[i_d], deriv_r[0:self.n_d], self.fvar, self.l2)
	 self.Cmat[self.n_f+i_d,self.n_f+i_d] += deriv_noise[i_d]
      concat_func_deriv_val = concatenate((func_val, deriv_val))
      self.vec = linalg.solve(self.Cmat, concat_func_deriv_val)
      return self

   def predict_f(self, r):
      k = concatenate(( self.kernel.ffkernel(self.func_r[0:self.n_f], r,  self.fvar, self.l2),  
                        self.kernel.dfkernel(self.deriv_r[0:self.n_d], r, self.fvar, self.l2) ))
      est = dot(k, self.vec)
      return est

   def predict_var(self, r):
      k = concatenate(( self.kernel.ffkernel(self.func_r[0:self.n_f], r,  self.fvar, self.l2),  
                        self.kernel.dfkernel(self.deriv_r[0:self.n_d], r, self.fvar, self.l2) ))
      Cmat_inv_k = linalg.solve(self.Cmat, k)
      est = self.fvar - dot(k, Cmat_inv_k)
      return est

def bracket_zero(f, min, step):
   x0 = min
   f0 = f(min)
   xc = min+step
   fc = f(xc)
   if (abs(fc) > abs(f0)):
      step *= -1.0
   i = 1
   while (i < 100):
      xc = min+i*step
      fc = f(xc)
      if (f0*fc < 0.0):
	 break
      x0 = xc
      f0 = fc
      i += 1
   if (i == 100):
      sys.stderr.write("Failed to bracket in 100 steps\n")
      sys.exit(1)
   return (x0, f0, xc, fc)

def bisect_zero(f, xmin, xmax, tol):
   fmin = f(xmin)
   fmax = f(xmax)
   if (fmax < fmin):
      t_f = fmax
      t_x = xmax
      fmax = fmin
      xmax = xmin
      fmin = t_f
      xmin = t_x
   while (abs(xmax-xmin) > tol):
      xmid = (xmin+xmax)/2.0
      fmid = f(xmid)
      if (fmid > 0):
	 xmax = xmid
	 fmax = fmid
      else:
	 xmin = xmid
	 fmin = fmid
   return (xmid, fmid)

transformed_len_scale = len_scale*coord_transform_d(min_collective+len_scale/2.0)
fvar = variance_prior

sys.stderr.write("teach d_from_d_gp\n")
d_from_d_gp = GP().teach(kernel=GP_Gaussian_kernel(),
   func_r = coord_transform(collective_vec), func_val = dE_dcollective_vec, func_noise = noise_variance_vec, 
   len_scale=transformed_len_scale, fvar = variance_prior)
sys.stderr.write("done\n")

(x0, f0, xc, fc) = bracket_zero( (lambda x : d_from_d_gp.predict_f(coord_transform(x))), min_bracket, (max_collective-min_collective)/float(n_collective))
print "# zero bracket between f(%f)=%f f(%f)=%f" % (x0, f0, xc, fc)

(xmin, fmin) = bisect_zero( (lambda x : d_from_d_gp.predict_f(coord_transform(x))), x0, xc, 0.00001)
print "# zero found at f(%f)=%f" % (xmin, fmin)

# learn function of transformed coord from transformed gradients
sys.stderr.write("teach func_from_d_gp\n")
func_from_d_gp = GP().teach(kernel=GP_Gaussian_kernel(),
                            func_r = coord_transform( array( [ xmin ] ) ), 
			    func_val = array( [ 0.0 ] ),
			    func_noise = array ( [ 0.0 ] ),
			    deriv_r = coord_transform(collective_vec),
			    deriv_val = dE_dcollective_vec/coord_transform_d(collective_vec),
			    deriv_noise = noise_variance_vec/(coord_transform_d(collective_vec)**2), 
   len_scale=transformed_len_scale, fvar = variance_prior)

# print func_from_d_gp.Cmat

print "#UI coll_val transformed_coll_val     gp_deriv gp_deriv_var     gp_func gp_func_var    numerical_int"
# print out GP gradient and function
integral_accum = 0.0
deriv_prev = d_from_d_gp.predict_f(coord_transform(min_collective))
for i in range(n_collective+1):
   # sys.stderr.write("%d\n" % i)
   collective_val = min_collective + float(i)/float(n_collective)*(max_collective-min_collective)
   deriv_cur = d_from_d_gp.predict_f(coord_transform(collective_val))
   integral_accum += ((deriv_cur+deriv_prev)/2.0*(max_collective-min_collective)/float(n_collective))
   print "UI %f %f    %f %f    %f %f     %f" % (collective_val, coord_transform(collective_val), deriv_cur, d_from_d_gp.predict_var(coord_transform(collective_val)), 
      func_from_d_gp.predict_f(coord_transform(collective_val)), func_from_d_gp.predict_var(coord_transform(collective_val)), integral_accum )
   deriv_prev = deriv_cur

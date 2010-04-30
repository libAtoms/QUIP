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

from pyatoms import *
from numpy import *
from math import *
from pylab import *

def stress_band_avg(at,component,dx=10.0):
   xr = arange(at.pos[:,0].min(),at.pos[:,0].max(),dx)
   res = []
   for i in range(len(xr)-1):
      band = logical_and(at.pos[:,0] > xr[i], at.pos[:,0] < xr[i+1])
      res.append(getattr(at,component)[band].mean())
   return xr, res


def scan_lines(at):
   order = range(at.n)
   order.sort(key=lambda i: (at.pos[i,1],at.pos[i,0]))

   # break up order into scan lines
   scans = []

   while order:
      scans.append([])
      while len(order >=2) and order[0] < order[1]:
         scans[-1].append(order.pop(0))

def stress_line(at,component):
   tipline = where(abs(at.pos[:,1]) < 3.0)[0]
   z = zip(at.pos[tipline,0],getattr(at,component)[tipline])

   # order by x coordinate
   z.sort()
   
   x, y = array([p[0] for p in z]), array([p[1] for p in z])

   # align max value with x=0
   x = x - x[y.argmax()]
   
   return x, y

def interp_stress_line(at, x, component):
   x0, y0 = stress_line(at, component)
   return interp(x,x0,y0)


def K_field_line(at, r, s):
   # K field solution for plane strain
   # s is compliance matrix, inverse of elastic constant matrix C_ij

   def sigma_xx(r,t):
      return 1.0/sqrt(2.0*pi*r)*cos(t/2.0)*(1.0 - sin(t/2.0)*sin(3.0*t/2.0))

   def sigma_yy(r,t):
      return 1.0/sqrt(2.0*pi*r)*cos(t/2.0)*(1.0 + sin(t/2.0)*sin(3.0*t/2.0))

   def sigma_zz(r,t,nu):
      return nu*(sigma_xx(r,t) + sigma_yy(r,t))

   def sigma_xy(r,t):
      return 1.0/sqrt(2.0*pi*r)*sin(t/2.0)*cos(t/2.0)*cos(3.0*t/2.0)

   nu = at.params['PoissonRatio']
   E = at.params['YoungsModulus']
   Ep = E#/(1.0-nu*nu)
   k = 3.0 - 4.0*nu
   mu = 1.0
   K_I = sqrt(at.params['G']*1.0e9*Ep)
   print 'K_I = %f' % K_I

   sigs = []
   epss = []

   for x in r:
      sig = zeros(6,'d')
      sig[0] = sigma_xx(x*1e-10,0.0)*K_I/1e9
      sig[1] = sigma_yy(x*1e-10,0.0)*K_I/1e9
      sig[2] = sigma_zz(x*1e-10,0.0,nu)*K_I/1e9
      sig[3] = sigma_xy(x*1e-10,0.0)*K_I/1e9
      sig[4] = 0.0
      sig[5] = 0.0

      eps = dot(s,sig)
      sigs.append(sig)
      epss.append(eps)

   sigs = array(sigs)
   epss = array(epss)

   # plane strain: exactly zero eps_zz, eps_xz, eps_yz
   epss[:,2] = 0.0
   epss[:,4] = 0.0
   epss[:,5] = 0.0
   
   return sigs, epss

def extrap_stress_line(atoms,x):

   recip_heights = [ 1.0/(at.pos[:,1].max()-at.pos[:,1].min()) for at in atoms ]
   val = {}
   extrap = {}
   
   for component in ('Sig_xx', 'Sig_yy', 'Sig_zz', 'Sig_xy', 'Sig_xz', 'Sig_yz', 
                     'S_xx_sub1', 'S_yy_sub1', 'S_zz_sub1', 'S_xy', 'S_xz', 'S_yz'):

      val[component] = array([ interp_stress_line(at, x, component) for at in atoms ])

      # fit straight line to (recip_heights, stress) data at each point

      extrap[component] = zeros(x.shape)
      A = ones((len(recip_heights),2))
      A[:,0] = recip_heights
      for p in range(x.size):
         coeffs, resids, rank, s = linalg.lstsq(A,val[component][:,p])
         # best fit line is coeffs[0]*x + coeffs[1], so
         # value at 1/h = 0 is coeffs[1]
#         print component, recip_heights, val[component][:,p], 'resids = ', resids
         extrap[component][p] = coeffs[1]

   sigs = array(zip(extrap['Sig_xx'], extrap['Sig_yy'], extrap['Sig_zz'],
                    extrap['Sig_xy'], extrap['Sig_xz'], extrap['Sig_yz']))
   epss = array(zip(extrap['S_xx_sub1'], extrap['S_yy_sub1'], extrap['S_zz_sub1'],
                    extrap['S_xy'], extrap['S_xz'], extrap['S_yz']))

   return sigs, epss

              
def mix_fields(xr, f1, f2, R1, R2):
   res = zeros(xr.shape)
   nf1 = zeros(xr.shape)
   nf1[:f1.size] = f1

   res[xr < R1] = nf1[xr < R1]
   
   for i in where(logical_and(xr >= R1, xr <= R2))[0]:
      f = (xr[i]-R1)/(R2-R1)
      res[i] = nf1[i]*(1.0-f) + f2[i]*f

   res[xr > R2] = f2[xr > R2]
   return res

def elastic_constant_matrix(C11,C12,C44):
   c = zeros((6,6),'d')
   c[0,0] = c[1,1] = c[2,2] = C11
   c[0,1] = c[1,0] = c[0,2] = c[2,0] = c[1,2] = c[2,1] = C12
   c[3,3] = c[4,4] = c[5,5] = C44
   return c

def load_atoms():
   us = [ Atoms('cl_unrec_%03d_elastic.xyz' % i) for i in (400,500,600,1000)]
   rs = [ Atoms('cl_recon_%03d_elastic.xyz' % i) for i in (400,500,600,1000)]
   return us, rs

def calc_fields(us, rs):
   c = elastic_constant_matrix(C11=151.3705973705, C12=76.395128559, C44=56.4263274427)
   s = linalg.inv(c)

   xr = arange(1.0,100000.0,1.0)

   usigs, uepss = extrap_stress_line(us,xr)
   rsigs, repss = extrap_stress_line(rs,xr)
   ksigs, kepss = K_field_line(us[0], xr, s)

   # hack - copy sig_yy over sig_xx and eps_yy over eps_xx
   usigs[:,0] = usigs[:,1]
   rsigs[:,0] = rsigs[:,1]
   uepss[:,0] = uepss[:,1]
   repss[:,0] = repss[:,1]

   umixsigs = []; rmixsigs = []
   for usig, rsig, ksig in zip(usigs.transpose(), rsigs.transpose(), ksigs.transpose()):
      umixsigs.append(mix_fields(xr, usig, ksig, 400, 700))
      rmixsigs.append(mix_fields(xr, rsig, ksig, 400, 700))
   umixsigs = array(umixsigs).transpose()
   rmixsigs = array(rmixsigs).transpose()

   umixepss = []; rmixepss = []
   for ueps, reps, keps in zip(uepss.transpose(), repss.transpose(), kepss.transpose()):
      umixepss.append(mix_fields(xr, ueps, keps, 400, 700))
      rmixepss.append(mix_fields(xr, reps, keps, 400, 700))
   umixepss = array(umixepss).transpose()
   rmixepss = array(rmixepss).transpose()

   write_file('unrec_10um_nm_res.dat', xr, umixepss)
   write_file('recon_10um_nm_res.dat', xr, rmixepss)
   return usigs, uepss, rsigs, repss, ksigs, kepss, umixsigs, umixepss, rmixsigs, rmixepss


def write_file(fname, xr, mixepss):
   f = open(fname, 'w')
   for x, eps, in zip(xr, mixepss)[::10]:
      f.write("%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f%16.8f\n" %
              (x/10.0,eps[0],eps[1],eps[2],eps[3],eps[4],eps[5]))

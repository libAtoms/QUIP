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

# K-means clustering

from numpy import *

def kmeans(data, K, means=None):
   data = array(data)
   N, M = data.shape

   if means is None:
      # Initialise centroids randomly within range of input data
      means = array([[random.uniform(low=data[:,i].min(), high=data[:,i].max()) for i in range(M)]
                     for j in range(K) ])
   else:
      means = array(means)

   assign = zeros(N, int)
   change = True
   while change:
      
      change = False

      # E-step
      for n in range(N):
         dmin = 9.99e99
         for k in range(K):
            d = dot((data[n,:] - means[k,:]), (data[n,:] - means[k,:]))
            if d < dmin:
               dmin = d
               kmin = k
               
         if assign[n] != kmin:
            assign[n] = kmin
            change = True

      # M-step
      means[:] = 0.0
      for k in range(K):
         if any(assign == k):
            means[k,:] = data[assign == k,:].mean(axis=0)

   err = 0.0
   for n in range(N):
      err += dot(data[n,:] - means[assign[n],:],data[n,:] - means[assign[n],:])

   return assign, means, err

def global_kmeans(data, Kmax=None):
   N, M = data.shape

   assign = {}
   means = {}

   assign[1], means[1], err = kmeans(data,1)
   kmin = 1
   errmin_k = err

   if Kmax is None:
      Kmax = N

   print 'K', 1, 'errmin_k', errmin_k
   print means
   
   for K in range(2,Kmax+1):

      errmin_n = 9.99e99
      nmeans = zeros((K,M))
      nmeans[0:K-1,:] = means[K-1] # copy from optimal solution for K-1
      
      for n in range(N):
         nmeans[K-1,:] = data[n,:]
         print 'K', K, 'n', n

         print 'nmeans'
         print nmeans

         tmp_assign, tmp_means, err = kmeans(data,K=K,means=nmeans)

         print 'tmp_means'
         print tmp_means
         
         print 'err', err
         print
         if err < errmin_n:
            errmin_n = err
            nmin = n
            assign[K] = tmp_assign.copy()
            means[K]  = tmp_means.copy()

      print 'K', K, 'nmin', nmin, 'err', errmin_n     
      if errmin_n < errmin_k:
         errmin_k = errmin_n
         kmin = K
         
   print 'kmin', kmin, 'errmin_k', errmin_k

   return kmin, assign[kmin], means[kmin], errmin_k

def kmeans_plot(data, K, assign, means):

   colours = ['r','g','b','c','m','y','k']

   for k in range(K):
      if not any(assign == k): continue
      scatter(data[assign == k,0], data[assign == k,1], c=colours[k], marker='o')
      scatter(means[k,0], means[k,1], c=colours[k], marker='s')

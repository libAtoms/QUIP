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

# Detection and classification of local minima

def valid_index(a, x):
   try:
      a[x]
      return True
   except IndexError:
      return False

def neighbours(a, x, avoid_zero=False):
   nvec = [[1,0,0],
           [0,1,0],
           [0,0,1],
           [-1,0,0],
           [0,-1,0],
           [0,0,-1]]

   for d in nvec:
      idx = tuple([i+j for i,j in zip(x,d)])
      if valid_index(a, idx):
         if avoid_zero and a[idx] == 0: continue
         yield idx


def steepest_descent(a, x0, avoid_zero=True):
      
   t = a[x0]
   change = True
   while change:
      change = False

      for x in neighbours(a, x0):
         t1 = a[x]
         if avoid_zero and t1 == 0: continue
         if t1 < t:
            x0 = x
            t  = t1
            change = True

   return x0


def find_local_minima(a, pot_zero=None, shallow_threshold=None):
   if pot_zero is None:
      pot_zero = a.max()

   local_minima = []
   basin_volume = []
   a_trans = fzeros(a.shape,int)

   for x0 in transpose((a != 0).nonzero()):
      min = steepest_descent(a, tuple(x0))
      if min not in local_minima:
         local_minima.append(min)
         basin_volume.append(0)

      a_trans[tuple(x0)] = a[min]
      basin_volume[local_minima.index(min)] += 1

   total_volume = float((a != 0).sum())
   
   if shallow_threshold is not None:
      for min in local_minima:
         well_depth = pot_zero - a[min]

         print min, well_depth, basin_volume[local_minima.index(min)], basin_volume[local_minima.index(min)]/total_volume
         
         basin_volume[local_minima.index(min)] = basin_volume[local_minima.index(min)]*well_depth


      local_minima = array(local_minima)
      basin_volume = array(basin_volume,float)

      #mask = (basin_volume - mean(basin_volume))/std(basin_volume) > shallow_threshold
      #local_minima = [ tuple(min) for min in local_minima[mask] ]

   return local_minima, a_trans


def metropolis(a, x0, kmax, T0=1.0, RT=0.85, avoid_zero=True):

   class Reject: pass

   x = x0
   k = 0
   n_accept = 0
   T = T0
   m = 5*len(x0)

   for k in range(kmax):

      cand = random.choice(list(neighbours(a, x, avoid_zero)))
      try:

         if a[cand]-a[x] < 0:
            # always go downhill
            x = cand
         else:
            # go uphill with probability exp(-(e2-e1)/kT)

            met = exp(-(a[cand]-a[x])/T)
            if rand() < met:
               x = cand
            else:
               raise Reject

         n_accept = n_accept + 1

      except Reject:
         pass

      if k % m == 0:
         T = T0*(1.0 - float(k)/kmax)**4.0
         print 'temp', T
      

   print 'acceptance rate', float(n_accept)/kmax
   return x
         

def metropolis_local_minima(a,nsample):
   local_minima = {}

   n = 0
   while n < nsample:

      x0 = (randint(a.shape[0])+1,randint(a.shape[1])+1,randint(a.shape[2])+1)
      if not valid_index(a, x0) or a[x0] == 0: continue

      n += 1
      
      min = metropolis(a, tuple(x0), 1000, 2.0, 0.99)
      
      if min not in local_minima:
         local_minima[min] = 0
         
      local_minima[min] += 1

   ordered_minima = sorted([(v,k) for (v,k) in zip(local_minima.values(), local_minima.keys())])

   print ordered_minima

   class Duplicate: pass

   while True:
      try:
         for (freq1,min1),(freq2,min2) in combinations(ordered_minima,2):
            if (farray(min1) - farray(min2)).norm() < sqrt(5.0)+0.1:
               print 'duplicate', (freq1, min1), (freq2, min2)
               raise Duplicate
         else:
            break
      except Duplicate:
         if freq1 < freq2:
            ordered_minima.remove((freq1,min1))
            ordered_minima[ordered_minima.index((freq2,min2))] = (freq1+freq2, min2)
         else:
            ordered_minima.remove((freq2,min2))
            ordered_minima[ordered_minima.index((freq1,min1))] = (freq1+freq2, min1)

   ordered_minima = sorted(ordered_minima)

   print ordered_minima

   res = []
   frac = 0.0
   while ordered_minima:
      freq, min = ordered_minima.pop()
      frac += float(freq)/nsample
      res.append(min)
      print res, frac
      if frac > 0.6: break

   return res
      

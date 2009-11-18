from quippy import *

def oxygenate(at):
   saved_cutoff, saved_use_uniform_cutoff = at.cutoff, at.use_uniform_cutoff
   at.set_cutoff_factor(1.2)
   at.calc_connect()

   add_pos = []
   rem_list = []
   
   for i in frange(at.n):
      if at.z[i] != 14: continue # Only consider silicon atoms

      neighb = at.neighbours[i]
      
      if len(neighb) == 4:
         continue
      elif len(neighb) == 3:
         # add single O atom to complete SiO4 tetrahedron
         
         p1 = neighb[1].diff - neighb[2].diff
         p2 = neighb[1].diff - neighb[3].diff

         length = (neighb[1].distance + neighb[2].distance + neighb[3].distance)/3.0

         n = cross(p1, p2)
         n = n/sqrt(dot(n,n))

         if dot(n,neighb[1].diff) > 0.0: n = -n

         add_pos.append(at.pos[:,i]+length*n)
      elif len(neighb) == 2:
         # add two O atoms to complete SiO4 tetrahedron

         length = (neighb[1].distance + neighb[2].distance)/2.0

         n1 = cross(neighb[1].diff, neighb[2].diff)
         n1 = n1/sqrt(dot(n1,n1))

         n2 = neighb[2].diff - neighb[1].diff
         n2 = n2/sqrt(dot(n2,n2))

         n3 = cross(n1,n2)

         o1 = at.pos[:,i] + length*( n1*sqrt(2./3)+n3*sqrt(1./3))
         o2 = at.pos[:,i] + length*(-n1*sqrt(2./3)+n3*sqrt(1./3))

         add_pos.append(o1)
         add_pos.append(o2)
      elif len(neighb) <= 1:
         rem_list.append(i)
         

   if len(add_pos) > 0:
      add_z = [8]*len(add_pos)
      add_pos = farray(add_pos)
      at.add_atoms(add_pos, add_z)

   if len(rem_list) > 0:
      at.remove_atoms(rem_list)

   at.cutoff, at.use_uniform_cutoff = saved_cutoff, saved_use_uniform_cutoff
   

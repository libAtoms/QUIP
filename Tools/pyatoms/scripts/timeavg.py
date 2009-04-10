#!/usr/bin/env python

# read in a bunch of frames and do a linear average of all their real properties

from pyatoms import *

def temperature(a):
   return a.v2.mean()*a.mass[0]/(3.0*BOLTZMANN_K)

def add_speed(frames):
   for a in frames:
      a.add_property('v',norm(a.velo))
      a.add_property('v2',norm2(a.velo))
      yield a

def timeavg(frames):
   first = frames.next()
   avg = first.copy()

   n = 1
   for a in frames:
      print n, temperature(a), temperature(avg)
      avg.real[:] = (n*avg.real + a.real)/float(n+1.0)
      n = n + 1

   return avg

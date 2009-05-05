from numpy import *
from quippy import *
import quippy._quippy

def f(x):
   return sqrt(dot((x-1.0),(x-1.0)))

def df(x):
   return [1.0]

def hook(x,dx,e,done,do_print,data):
   pass

a = array([0.0,0.0,0.0],order='F')

quippy._quippy.minmod.minim(a,f,df,"cg",1e-3,10,"FAST_LINMIN",hook)


#!/usr/bin/env python

# Generate a Monkhorst-Pack k-point mesh for QUIP, using CASTEP
# to find the irreducible k-point set.

from pyatoms import *
import sys, os

if len(sys.argv) not in  (5,8):
   print 'Usage: %s atoms.xyz nx ny nz [shiftx shifty shiftz]'
   sys.exit(1)

# Command used to execute castep, with a %s where seed name should go
if os.environ.has_key('CASTEP'):
   CASTEP = os.environ['CASTEP']
else:
   CASTEP = './castep %s'

# If there's no %s, put seed name at end of string
if CASTEP.find('%s') == -1:
   CASTEP = CASTEP + ' %s'

nx,ny,nz = map(int, sys.argv[2:5])

if len(sys.argv) == 8:
   sx,sy,sz = map(float, sys.argv[5:8])
else:
   sx = sy = sz = 0.0

cell_lines = ['KPOINT_MP_GRID   %d %d %d\n' % (nx, ny, nz),
              'KPOINT_MP_OFFSET %f %f %f\n' % (sx, sy, sz)]

param_lines = ['task : SinglePoint\n']

cell = castep.CastepCell(cell_lines)
param = castep.CastepParam(param_lines)
a = Atoms(sys.argv[1])
cell.update_from_atoms(a)

cell.write('tmp.cell')
param.write('tmp.param')

os.system((CASTEP % 'tmp') + ' -dryrun')

castep_lines = open('tmp.castep').readlines()

try:
   kp_start = castep_lines.index('             +  Number       Fractional coordinates        Weight  +\n')
except:
   raise ValueError('No k-points found in CASTEP output file')
   sys.exit(1)

i = kp_start + 2
kp_lines = []
while True:
   line = castep_lines[i]
   if line.strip().startswith('+++++++'): break
   kp_lines.append(line)
   i += 1

kpoints = []
for line in kp_lines:
   p1, n, kx, ky, kz, weight, pt = line.split()
   kpoints.append(map(float, (kx,ky,kz,weight)))

os.remove('tmp.cell')
os.remove('tmp.param')
os.remove('tmp.castep')

print '<KPoints N="%d">' % len(kpoints)

for kx, ky, kz, weight in kpoints:
   print '  <point weight="%d"> %.5f %.5f %.5f </point>' % (int(1.0/weight), kx, ky, kz)

print '</KPoints>'

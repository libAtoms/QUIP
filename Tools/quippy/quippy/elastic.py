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

# Elastic constant calculation.

# Code adapted from elastics.py script, available from
# http://github.com/djw/elastic-constants
#
# Copyright (c) 2008, Dan Wilson
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the copyright holder nor the
#       names of its contributors may be used to endorse or promote products
#       derived from this software without specific prior written permission.
# 
# THIS SOFTWARE IS PROVIDED BY DAN WILSON ''AS IS'' AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL DAN WILSON BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from scipy import stats
try:
   import pylab
except  ImportError:
   pass

from farray import *

def strain_matrix(strain_vector):
   e1, e2, e3, e4, e5, e6 = strain_vector
   return farray([[1.0+e1, 0.5*e6, 0.5*e5],
                  [0.5*e6, 1.0+e2, 0.5*e4],
                  [0.5*e5, 0.5*e4, 1.0+e3]])

def stress_matrix(stress_vector):
   s1, s2, s3, s4, s5, s6 = stress_vector
   return farray([[s1, s6, s5],
                  [s6, s2, s4],
                  [s5, s4, s3]])

def strain_vector(strain_matrix):
   return farray([strain_matrix[1,1] - 1.0,
                  strain_matrix[2,2] - 1.0,
                  strain_matrix[3,3] - 1.0,
                  2.0*strain_matrix[2,3],
                  2.0*strain_matrix[1,3],
                  2.0*strain_matrix[1,2]])

def stress_vector(stress_matrix):
   return farray([stress_matrix[1,1],
                  stress_matrix[2,2],
                  stress_matrix[3,3],
                  stress_matrix[2,3],
                  stress_matrix[1,3],
                  stress_matrix[1,2]])



Cij_symmetry = {
   'cubic':           farray([[1, 7, 7, 0, 0, 0],
                              [7, 1, 7, 0, 0, 0],
                              [7, 7, 1, 0, 0, 0],
                              [0, 0, 0, 4, 0, 0],
                              [0, 0, 0, 0, 4, 0],
                              [0, 0, 0, 0, 0, 4]]),
      
   'trigonal_high':   farray([[1, 7, 8, 9, 10, 0],
                              [7, 1, 8, 0,-9, 0],
                              [8, 8, 3, 0, 0, 0],
                              [9, -9, 0, 4, 0, 0],
                              [10, 0, 0, 0, 4, 0],
                              [0, 0, 0, 0, 0, 6]]),

   'trigonal_low':    farray([[1,  7,  8,  9,  10,  0 ],
                              [7,  1,  8, -9, -10,  0 ],
                              [8,  8,  3,  0,   0,  0 ],
                              [9, -9,  0,  4,   0, -10],
                              [10,-10, 0,  0,   4,  9 ],
                              [0,  0,  0, -10 , 9,  6 ]]),
   
   'tetragonal_high': farray([[1, 7, 8, 0, 0, 0],
                              [7, 1, 8, 0, 0, 0],
                              [8, 8, 3, 0, 0, 0],
                              [0, 0, 0, 4, 0, 0],
                              [0, 0, 0, 0, 4, 0],
                              [0, 0, 0, 0, 0, 6]]),
   
   'tetragonal_low':  farray([[1, 7, 8, 0, 0, 11],
                              [7, 1, 8, 0, 0, -11],
                              [8, 8, 3, 0, 0, 0],
                              [0, 0, 0, 4, 0, 0],
                              [0, 0, 0, 0, 4, 0],
                              [11, -11, 0, 0, 0, 6]]),

   'orthorhombic':    farray([[ 1,  7,  8,  0,  0,  0],
                              [ 7,  2, 12,  0,  0,  0],
                              [ 8, 12,  3,  0,  0,  0],
                              [ 0,  0,  0,  4,  0,  0],
                              [ 0,  0,  0,  0,  5,  0],
                              [ 0,  0,  0,  0,  0,  6]]),

   'monoclinic':      farray([[ 1,  7,  8,  0,  10,  0],
                              [ 7,  2, 12,  0, 14,  0],
                              [ 8, 12,  3,  0, 17,  0],
                              [ 0,  0,  0,  4,  0,  20],
                              [10, 14, 17,  0,  5,  0],
                              [ 0,  0,  0, 20,  0,  6]]),

   'triclinic':       farray([[ 1,  7,  8,  9,  10, 11],
                              [ 7,  2, 12,  13, 14, 15],
                              [ 8, 12,  3,  16, 17, 18],
                              [ 9, 13, 16,  4,  19, 20],
                              [10, 14, 17, 19,  5,  21],
                              [11, 15, 18, 20,  21, 6 ]]),
   }


strain_patterns = {

   'cubic': [
      # strain pattern e1+e4, yields C11, C21, C31 and C44, then C12 is average of C21 and C31
      [ farray([1,0,0,1,0,0]), [(1,1), (2,1), (3,1), (4,4)]]
   ],

   'trigonal_high': [
      # strain pattern e3 yield C13, C23 and C33
      [ farray([0,0,1,0,0,0]), [(1,3), (2,3), (3,3)]],

      # strain pattern e1+e4 yields C11 C21 C31 and C44
      [ farray([1,0,0,1,0,0]), [(1,1), (2,1), (3,1), (4,4)]],

      # strain pattern e1 yields C11 C21 C31 C41 C51
      [ farray([1,0,0,0,0,0]), [(1,1), (2,1), (3,1), (4,1), (5,1)]],

      # strain pattern e3+e4
      [ farray([0,0,1,1,0,0]), [(3,3), (4,4)]]
      
   ],
   
   'trigonal_low': [
     # strain pattern e1, yields C11, C21, C31, C41, C51
     [ farray([1,0,0,0,0,0]), [(1,1), (2,1), (3,1), (4,1), (5,1)]],
   
     # strain pattern e3 + e4, yields C33, C44
     [ farray([0,0,1,1,0,0]), [(3,3), (4,4)] ],

     [ farray([0,0,0,0,0,1]), [(6,6)] ]
   ],

   'tetragonal': [
     # strain pattern e1+e4
     [ farray([1,0,0,1,0,0]), [(1,1), (2,1), (3,1), (6,1), (4,4)] ],

     # strain pattern e3+e6
     [ farray([0,0,1,0,0,1]), [(3,3), (6,6)] ]
   ],

   'orthorhombic': [
      # strain pattern e1+e4
      [ farray([1,0,0,1,0,0]), [(1,1), (2,1), (3,1), (4,4)] ],

      # strain pattern e2+e5
      [ farray([0,1,0,0,1,0]), [(1,2), (2,2), (3,2), (5,5)] ],

      # strain pattern e3+e6
      [ farray([0,0,1,0,0,1]), [(1,3), (2,3), (3,3), (6,6)] ]
   ],

   'monoclinic': [
      # strain pattern e1+e4
      [ farray([1,0,0,1,0,0]), [(1,1), (2,1), (3,1), (4,4), (5,1), (6,4)] ],

      # strain pattern e3+e6
      [ farray([0,0,1,0,0,1]), [(1,3), (2,3), (3,3), (5,3), (4,6), (6,6)] ],

      # strain pattern e2
      [ farray([0,1,0,0,0,0]), [(1,2), (2,2), (3,2), (5,2)] ],

      # strain pattern e5
      [ farray([0,0,0,0,1,0]), [(1,5), (2,5), (3,5), (5,5)] ]
   ],

   'triclinic': [
      [ farray([1,0,0,0,0,0]), [(1,1), (2,1), (3,1), (4,1), (5,1), (6,1)]],
      [ farray([0,1,0,0,0,0]), [(1,2), (2,2), (3,2), (4,2), (5,2), (6,2)]],
      [ farray([0,0,1,0,0,0]), [(1,3), (2,3), (3,3), (4,3), (5,3), (6,3)]],
      [ farray([0,0,0,1,0,0]), [(1,4), (2,4), (3,4), (4,4), (5,4), (6,4)]],
      [ farray([0,0,0,0,1,0]), [(1,5), (2,5), (3,5), (4,5), (5,5), (6,5)]],
      [ farray([0,0,0,0,0,1]), [(1,6), (2,6), (3,6), (4,6), (5,6), (6,6)]],
   ]

   }

Cij_symmetry['hexagonal'] = Cij_symmetry['trigonal_high']
Cij_symmetry[None] = Cij_symmetry['triclinic']

strain_patterns['hexagonal'] = strain_patterns['trigonal_high']
strain_patterns['tetragonal_high'] = strain_patterns['tetragonal_low'] = strain_patterns['tetragonal']
strain_patterns[None] = strain_patterns['triclinic']


def generate_strained_configs(at0, symmetry='triclinic', N_steps=5, delta=1e-2):
   """Generate a sequence of strained configurations"""

   if not symmetry in strain_patterns:
      raise ValueError('Unknown symmetry %s. Valid options are %s' % (symmetry, strain_patterns.keys()))

   for pindex, (pattern, fit_pairs) in fenumerate(strain_patterns[symmetry]):
      for step in frange(N_steps):
         strain = numpy.where(pattern == 1, delta*(step-(N_steps+1)/2.0), 0.0)
         at = at0.copy()
         T = strain_matrix(strain)
         at.set_lattice(numpy.dot(T,at.lattice), scale_positions=False)
         at.pos[:] = numpy.dot(T,at.pos)
         at.params['strain'] = T
         yield at


def calc_stress(configs, pot, relax=False, relax_tol=1e-3, relax_steps=100):
   """Given a sequence of configs, calculate stress on each one"""
   from quippy import GPA
   for at in configs:
      at2 = at.copy()
      at2.set_cutoff(pot.cutoff())
      at2.calc_connect()
      if relax:
         pot.minim(at2, 'cg', relax_tol, relax_steps, do_pos=True, do_lat=False)
      pot.calc(at2, virial=True)
      at2.params['stress'] = -at2.params['virial']*GPA/at2.cell_volume()
      yield at2


def fit_elastic_constants(configs, symmetry=None, N_steps=5, verbose=True, graphics=True):
   """Given a sequence of configs with strain and stress parameters, fit elastic constants C_ij"""

   def do_fit(index1, index2, stress, strain, patt):
      if verbose:
         print 'Fitting C_%d%d' % (index1, index2)
         print 'Strain %r' % strain[:,index2]
         print 'Stress %r' % stress[:,index1]

      cijFitted,intercept,r,tt,stderr = stats.linregress(strain[:,index2],stress[:,index1])

      if verbose:
         # print info about the fit
         print     'Cij (gradient)          :    ',cijFitted
         print     'Error in Cij            :    ', stderr
         if abs(r) > 0.9:
            print 'Correlation coefficient :    ',r
         else:
            print 'Correlation coefficient :    ',r, '     <----- WARNING'

      if graphics:
         # position this plot in a 6x6 grid
         sp = pylab.subplot(6,6,6*(index1-1)+index2)
         sp.set_axis_on()

         # change the labels on the axes
         xlabels = sp.get_xticklabels()
         pylab.setp(xlabels,'rotation',90,fontsize=7)
         ylabels = sp.get_yticklabels()
         pylab.setp(ylabels,fontsize=7)

         # colour the plot depending on the strain pattern
         colourDict = {1: '#BAD0EF', 2:'#FFCECE', 3:'#BDF4CB', 4:'#EEF093',5:'#FFA4FF',6:'#75ECFD'}
         sp.set_axis_bgcolor(colourDict[patt])

         # plot the data
         pylab.plot([strain[1,index2],strain[-1,index2]],[cijFitted*strain[1,index2]+intercept,cijFitted*strain[-1,index2]+intercept])
         pylab.plot(list(strain[:,index2]),list(stress[:,index1]),'ro')

      return cijFitted, stderr

   if not symmetry in strain_patterns:
      raise ValueError('Unknown symmetry %s. Valid options are %s' % (symmetry, strain_patterns.keys()))
   
   # There are 21 independent elastic constants
   Cijs = {}
   Cij_err = {}

   # Construct mapping from (i,j) to index into Cijs in range 1..21
   # (upper triangle only to start with)
   Cij_map = {}
   Cij_map_sym = {}
   for i in frange(6):
      for j in frange(i,6):
         Cij_map[(i,j)] = Cij_symmetry[None][i,j]
         Cij_map_sym[(i,j)] = Cij_symmetry[symmetry][i,j]

   # Reverse mapping, index 1..21 -> tuple (i,j) with i, j in range 1..6
   Cij_rev_map = dict(zip(Cij_map.values(), Cij_map.keys()))

   # Add the lower triangle to Cij_map, e.g. C21 = C12
   for (i1,i2) in Cij_map.keys():
      Cij_map[(i2,i1)] = Cij_map[(i1,i2)]
      Cij_map_sym[(i2,i1)] = Cij_map_sym[(i1,i2)]


   N_pattern = len(strain_patterns[symmetry])
   configs = iter(configs)

   strain = fzeros((N_pattern, N_steps, 6))
   stress = fzeros((N_pattern, N_steps, 6))

   if graphics:
      fig = pylab.figure(num=1, figsize=(9.5,8),facecolor='white')
      fig.clear()
      fig.subplots_adjust(left=0.07,right=0.97,top=0.97,bottom=0.07,wspace=0.5,hspace=0.5)
      
      for index1 in range(6):
         for index2 in range(6):
            # position this plot in a 6x6 grid
            sp = pylab.subplot(6,6,6*(index1)+index2+1)
            sp.set_axis_off()
            pylab.text(0.4,0.4, "n/a")

   # Fill in strain and stress arrays from config Atoms list
   for pindex, (pattern, fit_pairs) in fenumerate(strain_patterns[symmetry]):
      for step in frange(N_steps):
         at = configs.next()
         strain[pindex, step, :] = strain_vector(at.params['strain'])
         stress[pindex, step, :] = stress_vector(at.params['stress'])

   # Do the linear regression
   for pindex, (pattern, fit_pairs) in fenumerate(strain_patterns[symmetry]):
      for (index1, index2) in fit_pairs:
         fitted, err = do_fit(index1, index2, stress[pindex,:,:], strain[pindex,:,:], pindex)

         index = abs(Cij_map_sym[(index1, index2)])

         if not index in Cijs:
            if verbose:
               print 'Setting C%d%d (%d) to %f +/- %f' % (index1, index2, index, fitted, err)
            Cijs[index] = [fitted]
            Cij_err[index] = [err]
         else:
            if verbose:
               print 'Updating C%d%d (%d) with value %f +/- %f' % (index1, index2, index, fitted, err)
            Cijs[index].append(fitted)
            Cij_err[index].append(err)
         if verbose: print '\n'


   C = fzeros((6,6))
   C_err = fzeros((6,6))
   C_labels = fzeros((6,6),dtype='S4')
   C_labels[:] = '    '

   # Convert lists to mean
   for k in Cijs:
      Cijs[k] = numpy.mean(Cijs[k])

   # Combine statistical errors
   for k, v in Cij_err.iteritems():
      Cij_err[k] = numpy.sqrt(numpy.sum(farray(v)**2))/numpy.sqrt(len(v))

   if symmetry.startswith('trigonal'):
      # Special case for trigonal lattice: C66 = (C11 - C12)/2
      Cijs[Cij_map[(6,6)]] = 0.5*(Cijs[Cij_map[(1,1)]]-Cijs[Cij_map[(1,2)]])
      Cij_err[Cij_map[(6,6)]] = numpy.sqrt(Cij_err[Cij_map[(1,1)]]**2 + Cij_err[Cij_map[(1,2)]]**2)

   # Generate the 6x6 matrix of elastic constants 
   # - negative values signify a symmetry relation
   for i in frange(6):
      for j in frange(6):
         index = Cij_symmetry[symmetry][i,j]
         if index > 0:   
            C[i,j] = Cijs[index]
            C_err[i,j] = Cij_err[index]
            C_labels[i,j] = ' C%d%d' % Cij_rev_map[index]
            C_err[i,j] = Cij_err[index]
         elif index < 0:
            C[i,j] = -Cijs[-index]
            C_err[i,j] = Cij_err[-index]
            C_labels[i,j] = '-C%d%d' % Cij_rev_map[-index]

   if verbose:
      print numpy.array2string(C_labels).replace("'","")
      print '\n = \n'
      print numpy.array2string(C, suppress_small=True, precision=2)
      print

      # Summarise the independent components of C_ij matrix
      printed = {}
      for i in frange(6):
         for j in frange(6):
            index = Cij_symmetry[symmetry][i,j]
            if index <= 0 or index in printed: continue
            print 'C_%d%d = %-4.2f +/- %-4.2f GPa' % (i, j, C[i,j], C_err[i,j])
            printed[index] = 1

   return C, C_err


def elastic_constants(pot, at, sym='cubic', relax=True, verbose=True, graphics=True):
   """Convenience function which returns the 6x6 matrix :math:`C_{ij}` of configuration `at`
   with Potential `pot` assumsing symmetry `sym` (default "cubic")."""
   
   strained_configs = generate_strained_configs(at, sym)
   stressed_configs = calc_stress(strained_configs, pot, relax=relax)
   C, C_err = fit_elastic_constants(stressed_configs, sym, verbose=verbose, graphics=graphics)

   return C


def atomic_strain(at, r0, crystal_factor=1.0):
   """Atomic strain as defined by JA Zimmerman in `Continuum and Atomistic Modeling of
   Dislocation Nucleation at Crystal Surface Ledges`, PhD Thesis, Stanford University (1999)."""

   strain = fzeros((3,3,at.n))

   for l in frange(at.n):
      for n in at.neighbours[l]:
         for i in frange(3):
            for j in frange(3):
               strain[i,j,l] += n.diff[i]*n.diff[j]/r0**2.0
               
   return strain/crystal_factor

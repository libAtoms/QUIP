#!/usr/bin/env python
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

import sys
from numpy import *
from pylab import *

data = open('fitmethods-results').readlines()

qm_force = -4.58694

res = {}
for line in data:
   method, embed, buffer, dum1, time, dum2, force, sig, rms = line.split()
   embed, buffer = map(int, (embed,buffer))
   try:
      sig, time, force = map(float, (sig,time,force))
      if not (embed,buffer) in res:
         res[(embed,buffer)] = []
      res[(embed,buffer)].append((-(force-qm_force),sig/sqrt(time)))
   except ValueError:
      print 'skipping: %s' % line
         

for k in res.keys():
   res[k] = array(res[k])

ind = arange(5)
width=0.5

locs = (7,8,9,4,5,6,1,2,3)

clf()

colors=('red','green','blue','cyan','yellow')

ax = gca()

ps = []

for i, k in enumerate(sorted(res.keys())):
   subplot(3,3,locs[i])
   
   ps.append(bar(ind, res[k][:,0],yerr=res[k][:,1],color=colors))
#   title(str(k))
   ylim(-0.2,0.5)
   xlim(-width/2,len(ind)+width/2)
   gca().get_xaxis().tick_bottom()
   xticks(ind+width, ('AP', 'CM', 'FMA', 'FMS', 'FMSS'))


subplots_adjust(left=0.15,bottom=0.15)

figtext(0.53,0.03,'Buffer Hops', size='x-large', horizontalalignment='center')
figtext(0.25,0.07,'3', size='large', horizontalalignment='center')
figtext(0.53,0.07,'4', size='large', horizontalalignment='center')
figtext(0.8,0.07,'5', size='large', horizontalalignment='center')


figtext(0.03,0.53,'Embed Hops', size='x-large', verticalalignment='center', rotation='vertical')
figtext(0.07,0.25,'2', size='large', verticalalignment='center', rotation='vertical')
figtext(0.07,0.53,'3', size='large', verticalalignment='center', rotation='vertical')
figtext(0.07,0.8, '4', size='large', verticalalignment='center', rotation='vertical')



show()

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

import sys, itertools, os.path

from numpy import *
from pylab import *
from pyatoms import *

def ema(s, n):
    """
    returns an n period exponential moving average
    s is a list ordered from oldest (index 0) to most recent (index
    -1) n is an integer
    returns a numeric array of the exponential moving average
 """
    s = array(s)
    ema = []
    j = 1
    #get n sma first and calculate the next n period ema
    sma = sum(s[:n]) / n
    multiplier = 2 / float(1 + n)
    ema.append(sma)
    #EMA(current) = ( (Price(current) - EMA(prev) ) xMultiplier) + EMA(prev)
    ema.append(( (s[n] - sma) * multiplier) + sma)
    #now calculate the rest of the values
    for i in s[n+1:]:
        tmp = ( (i - ema[j]) * multiplier) + ema[j]
        j = j + 1
        ema.append(tmp)
    return ema


def max_bolt(m,T,v):
   "Maxwell-Boltmann distribution of speeds at temperature T for particles of mass m"
   return 4*pi*(m/(2*pi*BOLTZMANN_K*T))**(3.0/2.0)*(v**2)*exp(-m*v**2/(2*BOLTZMANN_K*T))

def temperature(a,mask=None):
    if mask is None:
        mask = zeros((a.n,),bool)
        mask[:] = True
    assert all(a.mass == a.mass[0]) # only works for monatomic systems
    return norm2(a.velo[mask]).mean()*a.mass[0]/(3.0*BOLTZMANN_K)

def velo_hist(a,mask=None,bins=10):
    if mask is None:
        mask = zeros((a.n,),bool)
        mask[:] = True
    return histogram(norm(a.velo[mask]),bins=bins,normed=True)


def T_of_dist(a):
    return norm2(a.velo)*a.mass[0]/(3.0*BOLTZMANN_K)


def makeplot(fname):

    root, ext = os.path.splitext(fname)

    dist_atoms = Atoms(root+'_dist.xyz')

    fr = frame_reader(fname)

    a = fr.next() # read first frame

    temp_avg = array((temperature(a), temperature(a, a.embed==1)))

    sigma = sqrt(BOLTZMANN_K*300.0/a.mass[0])

    bins = linspace(0.0,5*sigma,200)
    vhist_avg, bins = velo_hist(a, bins=bins)

    qm_bins = linspace(0.0,10*sigma,200)
    qm_vhist_avg, qm_bins = velo_hist(a, mask=a.embed == 1, bins=qm_bins)

    T_of_d_avg = T_of_dist(a)

    temps = []

    dt = 0.1 # ps

    n = 1
    for a in fr:
        temp = array((temperature(a), temperature(a, a.embed==1)))
        temps.append((n*dt,temp[0],temp[1]))

        temp_avg[:] = (n*temp_avg + temp)/float(n+1.0)

        if n%1 == 0:
            print '%d: avg temp %.1f, QM avg temp %.1f' % (n, temp_avg[0], temp_avg[1])

        vhist, bins = velo_hist(a, bins=bins) # reuse previous bins
        vhist_avg[:] = (n*vhist_avg + vhist)/float(n+1.0)

        qm_vhist, qm_bins = velo_hist(a, mask=a.embed == 1, bins=qm_bins)
        qm_vhist_avg[:] = (n*qm_vhist_avg + qm_vhist)/float(n+1.0)

        T_of_d = T_of_dist(a)
        T_of_d_avg[:] = (n*T_of_d_avg + T_of_d)/float(n+1.0)

        n = n + 1

    total_time = n*dt
    avg_time = total_time/10.0

    f = figure(1)
    f.set_size_inches(11.0,8.0)
    clf()
    subplot(211)
    title(root)

    ta = array(temps)
    #plot(ta[:,0],ta[:,1],'b:')
    #plot(ta[:,0],ta[:,2],'r:')

    avg_t = ema(ta[:,1],int(avg_time/dt))
    plot(ta[-len(avg_t):,0],avg_t,'b:')

    avg_qm_t = ema(ta[:,2],int(avg_time/dt))
    plot(ta[-len(avg_qm_t):,0],avg_qm_t,'r:')

    plot(ta[:,0],[temp_avg[0] for i in ta[:,1]],'b-')
    plot(ta[:,0],[temp_avg[1] for i in ta[:,2]],'r-')

    xlabel('Time / ps')
    ylabel('Temp / K')

    subplot(223)

    bar(bins, vhist_avg, width=(bins[1]-bins[0]),fc='b',ec='b')
    vs = arange(bins.min(), bins.max(), 0.5*(bins[1]-bins[0]))
    plot(vs, [max_bolt(a.mass[0],temp_avg[0],v) for v in vs], 'k-', lw=2)
    title(r'Entire system $\bar{T}=%.1f$ K' % temp_avg[0])
    xlabel('Velocity')
    ylabel('Frequency')


    subplot(224)
    bar(qm_bins, qm_vhist_avg, width=(bins[1]-bins[0]),fc='r',ec='r')
    vs = arange(qm_bins.min(), qm_bins.max(), 0.5*(qm_bins[1]-qm_bins[0]))
    plot(vs, [max_bolt(a.mass[0],temp_avg[1],v) for v in vs], 'k-', lw=2)
    title(r'QM region $\bar{T}=%.1f$ K' % temp_avg[1])
    xlabel('Velocity')
    ylabel('Frequency')

    show()

    savefig(root+'.1.eps')

    f = figure(2)
    f.set_size_inches(11.0,8.0)
    clf()
#    h, hx, hy = histogram2d(dist_atoms.dist, T_of_d_avg, (20,20))
#
#    subplot(211)
#    imshow(h,cmap=cm.Greys,extent=[dist_atoms.dist.min(), dist_atoms.dist.max(), \
#                                   T_of_d_avg.min(), T_of_d_avg.max()],aspect='auto')
#    xlabel('Distance from QM region')
#    ylabel('Average temperature')
#
#   subplot(212)
    title(root)
    plot(dist_atoms.dist[a.embed==0], T_of_d_avg[a.embed==0], 'bx')
    plot(dist_atoms.dist[a.embed==1], T_of_d_avg[a.embed==1], 'rx')
    xlabel('Distance from QM region')
    ylabel('Average temperature')

    show()

    savefig(root+'.2.eps')

    

    

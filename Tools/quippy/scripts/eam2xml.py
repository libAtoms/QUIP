"""
Convert EAM tabulated functions to QUIP XML format

James Kermode <james.kermode@kcl.ac.uk> 2014

Input format is that used in Yuri Mishin's .plt files, e.g. those availabe
from http://www.ctcms.nist.gov/potentials/Ni

Usage: eam2xml.py <species1> <species2>... e.g. eam2xml.py Ni Al
"""

import sys
import os
import numpy as np
from quippy import ElementName

species = sys.argv[1:]
atomic_nums = [ElementName.index(sp) for sp in species]

ntypes = len(species)
nspline = 50

print '<EAM_ErcolAd_params n_types="%d" n_spline_V="%d" n_spline_rho="%d" n_spline_F="%d">' % (ntypes, nspline, nspline, nspline)

for ti,(spi,zi) in enumerate(zip(species, atomic_nums)):

    print '<per_type_data atomic_num="%d" type="%d">' % (zi, ti+1)
    
    rho_i = np.loadtxt('f%s.plt' % spi)
    print '<spline_rho>'
    for r, rho in rho_i[::len(rho_i)/nspline]:
        print '<point r="%f" y="%f" b="0.0" c="0.0" d="0.0" />' % (r, rho)
    print '</spline_rho>'


    F_i   = np.loadtxt('F_%s.plt' % spi)
    print '<spline_F>'
    for r, F in F_i[::len(F_i)/nspline]:
        print '<point r="%f" y="%f" b="0.0" c="0.0" d="0.0" />' % (r, F)
    print '</spline_F>'

    print '</per_type_data>'

for ti,(spi,zi) in enumerate(zip(species, atomic_nums)):
    for tj,(spj,zj) in enumerate(zip(species, atomic_nums)):

        if tj > ti:
            continue

        if spi == spj:
            sp = spi
        else:
            sp = spi+spj
            if not os.path.exists('p%s.plt' % sp):
                sp = spj+spi

        phi_ij = np.loadtxt('p%s.plt' % sp)

        rmin = phi_ij[:,0].min()
        rcut = phi_ij[:,0].max()
        print '<per_pair_data atomic_num_i="%d" atomic_num_j="%d" r_min="%f" r_cut="%f">' % (zi, zj, rmin, rcut)
        print '<spline_V>'
        for r, phi in phi_ij[::len(phi_ij)/nspline]:
            print '<point r="%f" y="%f"  b="0.0" c="0.0" d="0.0" />' % (r, phi)
        print '</spline_V>'
        print '</per_pair_data>'

print '</EAM_ErcolAd_params>'

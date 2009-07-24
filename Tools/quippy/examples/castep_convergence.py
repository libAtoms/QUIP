from quippy import *
from pylab import *
import sys, os

assert len(sys.argv[1:]) == 1, 'Usage: castep_convergence.py <castep_log_file>'
castep_file = sys.argv[1]

al = AtomsList(castep_file, save_params=True)

cutoff = [ float(a.cut_off_energy.split()[0]) for a in al ]
energy = [ a.energy for a in al ]
frms   = [ sqrt(a.force.norm2().mean()) for a in al ]
fone   = [ a.force[1,1] for a in al ]

if hasattr(al[0], 'virial'):
   sigma_xx = [ a.virial[1,1] for a in al ]
   sigma_yy = [ a.virial[2,2] for a in al ]
   sigma_zz = [ a.virial[3,3] for a in al ]

figure (1)
title('Energy convergence');
plot(cutoff, energy, 'x-')

figure(2)
title('RMS force convergence');
plot(cutoff, frms, 'x-')

figure(3)
title('Force (1,1) convergence')
plot(cutoff, fone, 'x-')

if hasattr(al[0], 'virial'):
   figure(4)
   clf()
   title(r'Stress convergence')
   plot(cutoff, sigma_xx, 'x-', label=r'$\sigma_{xx}$')
   plot(cutoff, sigma_yy, 'x-', label=r'$\sigma_{yy}$')
   plot(cutoff, sigma_zz, 'x-', label=r'$\sigma_{zz}$')
   legend()

from pylab import *
from quippy import *
from numpy import *
from scipy.special import erf

def pot_erf(at):

    energy = 0.
    sigma = 1.
    force = fzeros((3,at.n))
    
    for i in frange(at.n):
    	for neighb in at.neighbours[i]:

            r_ij = neighb.distance
            u_ij = neighb.cosines
            
            erf_val = erf(r_ij/(sqrt(2.)*sigma))
            erf_deriv = sqrt(2.)/(sqrt(pi)*sigma)*exp(-r_ij*r_ij/(2.0*sigma*sigma))

            energy = energy + 0.5*erf(r_ij/(sqrt(2.)*sigma))
            force[:,i] += erf_deriv*u_ij

    at.params['energy'] = energy
    at.params['virial'] = fzeros((3,3))
    at.add_property('force', force, overwrite=True)
            

p = Potential('CallbackPot')
p.set_callback(pot_erf)

a = Atoms(n=2, lattice=10.*fidentity(3))
a.set_atoms(14)
a.set_cutoff(2.0)

a.pos[:,1] = [0., 0., 0.]

xa = arange(-1, 1, 0.1)
e = []
f = []

for x in xa:
    a.pos[:,2] = [x, 0., 0.]
    a.calc_connect()
    p.calc(a, energy=True, force=True)
    e.append(a.energy)
    f.append(a.force[1,2])

e = array(e)
f = array(f)

clf()
plot(xa, e, 'o-')
plot(xa, f, 'o-')
legend(['energy', 'force'])

p.n_test_gradient(a)

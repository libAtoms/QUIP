import numpy as np
import itertools
import sys
from quippy.atoms import Atoms
from quippy.atomslist import AtomsList
from quippy.potential import Potential
from glob import glob
from pylab import figure, plot, legend, draw, clf, xlabel, ylabel

do_calc = False
do_plot = False

extra_args = {'molpro': {'template': 'dft'}}

codes = ['castep', 'cp2k', 'gpaw', 'molpro', 'vasp']

castep = Potential('FilePot', command='./castep-driver.sh')
cp2k = Potential('FilePot', command='./cp2k-driver.sh')
gpaw = Potential('FilePot', command='./gpaw-driver.sh')
molpro = Potential('FilePot', command='./molpro-driver.sh')
vasp = Potential('FilePot', command='./vasp-driver.sh')

potentials = [ globals()[code] for code in codes ]

def h2_molecule(a, vacuum=10.0):
    h2 = Atoms(n=2, lattice=np.diag([vacuum, vacuum, vacuum]))
    h2.set_atoms([1,1])
    h2.params['bond_length'] = a
    h2.pos[1,1] = -a/2.0
    h2.pos[1,2] = +a/2.0
    return h2

def fit_and_plot(molecules, code, color):
    energy = getattr(molecules, code+'_energy')
    energy = np.array(energy) - min(energy)
    plot(molecules.bond_length, energy, color+'o', label=code.upper()+' data')

    p = np.polyfit(molecules.bond_length, energy, 2)
    bond_length = -p[1]/(2*p[0])
    spring_constant = 2.0*p[0]
    a = np.linspace(min(molecules.bond_length), max(molecules.bond_length), 100)
    plot(a, np.polyval(p, a), color+'-', label=code.upper()+' fit')
    print '|%-10s|%10.3f|%10.1f|' % (code.upper(), bond_length, spring_constant)

def calc_all(molecules, codes):

    for code, pot in zip(codes, potentials):
        args = extra_args.get(code,{})
        for molecule in molecules:
            pot.calc(molecule, energy="%s_energy" % code, **args)

def plot_all(molecules, codes):
    colors = itertools.cycle(['r','g','b','c','m','y'])
    figure(1)
    clf()
    for code, color in zip(codes, colors):
        fit_and_plot(molecules, code, color)
    legend(loc='upper center')
    xlabel(r'Bond length / $\AA{}$')
    ylabel(r'Energy / eV')
    draw()

vacuum = 10.0
bond_lengths = np.linspace(0.7,0.8,5)
molecules = AtomsList([h2_molecule(a,vacuum) for a in bond_lengths])

if do_calc:
    calc_all(molecules, potentials, codes)
    molecules.write('results.xyz')
else:
    molecules = AtomsList('results.xyz')

if do_plot:
    plot_all(molecules)

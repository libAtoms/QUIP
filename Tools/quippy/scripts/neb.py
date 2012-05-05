#!/home/kermode/bin/python

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

from quippy import *
from ase.constraints import FixAtoms
from ase.neb import SingleCalculatorNEB
from ase.optimize import MDMin
from ase.optimize.fire import FIRE
import numpy as np
import sys
import os
import optparse

def write_band(neb, filename):
    out = AtomsWriter(filename)
    for image in neb.images:
        cp = image.copy()
        cp.params['energy'] = image.get_potential_energy()
        cp.add_property('force', image.get_forces().T, overwrite=True)
        out.write(cp)
    out.close()

def callback_calc(at):
    global global_pot, calc_args
    # dummy calculator which calls ForceMixing pot and sets energies to zero
    if at.calc_energy:
        at.params['energy'] = 0.0
    if at.calc_force:
        at.set_cutoff(global_pot.cutoff())
        at.calc_connect()
        print 'Calling global_pot.calc() with args_str=%s' % (calc_args[at.params['image']] + ' force')
        global_pot.calc(at, args_str=calc_args[at.params['image']]+" force")

p = optparse.OptionParser(usage='%prog [OPTIONS] INFILE')
p.add_option('-R', '--refine', action='store', type='int', help='Refine NEB pathway by this many steps')
p.add_option('-I', '--init-args', action='store', help='QUIP Potential init_args string')
p.add_option('-C', '--calc-args', action='store', help='QUIP Potential calc_args string', default='')
p.add_option('-p', '--param-file', action='store', help='Read QUIP potential XML from this file.', default='crack.xml')
p.add_option('-c', '--climb', action='store_true', help='If true, do climbing image NEB', default=False)
p.add_option('-k', action='store', type='float', help='Spring constant for NEB', default=1.0)
p.add_option('-O', '--optimizer', action='store', help='Optimizer to use: currently FIRE or MDMin')
p.add_option('-w', '--write-traj', action='store', type='int', help='Write trajectory every this many steps')
p.add_option('-m', '--minim', action='store_true', help="Minimize start and end points before starting NEB")
p.add_option('--minim-method', action='store', help="Method to use for initial minimisations: sd, cg, cg_n or fire", default="cg_n")
p.add_option('--minim-tol', action='store', type='float', help="Tolerance for initial minimisations - df^2, (eV/A)^2", default=1e-3)
p.add_option('--minim-max-steps', action='store', type='int', help="Maximum number of steps for initial minimisations", default=100)
#p.add_option('--minim-do-pos', action='store', help="If true, relax positions in initial minimisations", default=True)
#p.add_option('--minim-do-lat', action='store_true', help="If true, relax lattice in initial minimisations", default=False)
p.add_option('--verbosity', action='store', help="Set QUIP verbosity")
p.add_option('--run-suffix', action='store', help="run_suffix to be used for each image's calculator. %d will be replaced by number of image")
p.add_option('--callback', action='store_true', help="Use a callback potential")

opt, args = p.parse_args()

print 'init_args', opt.init_args
print 'calc_args', opt.calc_args

if opt.verbosity is not None:
    verbosity_push(verbosity_of_str(opt.verbosity))

try:
    infile, = args
except:
    p.error('exactly one input file is required')
    
basename = os.path.splitext(infile)[0]
images = AtomsList(infile)

if opt.minim:
    pot = Potential(opt.init_args, param_filename=opt.param_file)
    for at in (images[0], images[-1]):
        at.set_cutoff(pot.cutoff() + 1.0)
        pot.minim(at,
                  opt.minim_method,
                  opt.minim_tol,
                  opt.minim_max_steps,
                  do_pos=opt.minim_do_pos,
                  do_lat=opt.minim_do_lat)


images = list(images) # FIXME -- unclear why we have to convert from AtomsList to list

have_constraint = hasattr(images[0], 'move_mask')
if have_constraint:
    constraint = FixAtoms(mask=np.logical_not(images[0].move_mask.view(np.ndarray)))

neb = SingleCalculatorNEB(images, climb=opt.climb, k=opt.k)

if opt.refine:
    neb.refine(opt.refine)

if have_constraint:
    for image in images:
        image.set_constraint(constraint)

if opt.callback:
    global_pot = Potential(opt.init_args, param_filename=opt.param_file)
    global_pot.print_()

calculators = []
calc_args = []
for i, at in enumerate(neb.images):
    if opt.callback:
        p = Potential(callback=callback_calc)
    else:
        p = Potential(opt.init_args, param_filename=opt.param_file)
    calc_args.append(opt.calc_args)
    if opt.run_suffix is not None:
        run_suffix = opt.run_suffix % i
        calc_args[-1] += ' run_suffix=%s' % run_suffix
        at.add_property('hybrid_mark'+run_suffix, at.hybrid_mark)
    at.params['image'] = i
    calculators.append(p)
neb.set_calculators(calculators)

print 'Starting NEB run with %d images' % len(neb.images)

# Optimize:
if opt.optimizer == 'FIRE':
    optimizer = FIRE(neb)
elif opt.optimizer == 'MDMin':
    optimizer = MDMin(neb)
else:
    p.error('Unknown optimizer %s' % opt.optimizer)

if opt.write_traj is not None:
    n = 0
    def write_traj():
        global neb, n
        write_band(neb, '%s-band-%d.xyz' % (basename, n))
        n += 1
    optimizer.attach(write_traj, interval=opt.write_traj)
    
optimizer.run(fmax=0.03)
write_band(neb, '%s-optimized.xyz' % basename)

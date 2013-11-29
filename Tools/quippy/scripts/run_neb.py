#!/usr/bin/env python

from quippy import *
import numpy as np
import sys
import os
import shutil
import optparse

from math import sqrt
from math import atan, pi
from ase.parallel import world, rank, size
from ase.io import read

class SubdirCalculator:

    def __init__(self, pot, chdir, calc_args):
        self.pot = pot
        self.chdir = chdir
        self.calc_args = calc_args

    def get_potential_energy(self, at):
        if self.chdir:
            old_dir = os.getcwd()
            image_dir = os.path.join(old_dir, 'image-%d' % at.image)
            os.chdir(image_dir)

        calc_args = self.calc_args + " energy=energy"
        self.pot.calc(at, args_str=calc_args)

        if self.chdir:
            os.chdir(old_dir)

        return at.energy    

    def get_forces(self, at):
        at.set_cutoff(self.pot.cutoff())
        at.calc_connect()
        if self.chdir:
            old_dir = os.getcwd()
            image_dir = os.path.join(old_dir, 'image-%d' % at.image)
            os.chdir(image_dir)

        calc_args = self.calc_args + " force=force"
        if at.has_property('weight_region1') and not np.all(at.weight_region1 == 0.):
            calc_args = calc_args + " calc_weights=F"
        else:
            calc_args = calc_args + " calc_weights=T"
            
        if '%image' in calc_args:
            calc_args = calc_args.replace('%image', str(at.params['image']))
            
        #print 'rank %d running image %d with args_str=%s' % (rank, at.image, calc_args)
        self.pot.calc(at, args_str=calc_args)
        #at.write('callback_calc_log.xyz', append=True)

        if self.chdir:
            os.chdir(old_dir)

        return at.force.view(np.ndarray).T.copy()

p = optparse.OptionParser(usage='%prog [OPTIONS] INFILE')
p.add_option('-d', '--dry-run', action='store_true', help='Do nothing')
p.add_option('-R', '--refine', action='store', type='int', help='Refine NEB pathway by this many steps')
p.add_option('-I', '--init-args', action='store', help='QUIP Potential init_args string')
p.add_option('-C', '--calc-args', action='store', help='QUIP Potential calc_args string', default='')
p.add_option('-p', '--param-file', action='store', help='Read QUIP potential XML from this file.', default='crack.xml')
p.add_option('--climb', action='store_true', help='If true, do climbing image NEB', default=False)
p.add_option('-k', action='store', type='float', help='Spring constant for NEB', default=1.0)
p.add_option('-O', '--optimizer', action='store', help='Optimizer to use: currently FIRE or MDMin', default='FIRE')
p.add_option('-w', '--write-traj', action='store', type='int', help='Write trajectory every this many steps')
p.add_option('--verbosity', action='store', help="Set QUIP verbosity")
p.add_option('--fmax', action='store', type='float', help="Maximum force for NEB convergence", default=0.03)
p.add_option('--parallel', action='store_true', help="Parallelise over images using MPI")
p.add_option('--chdir', action='store_true', help="Calculate each image in its own subdirectory, named image-i, for i=1..N-1")
p.add_option('--integrate-forces', action='store_true', help="Compute energies by integrating forces along path")
p.add_option('--eval', action='store_true', help="Evaluate forces and energies, including first and last images")
p.add_option('--plot', action='store_true', help="Plot NEB miniumum energy profile")
p.add_option('--plot-color', action='store', help="Color for plot (default black, k)", default="k")
p.add_option('--plot-label', action='store', help="Label for plot")
p.add_option('--nneightol', action='store', type='float', help='Nearest neighbour tolerance', default=1.3)
p.add_option('--exec-code', action='store', help="Python code string to execute on all images before starting minimisation or NEB")
p.add_option('--exec-file', action='store', help="Python code file to execute on all images before starting minimisation or NEB")
p.add_option('--minim', action='store_true', help="Minimise initial and final states before starting NEB")
p.add_option('--minim-method', action='store', help="Method used for minimisation: cg, fire, sd, lbfgs, etc.", default='cg')
p.add_option('--minim-tol', action='store', help="Minim tolerance in (eV/A)**2", default=1e-6, type=float)
p.add_option('--minim-max-steps', action='store', help="Minim max number of steps", default=1000, type=int)
p.add_option('--basename', action='store', help="Stem used to generate output file names")

opt, args = p.parse_args()

if len(args) != 1:
    p.error('exactly 1 input file must be provided')

print 'MPI rank %d, size %d' % (rank, size)

if rank == 0:
    print 'init_args', opt.init_args
    print 'calc_args', opt.calc_args

if opt.verbosity is not None:
    verbosity_push(verbosity_of_str(opt.verbosity))

if opt.basename is not None:
    basename = opt.basename
else:
    basename = os.path.splitext(args[0])[0]
images = AtomsList(args)
images = list(images) # FIXME -- unclear why we have to convert from AtomsList to list

have_constraint = hasattr(images[0], 'move_mask')
if have_constraint:
    from ase.constraints import FixAtoms
    constraint = FixAtoms(mask=np.logical_not(images[0].move_mask.view(np.ndarray)))

if opt.exec_file is not None:
    for at in images:
        execfile(opt.exec_file)

if opt.exec_code is not None:
    for at in images:
        exec(opt.exec_code)

if opt.minim is not None:
    relax_pot = Potential(opt.init_args, param_filename=opt.param_file)

    if len(images) != 2:
        p.error("Number of images should be exactly 2 when --minim option present")
    for at in images:
        at.set_cutoff(relax_pot.cutoff() + 1.0)
        at.calc_connect()
        relax_pot.minim(at, opt.minim_method, opt.minim_tol, opt.minim_max_steps, do_pos=True, do_lat=False, args_str=opt.calc_args)

    del relax_pot

if not opt.dry_run:
    quip_pot = Potential(opt.init_args, param_filename=opt.param_file)
    if rank == 0:
        quip_pot.print_()
    images[0].set_calculator(SubdirCalculator(quip_pot, opt.chdir, opt.calc_args))
    images[-1].set_calculator(SubdirCalculator(quip_pot, opt.chdir, opt.calc_args))

if opt.refine is not None:
   tmp_neb = NEB(images)
   tmp_neb.refine(opt.refine)
   images = tmp_neb.images
   del tmp_neb

if opt.parallel:
    keys = os.environ.keys()
    for key in keys:
        if key.startswith('OMPI_'):
            del os.environ[key]

neb = NEB(images,
          climb=opt.climb,
          k=opt.k,
          parallel=opt.parallel,
          integrate_forces=opt.integrate_forces,
          analyse=opt.dry_run)

if have_constraint:
    for image in images:
        image.set_constraint(constraint)

for image in images:
    image.nneightol = opt.nneightol


for i in range(len(neb.images)):
    at = neb.images[i]
    at.params['image'] = i

    if not opt.dry_run:
        at.set_calculator(SubdirCalculator(quip_pot, opt.chdir, opt.calc_args))

if rank == 0:
    print 'Starting NEB run with %d images' % len(neb.images)
    if not os.path.exists('%s-initial.xyz' % basename):
        neb.write('%s-initial.xyz' % basename)

if opt.eval:
    # actually evaluate end forces as well
    neb.get_forces(all=True)
    neb.write('%s-eval.xyz' % basename)
elif not opt.dry_run:
    # Optimize:
    if opt.optimizer == 'FIRE':
        from ase.optimize.fire import FIRE
        optimizer = FIRE(neb)
    elif opt.optimizer == 'MDMin':
        from ase.optimize import MDMin
        optimizer = MDMin(neb)
    else:
        p.error('Unknown optimizer %s' % opt.optimizer)

    if opt.write_traj is not None:
        def write_traj():
            global neb
            n = 0
            while os.path.exists('%s-band-%d.xyz' % (basename, n)):
                n += 1
            neb.write('%s-band-%d.xyz' % (basename, n))
        optimizer.attach(write_traj, interval=opt.write_traj)

    optimizer.run(fmax=opt.fmax)
    if os.path.exists('%s-optimized.xyz' % basename):
        os.unlink('%s-optimized.xyz' % basename)
    neb.get_forces(all=True)
    neb.write('%s-optimized.xyz' % basename)

if opt.plot:
    neb.plot(color=opt.plot_color, label=opt.plot_label)

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
from ase.neb import NEB
from ase.optimize import MDMin
from ase.optimize.fire import FIRE
from ase.parallel import rank, size
import numpy as np
import sys
import os
import shutil
import optparse

class FracPosNEB(NEB):
    """Subclass of ase.neb.NEB which does interpolation using fractional positions
       to avoid problems with PBC."""
    
    def interpolate(self, initial=0, final=-1):
        """Interpolate linearly between initial and final images."""
        if final < 0:
            final = self.nimages + final
        n = final - initial
        pos1 = self.images[initial].get_scaled_positions()
        pos2 = self.images[final].get_scaled_positions()
        diff = pos2 - pos1
        step = diff / n
        for i in range(1, n):
            self.images[initial + i].set_scaled_positions(self.images[initial + i-1].get_scaled_positions() + step)
            
    def refine(self, steps=1, begin=0, end=-1):
        """Refine the NEB trajectory."""
        if end < 0:
            end = self.nimages + end
        j = begin
        n = end - begin
        for i in range(n):
            for k in range(steps):
                self.images.insert(j + 1, self.images[j].copy())
            self.nimages = len(self.images)
            self.interpolate(j, j + steps + 1)
            j += steps + 1

    def write(self, filename):
        if rank != 0:
            return
        out = AtomsWriter(filename)
        for image in self.images:
            cp = image.copy()
            calc_at = image.get_calculator().atoms
            if calc_at.has_property('force'):
                cp.add_property('force', calc_at.force, overwrite=True)
            out.write(cp)
        out.close()

def callback_calc(at):
    global quip_pot, opt, rank
    keys = os.environ.keys()
    for key in keys:
        if key.startswith('OMPI_'):
            del os.environ[key]
    # dummy calculator which calls ForceMixing pot and sets energies to zero
    if at.calc_energy:
        at.params['energy'] = 0.0
    if at.calc_force:
        at.set_cutoff(quip_pot.cutoff())
        at.calc_connect()
        if opt.chdir:
            old_dir = os.getcwd()
            image_dir = os.path.join(old_dir, 'image-%d' % at.image)
            os.chdir(image_dir)
            print 'rank %d running image %d in dir %s' % (rank, at.image, image_dir)
        calc_args = opt.calc_args + " force"
        if '%image' in calc_args:
            calc_args = calc_args.replace('%image', str(at.params['image']))
        print 'rank %d calling quip_pot.calc() with args_str=%s' % (rank, calc_args)
        quip_pot.calc(at, args_str=calc_args)
        if opt.chdir:
            os.chdir(old_dir)

p = optparse.OptionParser(usage='%prog [OPTIONS] INFILE')
p.add_option('-R', '--refine', action='store', type='int', help='Refine NEB pathway by this many steps')
p.add_option('-I', '--init-args', action='store', help='QUIP Potential init_args string')
p.add_option('-C', '--calc-args', action='store', help='QUIP Potential calc_args string', default='')
p.add_option('-p', '--param-file', action='store', help='Read QUIP potential XML from this file.', default='crack.xml')
p.add_option('--climb', action='store_true', help='If true, do climbing image NEB', default=False)
p.add_option('-k', action='store', type='float', help='Spring constant for NEB', default=1.0)
p.add_option('-O', '--optimizer', action='store', help='Optimizer to use: currently FIRE or MDMin', default='FIRE')
p.add_option('-w', '--write-traj', action='store', type='int', help='Write trajectory every this many steps')
p.add_option('-m', '--minim', action='store_true', help="Minimize start and end points before starting NEB")
p.add_option('--minim-method', action='store', help="Method to use for initial minimisations: sd, cg, cg_n or fire", default="cg_n")
p.add_option('--minim-tol', action='store', type='float', help="Tolerance for initial minimisations - df^2, (eV/A)^2", default=1e-3)
p.add_option('--minim-max-steps', action='store', type='int', help="Maximum number of steps for initial minimisations", default=100)
p.add_option('--verbosity', action='store', help="Set QUIP verbosity")
p.add_option('--fmax', action='store', type='float', help="Maximum force for NEB convergence", default=0.03)
p.add_option('--parallel', action='store_true', help="Parallelise over images using MPI")
p.add_option('--chdir', action='store_true', help="Calculate each image in its own subdirectory, named image-i, for i=1..N-1")

opt, args = p.parse_args()

print 'MPI rank %d, size %d' % (rank, size)

if rank == 0:
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
                  do_pos=True,
                  do_lat=False) # FIXME -- add options for do_pos and do_lat


images = list(images) # FIXME -- unclear why we have to convert from AtomsList to list

have_constraint = hasattr(images[0], 'move_mask')
if have_constraint:
    constraint = FixAtoms(mask=np.logical_not(images[0].move_mask.view(np.ndarray)))

if opt.refine is not None:
   tmp_neb = FracPosNEB(images)
   tmp_neb.refine(opt.refine)
   images = tmp_neb.images
   del tmp_neb

neb = FracPosNEB(images, climb=opt.climb, k=opt.k, parallel=opt.parallel)

if have_constraint:
    for image in images:
        image.set_constraint(constraint)

quip_pot = Potential(opt.init_args, param_filename=opt.param_file)    
if rank == 0:
    quip_pot.print_()

for i, at in enumerate(neb.images):
    pot = Potential(callback=callback_calc, inplace=False)
    at.params['image'] = i
    at.set_calculator(pot)
    at.map_into_cell()

if rank == 0:
    print 'Starting NEB run with %d images' % len(neb.images)

neb.write('%s-initial.xyz' % basename)

# Optimize:
if opt.optimizer == 'FIRE':
    optimizer = FIRE(neb)
elif opt.optimizer == 'MDMin':
    optimizer = MDMin(neb)
else:
    p.error('Unknown optimizer %s' % opt.optimizer)

if rank == 0 and opt.write_traj is not None:
    n = 0
    def write_traj():
        global neb, n
        neb.write('%s-band-%d.xyz' % (basename, n))
        n += 1
    optimizer.attach(write_traj, interval=opt.write_traj)
    
optimizer.run(fmax=opt.fmax)
if rank == 0:
    neb.write('%s-optimized.xyz' % basename)

#!/usr/bin/env python

from quippy import *
import numpy as np
import sys
import os
import shutil
import optparse

#from scipy.integrate import simps

from math import sqrt
from math import atan, pi

from ase.parallel import world, rank, size
from ase.io import read

def calc_df2(i):
    images = AtomsList(["start.xyz", "images-band-%d.*.xyz" % i, "end.xyz"])
    n = NEB(images, k=1.0, integrate_forces=True)
    for i in images:
        i.set_constraint(FixAtoms(mask=np.logical_not(i.arrays['move_mask'])))
        i.set_calculator(StoredValuesCalculator())
    f = n.get_forces()
    return (f**2).sum()

def pbc_diff(a, b):
    """
    Return difference between atomic positions of a and b, correcting for jumps across PBC
    """
    aa = a.copy()
    aa.undo_pbc_jumps() #     a.pos = a.pos
    bb = aa.copy()
    bb.pos[:] = b.pos
    bb.undo_pbc_jumps() #     b.pos = b.pos - (np.dot(b.lattice, np.floor(np.dot(b.g, (b.pos-a.pos))+0.5)))
    return aa.get_positions() - bb.get_positions()

def diff_min_image(a, b):
    diff = a.pos - b.pos + (np.dot(a.lattice, np.floor(np.dot(a.g, (b.pos-a.pos))+0.5)))
    return diff.view(np.ndarray).T


diff_min_image = pbc_diff


class StoredValuesCalculator:
    """Return forces and energies stored in Atoms' arrays dictionary"""

    def get_potential_energy(self, at):
        return float(at.params['energy'])

    def get_forces(self, at):
        f = at.arrays['force'].copy()
        if 'move_mask' in at.arrays:
            f[at.arrays['move_mask'] == 0,:] = 0
        return f


class SubdirCalculator:

    def __init__(self, pot, chdir, calc_args):
        self.pot = pot
        self.chdir = chdir
        self.calc_args = calc_args

    def get_potential_energy(self, at):
        calc_args = self.calc_args + " energy=energy"
        self.pot.calc(at, args_str=calc_args)
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
        


class NEB:
    def __init__(self, images, k=1.0, climb=False, parallel=False, integrate_forces=False, analyse=False):

        if isinstance(images, basestring):
            images = list(AtomsList(images))
        
        self.images = images

        # NB: Need support for variable spring constants
        self.k = k
        self.climb = climb
        self.parallel = parallel
        self.integrate_forces = integrate_forces
        self.natoms = len(images[0])

        # Make sure all the images are of even length.
        # NB: This test should be more elaborate and include species and
        #     possible strangeness in the path.
        assert [len(images[0]) for _ in images] == \
               [len(img) for img in images]

        self.nimages = len(images)
        self.emax = np.nan
        self.imax = None

        # Set the spring force implementation
        self.spring_force = 'norm'

        self.reset()
        if analyse:
            for at in self.images:
                at.set_calculator(StoredValuesCalculator())

            self.get_forces(all=True)
            self.fit()
            for i, at in enumerate(self.images):
                at.add_property('f_real', self.forces['real'][i].T, overwrite=True)
                at.add_property('f_spring', self.forces['spring'][i].T, overwrite=True)
                at.add_property('f_neb', self.forces['neb'][i].T, overwrite=True)
            

    def reset(self):
        # Set up empty arrays to store forces, energies and tangents
        self.forces = {}
        self.forces['real'] = np.zeros((self.nimages, self.natoms, 3))
        self.forces['neb'] = np.zeros((self.nimages, self.natoms, 3))
        self.forces['spring'] = np.zeros((self.nimages, self.natoms, 3))
        self.energies = np.zeros(self.nimages)
        self.tangents = np.zeros((self.nimages, self.natoms, 3))

        # Get the end point energies if they are available.
        if not self.integrate_forces:
            try:
                self.energies[0] = self.images[0].get_potential_energy()
            except:
                self.energies[0] = -np.inf
            try:
                self.energies[-1] = self.images[-1].get_potential_energy()
            except:
                self.energies[-1] = -np.inf
        

    def interpolate(self, initial=0, final=-1):
        """Interpolate linearly between two images. The end images are
           used by default"""
        if final < 0:
            final = self.nimages + final
        n = final - initial
        pos1 = self.images[initial].get_positions()
        d = diff_min_image(self.images[final], self.images[initial]) / n
        for i in range(1, n):
            self.images[initial + i].set_positions(pos1 + i * d)
            

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
        self.reset()            

    def get_positions(self):
        """Return the positions of all the atoms for all the images in
           a single array."""
        positions = np.zeros(((self.nimages - 2) * self.natoms, 3))
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            positions[n1:n2] = image.get_positions()
            n1 = n2
        return positions

    def set_positions(self, positions):
        """Set the positions of the images."""
        n1 = 0
        for image in self.images[1:-1]:
            n2 = n1 + self.natoms
            image.set_positions(positions[n1:n2])
            n1 = n2

            # Parallel NEB with Jacapo needs this:
            try:
                image.get_calculator().set_atoms(image)
            except AttributeError:
                pass

    def update_tangents(self):
        """Update the tangent estimates. Only a forward difference tangent,
           towards the neighboring top energy image, is currently
           supported."""
        images = self.images
        t_m = diff_min_image(images[1], images[0]) # images[1].get_positions() - images[0].get_positions()
        self.tangents[0] = t_m.copy()
        for i in range(1, self.nimages - 1):
            t_p = diff_min_image(images[i+1], images[i])
            #(images[i + 1].get_positions() - images[i].get_positions())
            e = self.energies[i]
            e_m = self.energies[i - 1]
            e_p = self.energies[i + 1]
            if (e < e_m and e > e_p) or \
               (i == self.nimages - 2 and e_p == -np.inf):
                t = t_m.copy()
            elif (e > e_m and e < e_p) or (i == 1 and e_m == -np.inf):
                t = t_p.copy()
            else:
                e_max = max(abs(e_p - e), abs(e_m - e))
                e_min = min(abs(e_p - e), abs(e_m - e))
                if e_p > e_m:
                    t = t_p * e_max + t_m * e_min
                else:
                    t = t_p * e_min + t_m * e_max
                t /= np.vdot(t, t)**0.5
                t *= (np.vdot(t_m, t_m)**0.5 + np.vdot(t_p, t_p)**0.5) / 2.0
            self.tangents[i] = t
            t_m = t_p
        self.tangents[-1] = t_m.copy()

    def calculate_image_forces(self, i):
        """Calculate and store the force for a single image."""
        self.forces['real'][i] = self.images[i].get_forces()

    def calculate_image_energies(self, i):
        if self.integrate_forces:
            # compute energies[i] as -\int{ F \dot dr} along path up to image i
            #de = [0.0]
            #r = [0.0]
            self.energies[i] = 0.0
            for j in range(1, i+1):
                dr = diff_min_image(self.images[j], self.images[j-1])
                f_dr = -np.vdot(self.forces['real'][j], dr)
                self.energies[i] += f_dr
                
                #mod_dr = np.sqrt((dr**2).sum())
                #f_dr = f_dr/mod_dr
                #de.append(f_dr)
                #r.append(r[-1] + mod_dr)

            #self.energies[i] = simps(de, r)

                #self.energies[i] += f_dr
                
                #if j == 1 or j == self.nimages-1:
                #    self.energies[i] += f_dr/3.0
                #elif j % 2 == 1:
                #    self.energies[i] += f_dr*2.0/3.0
                #else:
                #    self.energies[i] += f_dr*4.0/3.0
        else:
            self.energies[i] = self.images[i].get_potential_energy()

    def calculate_energies_and_forces(self, all=False):
        """Calculate and store the forces and energies for the band."""
        images = self.images

        if all:
            self.forces['real'][0] = images[0].get_forces()
            self.forces['real'][-1] = images[-1].get_forces()
            self.energies[0] = images[0].get_potential_energy()
            self.energies[-1] = images[-1].get_potential_energy()

        if not self.parallel:
            # Do all images - one at a time:
            for i in range(1, self.nimages - 1):
                self.calculate_image_forces(i)

            for i in range(1, self.nimages - 1):
                self.calculate_image_energies(i)
        else:
            # Parallelize over images: first the forces ...
            i = rank * (self.nimages - 2) // size + 1
            try:
                self.calculate_image_forces(i)
            except:
                # Make sure other images also fail:
                error = world.sum(1.0)
                raise
            else:
                error = world.sum(0.0)
                if error:
                    raise RuntimeError('Parallel NEB failed')
            for i in range(1, self.nimages - 1):
                root = (i - 1) * size // (self.nimages - 2)
                world.broadcast(self.forces['real'][i], root)

            # ... and now the energies
            try:
                self.calculate_image_energies(i)
            except:
                # Make sure other images also fail:
                error = world.sum(1.0)
                raise
            else:
                error = world.sum(0.0)
                if error:
                    raise RuntimeError('Parallel NEB failed')
            for i in range(1, self.nimages - 1):
                root = (i - 1) * size // (self.nimages - 2)
                world.broadcast(self.energies[i : i + 1], root)

        if self.integrate_forces:
            self.calculate_image_energies(len(self.images)-1)


    def get_forces(self, all=False):
        """Evaluate, modify and return the forces."""

        # Update the real forces and energies
        self.calculate_energies_and_forces(all=all)

        # Update the highest energy image
        self.imax = 1 + np.argsort(self.energies[1:-1])[-1]
        self.emax = self.energies[self.imax]

        # Calculate the tangents of all the images
        self.update_tangents()

        # Prjoect the forces for each image
        self.project_forces()

        return self.forces['neb'][1:self.nimages-1].reshape((-1, 3))


    def get_image_distances(self):
        s = [0]
        for i in range(len(self.images)-1):
            s.append(s[-1] + np.sqrt(((diff_min_image(self.images[i+1], self.images[i])**2).sum())))
        return s

    def get_norm_image_spring_force(self, i):
        """Calculate the 'norm' spring force for a single image."""
        t = self.tangents[i]
        nt = t / np.vdot(t, t)**0.5
        #p_m = self.images[i - 1].get_positions()
        #p = self.images[i].get_positions()
        #p_p = self.images[i + 1].get_positions()

        t_m = diff_min_image(self.images[i], self.images[i-1])  #p  - p_m
        t_p = diff_min_image(self.images[i+1], self.images[i])  #p_p - p
        
        nt_m = np.vdot(t_m, t_m)**0.5
        nt_p = np.vdot(t_p, t_p)**0.5
        return (nt_p - nt_m) * self.k * t

    def get_full_image_spring_force(self, i):
        """Calculate the 'full' spring force for a single image."""
        #p_m = self.images[i - 1].get_positions()
        #p = self.images[i].get_positions()
        #p_p = self.images[i + 1].get_positions()
        t_m = diff_min_image(self.images[i], self.images[i-1])  #p  - p_m
        t_p = diff_min_image(self.images[i+1], self.images[i])  #p_p - p
        return (t_p - t_m) * self.k

    def get_image_spring_force(self, i):
        """Calculate the spring force for a single image."""
        if self.spring_force == 'norm':
            return self.get_norm_image_spring_force(i)
        elif self.spring_force == 'full':
            return self.get_full_image_spring_force(i)
        else:
            e = 'The only supported spring force defintions are: "norm"' + \
                ' and "full".'
            raise NotImplementedError(e)

    def project_forces(self, sort='real'):
        """Project the forces, replace the force components along the path
           with the spring force. The input variable sort is included if
           previous force manipulations have been performed."""
        for i in range(1, self.nimages - 1):
            t = self.tangents[i]
            nt = t / np.vdot(t, t)**0.5
            f_r = self.forces[sort][i]
            f_r_para = np.vdot(f_r, nt) * nt
            f_r_perp = f_r - f_r_para
            if self.climb and i == self.imax:
                self.forces['neb'][i] = f_r - 2 * f_r_para
            else:
                f_s = self.get_image_spring_force(i)
                self.forces['spring'][i] = f_s
                self.forces['neb'][i] = f_r_perp + f_s

    def get_potential_energy(self):
        """Return the energy of the top energy image."""
        return self.emax

    def __len__(self):
        return (self.nimages - 2) * self.natoms

    def write(self, filename):
        for at in self.images:
            if size != 1 and at.image-1 != rank:
                continue
            if size == 1:
                at.write(filename, append=True)
            else:
                filename, ext = os.path.splitext(filename)
                filename = '%s.%02d%s' % (filename, rank, ext)
                at.write(filename)

    def fit(self):
        from ase.neb import fit0
        
        E = self.energies.copy()
        F = self.forces['real'].copy()
        R = [i.get_positions() for i in self.images]
        
        s, E, Sfit, Efit, lines = fit0(E, F, R)
        self.s = s
        self.barrier = max(Efit)
        return s, E, Sfit, Efit, lines

 
    def plot(self, normalize=False, smin=0., smax=1., dE=0., plot_images=True,
             plot_forces=True, plot_ts=True, Efac=1.0, color='k', label=None):

        from pylab import plot, xlim, ylim, xlabel, ylabel, scatter, draw, gca, hlines, subplot, legend, text
        s, E, Sfit, Efit, lines = self.fit()

        if normalize:
            s = np.array(s)
            E = (E + dE)*Efac
            Efit = (Efit + dE)*Efac
            max_s = s.max()
            s = (smax-smin)*s/max_s + smin
            Sfit = (smax-smin)*Sfit/max_s + smin
            lines = [ ((smax-smin)*x/max_s + smin, y*Efac + dE) for (x, y) in lines ]
        else:
            smin = np.min(s)
            smax = np.max(s)

        if plot_images:
            plot(s[1:-1], E[1:-1], color+'x')
        plot (Sfit, Efit, color+'-', label=label)
        if plot_forces:
            for x,y in lines:
                plot(x, y, 'g-')
        plot([s[0]], [E[0]], color+'o')
        plot([s[-1]], [E[-1]], color+'o')
        if plot_ts:
            plot([Sfit[Efit.argmax()]], [Efit.max()], 'ro')
        xlabel('Reaction coordinate $s$')
        ylabel(r'$\Delta E$ / eV')
        xlim(smin-0.05*(smax-smin), smax+0.05*(smax-smin))



if __name__ == '__main__':

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
            relax_pot.minim(at, 'cg', 1e-6, 1000, do_pos=True, do_lat=False, calc_args=opt.calc_args)

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
        neb.write('%s-optimized.xyz' % basename)

    if opt.plot:
        neb.plot(color=opt.plot_color, label=opt.plot_label)

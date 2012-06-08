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

from quippy import available_modules
from quippy.system import print_title
from pylab import plot, xlim, ylim, xlabel, ylabel, scatter, draw, gca, hlines, subplot, legend, text, figure
from quippy.farray import convert_farray_to_ndarray
from quippy.atomslist import AtomsList
import numpy as np
import itertools

__all__ = ['plot', 'plot_energy_error', 'plot_max_force_error',
           'plot_rms_force_error', 'plot_max_stress_error',
           'plot_rms_stress_error', 'scatter_force_error',
           'force_error_statistics', 'plot_force_error']

# Wrap pylab plot() function to automatically convert FortanArray to standard numpy arrays
plot = convert_farray_to_ndarray(plot)

def plot_energy_error(configs, ref_configs, energy_name='energy', energy_ref_name='energy', scale=None, *plot_args, **plot_kwargs):
    p = plot([abs(getattr(at, energy_name) - getattr(ref_at, energy_ref_name))/at.n for (at, ref_at) in zip(configs, ref_configs)], *plot_args, **plot_kwargs)
    xlim(0,len(configs)-1)
    ylabel('Energy error per atom / eV')
    if scale: ylim(0,scale)
    return p

def plot_max_force_error(configs, ref_configs, force_name='force', force_ref_name='force', scale=None, *plot_args, **plot_kwargs):
    p = plot([abs(getattr(at, force_name) - getattr(ref_at, force_ref_name)).max() for (at, ref_at) in zip(configs, ref_configs)], *plot_args, **plot_kwargs)
    xlim(0,len(configs)-1)
    ylabel('Max force error / eV/A')
    if scale: ylim(0,scale)
    return p

def plot_rms_force_error(configs, ref_configs, force_name='force', force_ref_name='force', scale=None, *plot_args, **plot_kwargs):
    p = plot([np.sqrt(((getattr(at, force_name) - getattr(ref_at, force_ref_name)).reshape(3*at.n)**2).mean())
          for (at, ref_at) in zip(configs, ref_configs)], *plot_args, **plot_kwargs)
    xlim(0,len(configs)-1)
    ylabel('RMS force error / eV/A')
    if scale: ylim(0,scale)
    return p

def plot_max_stress_error(configs, ref_configs, virial_name='virial', virial_ref_name='virial', scale=None, *plot_args, **plot_kwargs):
    max_stress_error = []
    for (at, ref_at) in zip(configs, ref_configs):
        if not hasattr(ref_at, virial_name):
            max_stress_error.append(0)
        else:
            max_stress_error.append((getattr(at, virial_name) - getattr(ref_at, virial_ref_name)).max()/at.cell_volume()*GPA)

    p = plot(max_stress_error, *plot_args, **plot_kwargs)
    xlim(0,len(configs)-1)
    ylabel('Max stress error / GPa')
    if scale: ylim(0,scale)
    return p

def plot_rms_stress_error(configs, ref_configs, virial_name='virial', virial_ref_name='virial', scale=None, *plot_args, **plot_kwargs):
    rms_stress_error = []
    for (at, ref_at) in zip(configs, ref_configs):
        if not hasattr(ref_at, virial_name):
            rms_stress_error.append(0)
        else:
            rms_stress_error.append(np.sqrt(((getattr(at, virial_name) - getattr(ref_at, virial_ref_name)).reshape(9)**2).mean())/at.cell_volume()*GPA)

    p = plot(rms_stress_error, *plot_args, **plot_kwargs)
    xlim(0,len(configs)-1)
    ylabel('RMS stress error / GPa')
    if scale: ylim(0,scale)
    return p

def scatter_force_error(configs, ref_configs, force_name='force', force_ref_name='force', *plot_args, **plot_kwargs):
    ref_force = np.hstack(getattr(ref_configs, force_ref_name))
    ref_force = ref_force.reshape(ref_force.size, order='F')

    force = np.hstack(getattr(configs, force_name))
    force = force.reshape(force.size, order='F')

    s = scatter(abs(ref_force), abs(ref_force - force), *plot_args, **plot_kwargs)
    xlim(0, abs(ref_force).max())
    ylim(0, abs(ref_force - force).max())
    xlabel('Reference forces / eV/A')
    ylabel('Force error / eV/A')

    rms_error = (((force - ref_force)**2).mean())**0.5
    hlines(rms_error, 0, abs(ref_force).max(), lw=2, color='r')

    return s


def force_error_statistics(configs, ref_configs, force_name='force', force_ref_name='force'):
    ref_force = np.hstack(getattr(ref_configs, force_ref_name))
    ref_force = ref_force.reshape(ref_force.size, order='F')

    force = np.hstack(getattr(configs, force_name))
    force = force.reshape(force.size, order='F')

    names = ['default']
    if hasattr(ref_configs, 'config_type'):
        ct = list(np.hstack([[at.config_type]*3*at.n for at in ref_configs]))
        names.extend(sorted(set(ct)))

    res = {}
    for name in names:

        if name == 'default':
            start = 0
            stop  = len(force)
        else:
            start = ct.index(name)
            stop  = len(ct)-ct[::-1].index(name)

        print_title('config_type "%s" indices [%d:%d]' % (name, start, stop))

        this_force = force[start:stop]
        this_ref_force = ref_force[start:stop]

        max_error = abs(this_force - this_ref_force).max()
        rms_error = (((this_force - this_ref_force)**2).mean())**0.5
        rm4_error = (((this_force - this_ref_force)**4).mean())**0.25
        rm8_error = (((this_force - this_ref_force)**8).mean())**0.125

        print 'Max force error %.3f eV/A' % max_error
        print 'RMS force error %.3f eV/A' % rms_error
        print 'RM4 force error %.3f eV/A' % rm4_error
        print 'RM8 force error %.3f eV/A' % rm8_error
        print

        res[name] = (max_error, rms_error, rm4_error, rm8_error)

    return res


def plot_force_error(configs, ref_configs, force_name='force', force_ref_name='force', *plot_args, **plot_kwargs):
    ref_force = np.hstack(getattr(ref_configs, force_ref_name))
    ref_force = ref_force.reshape(ref_force.size, order='F')

    force = np.hstack(getattr(configs, force_name))
    force = force.reshape(force.size, order='F')

    plot(abs(force - ref_force), *plot_args, **plot_kwargs)
    xlim(0, len(force))

    if hasattr(ref_configs, 'config_type'):
        ct = list(np.hstack([[at.config_type]*3*at.n for at in ref_configs]))
        label_axes_with_config_types(ct)

        names = sorted(set(ct))
        for name in names:
            start = ct.index(name)
            stop  = len(ct)-ct[::-1].index(name)

            this_force = force[start:stop]
            this_ref_force = ref_force[start:stop]

            rms_error = ((this_force - this_ref_force)**2).mean()**0.5
            hlines(rms_error, start, stop, lw=3, color='r')


def label_axes_with_config_types(config_types):
    names = sorted(set(config_types))
    major_tics = []
    minor_tics = []

    idx1 = 0
    for name1, name2 in zip(names, names[1:] + [None]):
        idx1 = config_types.index(name1)
        if name2 is not None:
            idx2 = config_types.index(name2)
        else:
            idx2 = len(config_types)
        half = (idx2+idx1)/2
        minor_tics.append(half)
        major_tics.append(idx2)

    import matplotlib.ticker as ticker
    ax = gca()
    ax.xaxis.set_major_locator(ticker.FixedLocator(major_tics))
    ax.xaxis.set_minor_locator(ticker.FixedLocator(minor_tics))

    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(names))

    for tick in ax.xaxis.get_minor_ticks():
        tick.tick1line.set_markersize(0)
        tick.tick2line.set_markersize(0)
        tick.label1.set_horizontalalignment('center')
        tick.label1.set_rotation(90)

    draw()


if 'ase' in available_modules:

    from ase.neb import fit
    from ase.constraints import FixAtoms
    from ase.calculators.singlepoint import SinglePointCalculator

    __all__.extend(['neb_plot', 'neb_plot_multiple'])

    def neb_plot(images, normalize=False, smin=0., smax=1., dE=0., plot_images=True,
                 plot_forces=True, plot_ts=True, Efac=1.0, color='k', label=None):
        s, E, Sfit, Efit, lines = fit(images)

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
        return s, E, Sfit, Efit, lines


    def neb_plot_multiple(filenames=None, images_list=None, offset=0, dE=0.0, Efac=1.0, cascade=False,
                          plot_barrier_E=False, barrier_xfunc=None, barrier_xlabel='', barrier_polydegree=1,
                          label_func=None, pathway_fig=None, barrier_fig=None, barrier_Escale=None, barrier_label=None,
                          barrier_color='k',**plot_args):
        barrier_E = []
        atoms = []

        if (filenames is None) + (images_list is None) != 1:
            raise ValueError('either filenames or images must be present, but not both')

        if filenames is not None:
            images_list = [ AtomsList(filename) for filename in filenames ]

        if plot_barrier_E:
            if barrier_xfunc is None:
                raise ValueError('xfunc must be present to plot barrier_E')

        if pathway_fig is not None:
            figure(pathway_fig)
        elif plot_barrier_E:
            subplot(121)

        colors = itertools.cycle(['r','g','b','c','m','y','k'])
        for i, (images, color) in enumerate(zip(images_list,colors)):
            initial = images[0]
            constraint = FixAtoms(mask=np.logical_not(initial.move_mask.view(np.ndarray)))

            for image in images:
                energy = image.energy
                forces = image.force.view(np.ndarray).T
                image.set_calculator(SinglePointCalculator(energy, forces, None, None, image))
                image.set_constraint(constraint)

            smin, smax = offset, 1. + offset
            if cascade:
                smin += i
                smax += i

            label = None
            if label_func is not None:
                label = label_func(images[0])
            s, E, Sfit, Efit, lines = neb_plot(images,
                                               True,
                                               smin=smin,
                                               smax=smax,
                                               dE=dE,
                                               Efac=Efac,
                                               color=color,
                                               label=label,
                                               **plot_args)
            if cascade:
                dE = E[-1]
            draw()
            barrier_E.append(Efit.max())

        if label_func is not None:
            legend(loc='lower left')

        if plot_barrier_E:
            if barrier_fig is not None:
                figure(barrier_fig)
            else:
                subplot(122)
            barrier_x = np.array([ barrier_xfunc(images[0]) for images in images_list ])
            barrier_E = np.array(barrier_E)

            if barrier_Escale is not None:
                barrier_E = barrier_E*barrier_Escale
            
            plot(barrier_x, barrier_E, barrier_color+'o', mec=barrier_color)

            p = np.polyfit(barrier_x, barrier_E, barrier_polydegree)
            r = np.roots(p)
            print 'Polynomial fit of degree', barrier_polydegree, p
            print 'Roots', r

            r_0 = r[r > barrier_x.max()].min() # first root after data
            
            xFit = np.linspace(barrier_x.min(), r_0)
            Efit = np.polyval(p, xFit)
            plot(xFit, Efit, barrier_color+'-', label=barrier_label)

            ylim(-barrier_E.max()*0.1)

            name, unit = barrier_xlabel.split('/',1)
            name = name.strip()
            unit = unit.strip()
            #text(r_0, 
            #     barrier_E.max()*0.05,
            #     r'%s$^{(0)}$ = %.2f % s' % (name, r_0, unit),
            #     horizontalalignment='left')

            xlabel(barrier_xlabel)
            ylabel(r'$\Delta E_\mathrm{act}$ / eV')
            x_min, x_max = xlim()
            hlines(0., x_min, x_max, linestyle='dashed')
            xlim(x_min, x_max)

        return images_list

#!/usr/bin/env python
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

from quippy import available_modules, print_title
from pylab import plot, xlim, ylim, xlabel, ylabel, scatter, draw, gca, hlines
from quippy.farray import convert_farray_to_ndarray
import numpy

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
    p = plot([numpy.sqrt(((getattr(at, force_name) - getattr(ref_at, force_ref_name)).reshape(3*at.n)**2).mean())
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
            rms_stress_error.append(numpy.sqrt(((getattr(at, virial_name) - getattr(ref_at, virial_ref_name)).reshape(9)**2).mean())/at.cell_volume()*GPA)
            
    p = plot(rms_stress_error, *plot_args, **plot_kwargs)
    xlim(0,len(configs)-1)
    ylabel('RMS stress error / GPa')
    if scale: ylim(0,scale)
    return p

def scatter_force_error(configs, ref_configs, force_name='force', force_ref_name='force', *plot_args, **plot_kwargs):
    ref_force = numpy.hstack(getattr(ref_configs, force_ref_name))
    ref_force = ref_force.reshape(ref_force.size, order='F')

    force = numpy.hstack(getattr(configs, force_name))
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
    ref_force = numpy.hstack(getattr(ref_configs, force_ref_name))
    ref_force = ref_force.reshape(ref_force.size, order='F')

    force = numpy.hstack(getattr(configs, force_name))
    force = force.reshape(force.size, order='F')

    names = ['default']
    if hasattr(ref_configs, 'config_type'):
        ct = list(numpy.hstack([[at.config_type]*3*at.n for at in ref_configs]))
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
    ref_force = numpy.hstack(getattr(ref_configs, force_ref_name))
    ref_force = ref_force.reshape(ref_force.size, order='F')

    force = numpy.hstack(getattr(configs, force_name))
    force = force.reshape(force.size, order='F')

    plot(abs(force - ref_force), *plot_args, **plot_kwargs)
    xlim(0, len(force))

    if hasattr(ref_configs, 'config_type'):
        ct = list(numpy.hstack([[at.config_type]*3*at.n for at in ref_configs]))
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

        
        

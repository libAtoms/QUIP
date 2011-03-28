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

from quippy import available_modules
from pylab import plot, xlim, ylim, xlabel, ylabel, scatter
from quippy.farray import convert_farray_to_ndarray
import numpy

# Wrap pylab plot() function to automatically convert FortanArray to standard numpy arrays
plot = convert_farray_to_ndarray(plot)

def plot_energy_error(configs, ref_configs, energy_name='energy', energy_ref_name='energy', scale=None, **plot_args):
    p = plot([abs(getattr(at, energy_name) - getattr(ref_at, energy_ref_name))/at.n for (at, ref_at) in zip(configs, ref_configs)], **plot_args)
    xlim(0,len(configs)-1)
    ylabel('Energy error per atom / eV')
    if scale: ylim(0,scale)
    return p

def plot_max_force_error(configs, ref_configs, force_name='force', force_ref_name='force', scale=None, **plot_args):
    p = plot([abs(getattr(at, force_name) - getattr(ref_at, force_ref_name)).max() for (at, ref_at) in zip(configs, ref_configs)], **plot_args)
    xlim(0,len(configs)-1)
    ylabel('Max force error / eV/A')
    if scale: ylim(0,scale)
    return p

def plot_rms_force_error(configs, ref_configs, force_name='force', force_ref_name='force', scale=None, **plot_args):
    p = plot([numpy.sqrt(((getattr(at, force_name) - getattr(ref_at, force_ref_name)).reshape(3*at.n)**2).mean())
          for (at, ref_at) in zip(configs, ref_configs)], **plot_args)
    xlim(0,len(configs)-1)
    ylabel('RMS force error / eV/A')
    if scale: ylim(0,scale)
    return p

def plot_max_stress_error(configs, ref_configs, virial_name='virial', virial_ref_name='virial', scale=None, **plot_args):
    max_stress_error = []
    for (at, ref_at) in zip(configs, ref_configs):
        if not hasattr(ref_at, opt.virial_name):
            max_stress_error.append(0)
        else:
            max_stress_error.append((getattr(at, opt.virial_name) - getattr(ref_at, opt.virial_ref_name)).max()/at.cell_volume()*GPA)
            
    p = plot(max_stress_error, **plot_args)
    xlim(0,len(configs)-1)
    ylabel('Max stress error / eV')
    if scale: ylim(0,scale)
    return p

def plot_rms_stress_error(configs, ref_configs, virial_name='virial', virial_ref_name='virial', scale=None, **plot_args):
    rms_stress_error = []
    for (at, ref_at) in zip(configs, ref_configs):
        if not hasattr(ref_at, opt.virial_name):
            rms_stress_error.append(0)
        else:
            rms_stress_error.append(numpy.sqrt(((getattr(at, opt.virial_name) - getattr(ref_at, opt.virial_ref_name)).reshape(9)**2).mean())/at.cell_volume()*GPA)
            
    p = plot(rms_stress_error, **plot_args)
    xlim(0,len(configs)-1)
    ylabel('RMS stress error / eV')
    if scale: ylim(0,scale)
    return p

def scatter_force_error(configs, ref_configs, force_name='force', force_ref_name='force', **plot_args):
    ref_force = getattr(ref_configs, force_ref_name)
    ref_force.reshape(ref_force.size)

    force = getattr(configs, force_name)
    force.reshape(force.size)
                            
    s = scatter(abs(ref_force), abs(ref_force - force), **plot_args)
    xlim(0, abs(ref_force).max())
    ylim(0, abs(ref_force - force).max())
    xlabel('Reference forces / eV/A')
    ylabel('Force error / eV/A')
    return s


def force_error_statistics(configs, ref_configs, force_name='force', force_ref_name='force'):
    ref_force = getattr(ref_configs, force_ref_name)
    ref_force.reshape(ref_force.size)

    force = getattr(configs, force_name)
    force.reshape(force.size)
    
    max_error = abs(force - ref_force).max()
    rms_error = (((force - ref_force)**2).mean())**0.5
    rm4_error = (((force - ref_force)**4).mean())**0.25
    rm8_error = (((force - ref_force)**8).mean())**0.125

    print 'Max force error %.3f eV/A' % max_error
    print 'RMS force error %.3f eV/A' % rms_error
    print 'RM4 force error %.3f eV/A' % rm4_error
    print 'RM8 force error %.3f eV/A' % rm8_error

    return (max_error, rms_error, rm4_error, rm8_error)
    

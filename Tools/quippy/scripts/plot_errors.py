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

from pylab import *
from quippy import *
import optparse

p = optparse.OptionParser(usage='%prog [options] <input files>...')

p.add_option('-r', '--ref_file', action='store', help="""Separate XYZ file containing reference configs (i.e. teaching data)""")
p.add_option('-E', '--energy', action='store_true', help="""Plot energy errors.""", default=False)
p.add_option('-F', '--max-force', action='store_true', help="""Plot maximum force errors.""", default=False)
p.add_option('-f', '--rms-force', action='store_true', help="""Plot RMS force errors.""", default=False)
p.add_option('-S', '--max-stress', action='store_true', help="""Plot maximum stress errors.""", default=False)
p.add_option('-s', '--rms-stress', action='store_true', help="""Plot RMS stress errors.""", default=False)
p.add_option('-m', '--sparse-points', action='store_true', help="""Plot number of sparse points""", default=False)
p.add_option('-p', '--prefix', action='store', help="""Prefix for force, energy and virial property names""")
p.add_option('-P', '--ref-prefix', action='store', help="""Prefix for reference force, energy and virial property names""")
p.add_option('--force-name', action='store', help="""Name of property containing forces (default "force")""", default="force")
p.add_option('--energy-name', action='store', help="""Name of parameter containing energies (default "energy")""", default="energy")
p.add_option('--virial-name', action='store', help="""Name of parameter containing virials (default "virial")""", default="virial")
p.add_option('--force-ref-name', action='store', help="""Name of property containing reference forces (default "force")""", default="force")
p.add_option('--energy-ref-name', action='store', help="""Name of parameter containing reference energies (default "energy")""", default="energy")
p.add_option('--virial-ref-name', action='store', help="""Name of parameter containing reference virials (default "virial")""", default="virial")
p.add_option('--energy-scale', action='store', type='float')
p.add_option('--max-force-scale', action='store', type='float')
p.add_option('--rms-force-scale', action='store', type='float')
p.add_option('--max-stress-scale', action='store', type='float')
p.add_option('--rms-stress-scale', action='store', type='float')
p.add_option('-L', '--legend', action='store', help="Text to use for figure legend.")
p.add_option('-o', '--output', action='store', help="""Output file name for plot.""")

opt, args = p.parse_args()

if len(args) < 1:
    p.error('At least one input file expected.')

if opt.ref_file is not None:
    ref_configs = AtomsList(opt.ref_file)
    ref_configs.loadall()

if opt.prefix is not None:
    opt.force_name = opt.prefix + opt.force_name
    opt.energy_name = opt.prefix + opt.energy_name
    opt.virial_name = opt.prefix + opt.virial_name

if opt.ref_prefix is not None:
    opt.force_ref_name = opt.ref_prefix + opt.force_ref_name
    opt.energy_ref_name = opt.ref_prefix + opt.energy_ref_name
    opt.virial_ref_name = opt.ref_prefix + opt.virial_ref_name

n_plots = opt.energy + opt.max_force + opt.rms_force + opt.max_stress + opt.rms_stress + opt.sparse_points

if n_plots == 0:
    p.error('Nothing to do - at least one of --energy, --max-force or --rms-force required')

figure(figsize=(8,10))

for pot_file in args:
    configs = AtomsList(pot_file)
    configs.loadall()

    if opt.ref_file is not None:
        if len(ref_configs) != len(configs):
            p.error('Number of configs in file %s (%d) != number of reference configs (%d)' % (pot_file, len(configs), len(ref_configs)))
    else:
        ref_configs = configs # no separate reference file, so we alias ref_configs to configs
    
    i = 1
    if opt.energy:
        subplot(n_plots, 1, i)
        plot([abs(getattr(at, opt.energy_name) - getattr(ref_at, opt.energy_ref_name))/at.n for (at, ref_at) in zip(configs, ref_configs)])
        xlim(0,len(configs)-1)
        ylabel('Energy error per atom / eV')
        if opt.energy_scale: ylim(0,opt.energy_scale)
        i += 1

    if opt.max_force:
        subplot(n_plots, 1, i)
        plot([abs(getattr(at, opt.force_name) - getattr(ref_at, opt.force_ref_name)).max() for (at, ref_at) in zip(configs, ref_configs)])
        xlim(0,len(configs)-1)
        ylabel('Max force error / eV/A')
        if opt.max_force_scale: ylim(0,opt.max_force_scale)        
        i += 1

    if opt.rms_force:
        subplot(n_plots, 1, i)
        plot([sqrt(((getattr(at, opt.force_name) - getattr(ref_at, opt.force_ref_name)).reshape(3*at.n)**2).mean()) for (at, ref_at) in zip(configs, ref_configs)])
        xlim(0,len(configs)-1)
        ylabel('RMS force error / eV/A')
        if opt.rms_force_scale: ylim(0,opt.rms_force_scale)        
        i += 1

    if opt.max_stress or opt.rms_stress:
        max_stress_error = []
        rms_stress_error = []
        for (at, ref_at) in zip(configs, ref_configs):
            if not hasattr(ref_at, opt.virial_name):
                max_stress_error.append(0)
                rms_stress_error.append(0)
            else:
                max_stress_error.append((getattr(at, opt.virial_name) - getattr(ref_at, opt.virial_ref_name)).max()/at.cell_volume()*GPA)
                rms_stress_error.append(sqrt(((getattr(at, opt.virial_name) - getattr(ref_at, opt.virial_ref_name)).reshape(9)**2).mean())/at.cell_volume()*GPA)

                
        if opt.max_stress:
            subplot(n_plots, 1, i)            
            plot(max_stress_error)
            xlim(0,len(configs)-1)
            ylabel('Max stress error / eV')
            if opt.max_stress_scale: ylim(0,opt.max_stress_scale)                    
            i += 1

        if opt.rms_stress:
            subplot(n_plots, 1, i)            
            plot(rms_stress_error)
            xlim(0,len(configs)-1)
            ylabel('RMS stress error / eV')
            if opt.rms_stress_scale: ylim(0,opt.rms_stress_scale)                    
            i += 1

    if opt.sparse_points:
        subplot(n_plots, 1, i)
        n_sparse = float(sum([at.sparse.count() for at in configs]))
        plot([at.sparse.count()/float(at.n) for at in configs])
        xlim(0,len(configs)-1)
        ylabel('N_sparse / N_at')
    

xlabel('Configuration number')

if opt.legend is not None:
    subplots_adjust(bottom=0.15)
    labels = opt.legend.split()
    figlegend(gca().lines, labels, 'lower center', ncol=len(labels))
        
if opt.output is not None:
    savefig(opt.output)
else:
    if not rcParams['interactive']: show()

        

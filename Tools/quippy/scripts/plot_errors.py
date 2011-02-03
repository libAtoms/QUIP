from pylab import *
from quippy import *
import optparse

p = optparse.OptionParser(usage='%prog [options] <ref file> <pot file>...')

p.add_option('-E', '--energy', action='store_true', help="""Plot energy errors.""", default=False)
p.add_option('-F', '--max-force', action='store_true', help="""Plot maximum force errors.""", default=False)
p.add_option('-f', '--rms-force', action='store_true', help="""Plot RMS force errors.""", default=False)
p.add_option('-S', '--max-stress', action='store_true', help="""Plot maximum stress errors.""", default=False)
p.add_option('-s', '--rms-stress', action='store_true', help="""Plot RMS stress errors.""", default=False)
p.add_option('--force-name', action='store', help="""Name of property containing forces (default "force")""", default="force")
p.add_option('--energy-name', action='store', help="""Name of parameter containing energies (default "energy")""", default="energy")
p.add_option('--virial-name', action='store', help="""Name of parameter containing virials (default "virial")""", default="virial")
p.add_option('--energy-scale', action='store', type='float')
p.add_option('--max-force-scale', action='store', type='float')
p.add_option('--rms-force-scale', action='store', type='float')
p.add_option('--max-stress-scale', action='store', type='float')
p.add_option('--rms-stress-scale', action='store', type='float')
p.add_option('-l', '--legend', action='store_true', help="If given, draw a legend at bottom of plot")
p.add_option('-o', '--output', action='store', help="""Output file name for plot.""")

opt, args = p.parse_args()

if len(args) < 2:
    p.error('At least two input files expected.')

ref_file, pot_files = args[0], args[1:]

ref_configs = AtomsList(ref_file)
ref_configs.loadall()

n_plots = opt.energy + opt.max_force + opt.rms_force + opt.max_stress + opt.rms_stress

if n_plots == 0:
    p.error('Nothing to do - at least one of --energy, --max-force or --rms-force required')

figure(figsize=(8,10))

for pot_file in pot_files:
    configs = AtomsList(pot_file)
    configs.loadall()

    if len(ref_configs) != len(configs):
        p.error('Number of configs in file %s (%d) != number of reference configs (%d)' % (pot_file, len(configs), len(ref_configs)))

    i = 1
    if opt.energy:
        subplot(n_plots, 1, i)
        plot([abs(getattr(at, opt.energy_name) - getattr(ref_at, opt.energy_name)) for (at, ref_at) in zip(configs, ref_configs)])
        ylabel('Energy error / eV')
        if opt.energy_scale: ylim(0,opt.energy_scale)
        i += 1

    if opt.max_force:
        subplot(n_plots, 1, i)
        plot([abs(getattr(at, opt.force_name) - getattr(ref_at, opt.force_name)).max() for (at, ref_at) in zip(configs, ref_configs)])
        ylabel('Max force error / eV/A')
        if opt.max_force_scale: ylim(0,opt.max_force_scale)        
        i += 1

    if opt.rms_force:
        subplot(n_plots, 1, i)
        plot([sqrt(((getattr(at, opt.force_name) - getattr(ref_at, opt.force_name)).reshape(3*at.n)**2).mean()) for (at, ref_at) in zip(configs, ref_configs)])
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
                max_stress_error.append((getattr(at, opt.virial_name) - getattr(ref_at, opt.virial_name)).max()/at.cell_volume()*GPA)
                rms_stress_error.append(sqrt(((getattr(at, opt.virial_name) - getattr(ref_at, opt.virial_name)).reshape(9)**2).mean())/at.cell_volume()*GPA)

                
        if opt.max_stress:
            subplot(n_plots, 1, i)            
            plot(max_stress_error)
            ylabel('Max stress error / eV')
            if opt.max_stress_scale: ylim(0,opt.max_stress_scale)                    
            i += 1

        if opt.rms_stress:
            subplot(n_plots, 1, i)            
            plot(rms_stress_error)
            ylabel('RMS stress error / eV')
            if opt.rms_stress_scale: ylim(0,opt.rms_stress_scale)                    
            i += 1

xlabel('Configuration number')

if opt.legend:
    subplots_adjust(bottom=0.2)
    legend(pot_files, bbox_to_anchor=(0.85, 0.15), bbox_transform=gcf().transFigure)
        
if opt.output is not None:
    savefig(opt.output)
else:
    if not rcParams['interactive']: show()

        

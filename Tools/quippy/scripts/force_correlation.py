#!/usr/bin/env python

from pylab import *
from quippy import *
import optparse

p = optparse.OptionParser(usage='%prog [options] <input file>')

p.add_option('-c', '--calc', action='store_true', help="""If true, calculate forces for each config using potential.""")
p.add_option('-i', '--init-args', action='store', help="""Initialisation arguments for Potential.""")
p.add_option('-p', '--param-file', action='store', help="""XML parameter file.""")
p.add_option('-D', '--dft-force', action='store', help="""Name of property containing reference forces (default "force")""", default="force")
p.add_option('-P', '--pot-force', action='store', help="""Name of property containing potential forces (default "pot_force")""", default="pot_force")
p.add_option('-a', '--atomic-units', action='store_true', help="""Convert reference forces from atomic units (Hartree/Bohr).""")
p.add_option('-o', '--output', action='store', help="""Output file name for plot.""")
p.add_option('-y', '--ymax', action='store', type='float', help="""Maximum limit for y-axis""")
p.add_option('-t', '--title', action='append', help="""Title for plot. If multiple input files, option should be given multiple times.""")
p.add_option('-H', '--strip-hydrogen', action='store_true', help="""Remove hydrogen atoms from input files""")
p.add_option('-m', '--mask', action='store', help="""Only include atoms which sastisfy mask, e.g. hybrid_mark == HYBRID_ACTIVE_MARK""")

opt, args = p.parse_args()

if len(args) < 1:
    p.error('At least one input file expected.')

if opt.calc:
    if opt.param_file is not None:
        param_str = param_str=open(opt.param_file).read()
    if opt.init_args is None:
        if opt.param_file is not None:
            import xml.dom.minidom
            xml = xml.dom.minidom.parse(opt.param_file)
            opt.init_args = 'xml_label=%s' % str(xml.documentElement.tagName)
        else:
            p.error('--init-args (-i) argument mandatory when --param-file (-p) not present')
    else:
        param_str = quip_xml_parameters(opt.init_args)
    p = Potential(opt.init_args, param_str=param_str)

clf()
for i,arg in enumerate(args):
    subplot(1,len(args),i+1)
    dft_forces = []
    force_errors = []

    for at in AtomsList(arg):
        if opt.calc:
            if opt.strip_hydrogen:
                at = at.select(at.z != 1)
            at.set_cutoff(p.cutoff())
            at.calc_connect()
            p.calc(at, args_str="force=%s" % opt.pot_force)            

        dft_force = getattr(at, opt.dft_force)

        if opt.atomic_units:
            dft_force = dft_force*(HARTREE/BOHR)

        pot_force = getattr(at, opt.pot_force)

        if opt.mask is not None:
            dft_force = dft_force[:,eval(opt.mask)]
            pot_force = pot_force[:,eval(opt.mask)]
            
        force_error = abs(dft_force - pot_force)
        dft_forces.extend(list(dft_force.reshape([3*dft_force.shape[1]])))
        force_errors.extend(list(force_error.reshape([3*dft_force.shape[1]])))


    dft_forces = array(dft_forces)
    force_errors = array(force_errors)

    scatter(abs(dft_forces), force_errors, s=1)
    xlim(0, max(dft_forces))
    if opt.ymax is None:
        ylim(0, max(force_errors))
    else:
        ylim(0, opt.ymax)
    if i == 0:
        ylabel('Force error / eV/A')
    else:
        first_yticks = gcf().axes[0].get_yticks()
        yticks(first_yticks, ['']*len(first_yticks))
    if opt.title is not None:
        title(opt.title[i])

subplots_adjust(bottom=0.15)
text(0.5, 0.075, 'DFT force / eV/A', horizontalalignment='center',
     verticalalignment='center', transform=gcf().transFigure)

if opt.output is not None:
    savefig(opt.output)
else:
    if not rcParams['interactive']: show()

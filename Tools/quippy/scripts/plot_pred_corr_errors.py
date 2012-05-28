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

import os
import sys
import itertools
import matplotlib.pyplot as plt
import numpy as np
import optparse
from quippy.util import is_interactive_shell
from matplotlib import rc

p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-l', '--label', action='append', help='Labels for figure legend. One per input filename')
p.add_option('-o', '--output', action='store', help='Output file', default='pred-corr-errors.pdf')
p.add_option('-D', '--no-display', action='store_false', help='Do not display figure with "open" command')
p.add_option('-t', '--title', action='store', help='Title for plot', default='Predictor-corrector Force Errors')
p.add_option('-T', '--tex', action='store_true', help='Use LaTeX to render all text')
p.add_option('-n', '--n-cycles', action='store_true', help='Show the number of predictor-corrector cycles in legend')
p.add_option('-L', '--legend', action='store_true', help='Draw a figure legend')

opt, filenames = p.parse_args()

if opt.tex:
    print 'Using LaTeX fonts'
    rc('text', usetex=True)
    rc('font',**{'family':'serif','serif':['Computer Modern']})

titles = np.array([['RMS error - extrapolation', 'RMS error - interpolation'],
                   ['Max error - extrapolation', 'Max error - interpolation']])

if opt.label is None:
    opt.label = filenames
else:
    if len(opt.label) != len(filenames):
        p.error('Number of labels does not match number of input files')


fig = plt.figure()
fig.text(0.5, 0.95, opt.title, ha='center', size='x-large')

axes_grid = np.zeros((2,2), dtype=object)
for col in (0,1):
    for row in (0,1):
        axes_grid[row,col] = plt.subplot2grid((2,2), (row,col))

colors = itertools.cycle(['r', 'g', 'b', 'c', 'm', 'y', 'k'])
leg_lines = []
leg_text = []

for filename, label, color in zip(filenames, opt.label, colors):
    print "processing file", filename
    print "label", label

    extrap_steps = int(os.popen('grep extrapolate_steps %s' % filename).read().split()[2])
    extrap_data = np.loadtxt(os.popen('grep "^E err" %s' % filename), usecols=(2,3,4))
    interp_data = np.loadtxt(os.popen('grep "^I err" %s' % filename), usecols=(2,3,4))

    for col, (data, cycle_label) in enumerate(zip((extrap_data, interp_data),
                                                  ("extrap", "interp"))):
        time = data[:extrap_steps,0] - data[0,0]

        cycles = []
        for start in np.arange(0,len(data),extrap_steps):
            cycle = data[start:start+extrap_steps,1:]
            cycles.append(cycle)

        if len(cycles[-1]) != extrap_steps:
            cycles = cycles[:-1]

        print  ("n_%s = %d " % (cycle_label, len(cycles)))

        if len(cycles) == 0:
            continue

        all_data = np.dstack(cycles)

        max_data = all_data.max(axis=2)
        min_data = all_data.min(axis=2)
        mean_data = all_data.mean(axis=2)
        std_data  = all_data.std(axis=2)/np.sqrt(len(cycles))  # error in mean, so divide by sqrt(N)

        for row in (0, 1):
            ax = axes_grid[row,col]
            mean_plt = ax.errorbar(time, mean_data[:,row], std_data[:,row], None, color+'.-', lw=1, mec=color)

            if row == 0 and col == 0:
                leg_lines.append(mean_plt[0])
                if opt.n_cycles:
                    label += ' (%d cycles)' % len(cycles)
                leg_text.append(label)

            ax.set_title(titles[row,col], size='medium')

            if row == 1:
                ax.set_xlabel('Time / fs')

            if col == 0:
                ax.set_ylabel('Force error / eV/ A')

    print

if opt.legend:
    fig.subplots_adjust(top=0.87, bottom=0.05+0.05*len(leg_lines)/2, hspace=0.24)
    leg = fig.legend(leg_lines, leg_text, 'lower center', ncol=2, prop={'size':'small'})
    leg.draw_frame(False)

if is_interactive_shell():
    plt.show()
else:
    plt.savefig(opt.output)
    if not opt.no_display:
        os.system('open %s' % opt.output)

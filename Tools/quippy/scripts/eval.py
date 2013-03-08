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

import sys
import optparse
import os
import multiprocessing
from quippy import *
from numpy import *

p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-i', '--init-args', action='store', help="""Arguments used to initialise Potential""")
p.add_option('-p', '--param-file', action='store', default=None, help="""XML parameter file""")
p.add_option('-c', '--calc-args', action='store', default='', help="""Calc args to pass to Potential.calc()""")
p.add_option('-o', '--output', action='store', default='', help="""File to write results to""")
p.add_option('-P', '--parallel', action='store_true', help="""Launch calculations in parallel, with one process per input frame""")
p.add_option('-C', '--chdir', action='store_true', help="""Create a new subdirectory for each configuration""")

opt, infiles = p.parse_args()

pot = Potential(opt.init_args, param_filename=opt.param_file)

def do_calc(at, basename='eval'):
    print 'Preparing job %d' % at.i
    if opt.chdir:
        rundir = '%s-%05d' % (basename, at.i)
        origdir = os.getcwd()
        if not os.path.exists(rundir):
            os.mkdir(rundir)
        os.chdir(rundir)
    try:
        pot.calc(at, args_str=opt.calc_args)
        print 'Job %i completed' % at.i
        return at
    finally:
        if opt.chdir:
            os.chdir(origdir)

if __name__ == '__main__':

    configs = AtomsList(infiles)
    for i, at in enumerate(configs):
        at.params['i'] = i

    processes = 1
    if opt.parallel:
        processes = len(configs)

    pool = multiprocessing.Pool(processes=processes)
    results = pool.map(do_calc, configs)

    pool.close()
    pool.join()

    results = AtomsList(results)
    results.write(opt.output)

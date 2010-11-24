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

from quippy import *
from numpy import *
import sys, optparse, shutil, os
import atomeye

p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-f', '--frame', action='store', help='Frame number', type='int', default=-1)
p.add_option('-o', '--outfile', action='store', help='Output file', default='snap.png')
p.add_option('-w', '--window', action='store_true', help='Show AtomEye viewer window', default=False)
p.add_option('-l', '--load-view', action='store', help='Load view from AtomEye command script')
p.add_option('-p', '--property', action='store', help="""Property to use to colour atoms (default none)""")
p.add_option('-a', '--arrows', action='store', help="""Property to use to draw arrows (default none)""")
p.add_option('-e', '--exec_code', action='store', help="""Python code to execute before writing output file. Atoms object is available as `at`.""")
p.add_option('-W', '--width', action='store', help="""Width of output movie, in pixels.""", type='int')
p.add_option('-H', '--height', action='store', help="""Height of output movie, in pixels.""", type='int')
p.add_option('-A', '--aspect', action='store', help="""Aspect ratio. Used if only one of --width or --height is given. Default 0.75.""", default=0.75, type='float')
p.add_option('-c', '--centre', action='store', help="Atom index or position on which to centre view")
p.add_option('-s', '--shift', action='store', help="Shift to apply to crystal")

opt, args = p.parse_args()

if len(args) == 0:
   p.error('No input files specified')

at = Atoms(args[0], frame=opt.frame)

if opt.exec_code is not None:
    exec(opt.exec_code)

if opt.centre is not None:
    opt.centre = eval(opt.centre)

if opt.shift is None:
    opt.shift = [.5, .5, .5]
else:
    opt.shift = eval(opt.shift)

at.write(opt.outfile, width=opt.width, height=opt.height, aspect=opt.aspect,
         script=opt.load_view, centre=opt.centre, shift=opt.shift, property=opt.property,
         arrows=opt.arrows, nowindow=not opt.window)

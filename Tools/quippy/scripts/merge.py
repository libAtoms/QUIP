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
import sys, optparse, shutil, os, itertools

p = optparse.OptionParser(usage='%prog [options] <input file>...')

p.add_option('-o', '--outfile', action='store', help='Output file (default stdout)', default='stdout')
p.add_option('-f', '--format', action='store', help='Output format (inferred from filename by default)')
p.add_option('-O', '--overwrite', action='store_true', help="""If given, overwrite properties and parameters (moving through input files from left to right).
                                                                Default is off, so first property or parameter found takes precedence in name collisions
                                                                (unless --prefixes argument is given).""")
p.add_option('-p', '--prefixes', action='append', help="Prefixes to prepend to properties and parameters. One argument per input file.")
p.add_option('-x', '--exclude', action='append', help="Property or parameter names to exclude. Filtering is done after addition of prefixes.")
p.add_option('-P', '--prefix-names', action='append', help="Names of properties and parameters to add prefixes to. Default is [force, energy, virial].")


opt, args = p.parse_args()

if len(args) < 1:
    p.error("At least one input file is required")

if opt.format is None and opt.outfile is not None:
      opt.format = os.path.splitext(opt.outfile)[1][1:]

if opt.format is None or opt.format == '':
   opt.format = 'xyz'

sources = [AtomsList(f, store=False) for f in args]
outfile = AtomsWriter(opt.outfile, format=opt.format)

prefixes = ['' for s in sources]
if opt.prefixes is not None:
    if len(opt.prefixes) > len(args):
        p.error("Number of prefixes (%d) cannot exceed number of input files (%d)" % (len(prefixes), len(args)))
    prefixes[:len(opt.prefixes)] = opt.prefixes

if opt.prefix_names is None:
    opt.prefix_names = ['force', 'energy', 'virial']

for frame, confs in enumerate(itertools.izip(*sources)):

    at_out = Atoms(n=confs[0].n, lattice=confs[0].lattice, properties={})
    
    for prefix, at in zip(prefixes, confs):

        if at.n != at_out.n:
            p.error("Mismatch in number of atoms at frame %d", frame)

        if abs(at.lattice - at_out.lattice).max() > 1e-8:
            p.error("Mismatch in lattice at frame %d", frame)

        for k in at.params.keys():
            k_new = k.lower()
            if k_new in opt.prefix_names: k_new = prefix+k_new
            if k_new in opt.exclude: continue
            if not opt.overwrite and k_new in at_out.params.keys(): continue
            at_out.params[k_new] = at.params[k]

        for k in at.properties.keys():
            k_new = k.lower()
            if k_new in opt.prefix_names: k_new = prefix+k_new
            if k_new in opt.exclude: continue
            at_out.add_property(k_new, at.properties[k], property_type=at.properties.get_type(k), overwrite=opt.overwrite)

    outfile.write(at_out)

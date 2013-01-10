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
import sys, optparse, os


p = optparse.OptionParser(usage='%prog [options] ( <input file> [ <output file> ] | <input file> [ <input file> ... ] [ (-o|--output) <output file> ] )')

p.add_option('-r', '--range', action='store', help="""Range of frames to include. Should be either a single frame
number or a slice [start]:[stop][:step]. If -r is omitted,
default is all frames. Frames start from zero. Negative indices count backwards from the end of the file,
with -1 being the last frame.""")
p.add_option('-L', '--lattice', action='store', help="""Override lattice with LATTICE, given in
Fortran ordering as "R11 R21 R31 R12 R22 R32 R13 R23 R33". If
three fields are given then the lattice is assumed to be cubic with
a=R11, b=R22, c=R33 and all other components zero. One field implies
a cubic cell with a=b=c=R11.""")
p.add_option('-p', '--properties', action='store', help="""Comma or colon separated list of properties (i.e. per-atom variables)
that should be included in output file, e.g.
  -p species,pos,velo
The default is to print all properties.
For a valid XYZ file, the first should be "species" and the second
should be something position-like. The order of properties is as
specified in this command. If list is the string 'SAME', output
properties of each frame are same as input properties.""")
p.add_option('-P', '--params', action='store', help="""Comma or colon separated list of parameters (per-frame variables) to
include in output file. For example:
  -P Energy,MaxForce
The default is to output all parameters.
For XYZ output these parameters are in addition to the special
parameters "Lattice" and "Properties" which are always written.
The original order of parameters in the input file is preserved.""")
p.add_option('-f', '--format', action='store', help="""Explicitly specify output format, e.g. --format=xyz
Supported formats: %s.""" % ', '.join([s for s in AtomsWriters.keys() if isinstance(s, str)]))
p.add_option('-I', '--in-format', action='store', help="""Explicitly specify input format, e.g. --in-format=xyz
Supported formats: %s.""" % ', '.join([s for s in AtomsReaders.keys() if isinstance(s, str)]))
p.add_option('-m', '--merge', action='store', help="""Merge two input files. An auxilary input file name should be given.""")
p.add_option('-M', '--merge-properties', action='store', help="""List of properties to overwrite from MERGE file. Default is all properties.""")
p.add_option('-g', '--merge-params', action='store_true', help="""Merge params from MERGE file into output file.""", default=False)
p.add_option('-x', '--extract-params', action='store_true', help="""Read only the parameters from input files. If -N is not given, prints a table of params, one frame per line.""", default=False)
p.add_option('-F', '--extract-format', action='store', help="""Format string used to print parameters if the --extract-params option is used.""")
p.add_option('-e', '--exec-code', action='store', help="""Python code to execute on each frame before writing it to output file. Atoms object is
available as `at`, and do_print is set to True. If the user-supplied code sets do_print to False, the frame is not printed.""")
p.add_option('-E', '--exec-code-file', action='store', help="""File with python code to execute, just like -e/--exec-code.""")
p.add_option('-B', '--exec-before', action='store', help="""Python code to execute after looping over frames""")
p.add_option('-A', '--exec-after', action='store', help="""Python code to execute after looping over frames""")
p.add_option('--atoms-ref', action='store', help="""Reference configuration for reordering atoms. Applies to CASTEP file formats only.""")
p.add_option('-v', '--verbose', action='store_true', help="""Verbose output (first frame only)""", default=False)
p.add_option('-n', '--rename', action='append', help="""Old and new names for a property or parameter to be renamed. Can appear multiple times.""", nargs=2)
p.add_option('-R', '--remove', action='append', help="Remove a property or param. Can appear multiple times")
p.add_option('-s', '--select', action='store', help="""Output only a subset of the atoms in input file. Argument should resolve to logical mask.""")
p.add_option('--int-format', action='store', help="""Format string to use when writing integers in XYZ format.""")
p.add_option('--real-format', action='store', help="""Format string to use when writing real numbers in XYZ format.""")
p.add_option('-N', '--no-print-at', action='store_true', help="""Suppress printing of Atoms object (useful when also using -e argument).""", default=False)
p.add_option('-o', '--output', action='store', help="""File to output to (required when more than 1 input file is listed)""")
p.add_option('--extra-args', action='store', help="""Extra arguments to be passed when constructing AtomsWriter""")
p.add_option('--write-args', action='store', help="""Extra arguments to be passed to write() routine.""")
p.add_option('--read-args', action='store', help="""Extra arguments to be passed to read() routine.""")
p.add_option('-a', '--atom-range', action='store', help="""Range of atoms to include: should be followed by two arguments, min
and max 1-based atom indices. Using a reduced atom range dramatically
speeds up reading of large XYZ trajectories.""", nargs=2)
p.add_option('-l', '--load-all', action='store_true', help="""Read all frames before starting processing. Allows frame indexing for file types which do not support random access (e.g. .castep)""")
p.add_option('-t', '--time-ordered-series', action='store_true', help="""Join all input files, ordering by time and discarding duplicates""")

# Options related to rendering of images with AtomEye
p.add_option('-V', '--view', action='store', help='Load view from AtomEye command script')
p.add_option('--property', action='store', help="""Property to use to colour atoms (default none)""")
p.add_option('--arrows', action='store', help="""Property to use to draw arrows (default none)""")
p.add_option('-W', '--width', action='store', help="""Width of output movie, in pixels.""", type='int')
p.add_option('-H', '--height', action='store', help="""Height of output movie, in pixels.""", type='int')
p.add_option('--aspect', action='store', help="""Aspect ratio. Used if only one of --width or --height is given. Default 0.75.""", type='float')
p.add_option('-c', '--centre', action='store', help="Atom index or position on which to centre view")

# CASTEP specific output options
p.add_option('--cell-template', action='store', help='Template .cell file, to be used when writing CASTEP .cell files')
p.add_option('--param-template', action='store', help='Template .param file, to be used when writing CASTEP .cell files')

opt, args = p.parse_args()

# check for conflicts with -o|--output
if opt.no_print_at or opt.extract_params:
    if opt.output is not None:
        p.error('--no_print_at/--extract-params can not coexist with --output')
    if len(args) < 1:
        p.error('At least one input file must be specified')
    infiles = args
else: # didn't specify no_print_at
    if opt.output is None:
        if len(args) != 2:
            p.error('With no --output, exactly one input and one output file must be specified (use /dev/null or NONE for no output).')
        outfile = args.pop()
        infiles = args
    else:
        outfile = opt.output
        infiles = args

# make no-output outfile name standard
if opt.extract_params or opt.no_print_at or outfile.upper() == 'NONE' or outfile == '/dev/null' or outfile == 'dev_null':
    outfile = None

# convert - to stdin/stdout
infiles = [ f == '-' and 'stdin' or f for f in infiles ]
if outfile == '-': outfile = 'stdout'

# check for existing outfile
if opt.output is None and not (opt.no_print_at or opt.extract_params):
    if os.path.exists(outfile):
        p.error('Output file %s specified without -o|--output already exists. Use (-o|--output) filename to overwrite.' % outfile)

# check for proper format of --range
if opt.range is not None:
    try:
        opt.range = parse_slice(opt.range)
    except:
        p.error('Cannot parse slice "%s" - should be in format [start]:[stop][:step]')
else:
    # Default is all frames
    opt.range = slice(0, None, None)

if isinstance(opt.range, int):
    if opt.range >= 0:
        opt.range = slice(opt.range, opt.range+1,+1)
    else:
        opt.range = slice(opt.range, opt.range-1,-1)


if opt.lattice is not None:
    opt.lattice = [ float(x) for x in opt.lattice.split() ]
    if len(opt.lattice) == 9:
        opt.lattice = farray(opt.lattice).reshape((3,3),order='F').T
    elif len(opt.lattice) == 3:
        opt.lattice = farray(diag(opt.lattice))
    elif len(opt.lattice) == 1:
        opt.lattice = opt.lattice*fidentity(3)
    else:
        p.error('LATTICE should consist of 1, 3 or 9 numbers -- got %r' % opt.lattice)


print_same_properties=False
if opt.properties is not None:
    if opt.properties == 'SAME':
        print_same_properties=True
    else:
        opt.properties = parse_comma_colon_list(opt.properties)

if opt.params is not None:
    opt.params = parse_comma_colon_list(opt.params)

if opt.atoms_ref is not None:
    opt.atoms_ref = Atoms(opt.atoms_ref)

if opt.merge is not None:
    if opt.atoms_ref is not None:
        merge_configs = AtomsReader(opt.merge, atoms_ref=opt.atoms_ref)
    else:
        merge_configs = AtomsReader(opt.merge)

    if opt.merge_properties is not None:
        opt.merge_properties = parse_comma_colon_list(opt.merge_properties)
    else:
        at_merge = merge_configs[0]
        opt.merge_properties = at_merge.properties.keys()

def process(at, frame):
    if print_same_properties:
        write_args['properties'] = at.properties.keys()

    # filter atoms
    if opt.select is not None:
        at2 = at.select(mask=eval(opt.select))
        at = at2

    # Override lattice
    if opt.lattice is not None:
        at.set_lattice(opt.lattice, False)

    # Merge from merge_config
    if opt.merge:
        try:
            at_merge = merge_configs[frame]
        except IndexError:
            at_merge = merge_configs[0]
        for k in opt.merge_properties:
            at.add_property(k, at_merge.properties[k], property_type=at_merge.properties.get_type(k), overwrite=True)

        if opt.merge_params is not None:
            at.params.update(at_merge.params)

    # Execute user code
    do_print = True
    if opt.exec_code_file is not None:
        execfile(opt.exec_code_file)
    if opt.exec_code is not None:
        exec(opt.exec_code)

    if opt.centre is not None:
        write_args['centre'] = eval(opt.centre)

    # Filter parameters
    if opt.params is not None:
        for k in at.params.keys():
            k = k.lower()
            if not k in opt.params:
                del at.params[k]

    # Rename properties and parameters
    if opt.rename is not None:
        for (old, new) in opt.rename:
            if old in at.properties:
                at.properties[new] = at.properties[old]
                del at.properties[old]
            elif old in at.params:
                at.params[new] = at.params[old]
                del at.params[old]
            else:
                raise AttributeError('Cannont rename: no property or parameter named "%s" exists' % old)

    # Remove properties and parameters
    if opt.remove is not None:
        for remove in opt.remove:
            if remove in at.properties:
                del at.properties[remove]
            elif remove in at.params:
                del at.params[remove]
            else:
                raise AttributeError('Cannot remove: no property or parameter named "%s" exists' % old)

    # Verbose output
    if opt.verbose and frame == 0:
        print 'N_ATOMS', at.n
        print 'PROPERTIES:', at.properties.keys()
        print 'PARAMS:', at.params.keys()

    # Do the writing
    if opt.extract_params:
        if not opt.no_print_at:
            if frame == 0:
                print '#' + ' '.join(at.params.keys())
            for k in at.params.keys():
                if opt.extract_format:
                    print(opt.extract_format % at.params[k]),
                else:
                    print at.params[k],

            print
        else:
            pass
    elif do_print and outfile is not None:

        if opt.format in ('eps', 'png', 'jpg') and isinstance(opt.range, slice):
            write_args['frame'] = frame

        if opt.format == 'cell':
            write_args['cell_template'] = opt.cell_template
            write_args['param_template'] = opt.param_template

        if opt.properties is None:
            outfile.write(at, **write_args)
        else:

            # Convert from frac_pos to pos
            if 'frac_pos' in opt.properties and not at.has_property('frac_pos'):
                at.add_property('frac_pos', 0.0, n_cols=3)
                at.frac_pos[:] = dot(at.g, at.pos)

            # or vice-versa
            if 'pos' in opt.properties and not at.has_property('pos'):
                at.add_property('pos', 0.0, n_cols=3)
                at.pos[:] = dot(at.lattice, at.frac_pos)

            try:
                # Try to do the filtering at the writing stage
                outfile.write(at, **write_args)
            except TypeError:
                p.error('Cannot specify property filtering when writing to file "%s"' % outfile)


# Build dictionaries of arguments for AtomsWriter constructor
# and for write() method
init_arg_rename = {'view': 'script'}
init_args = {}
for arg in ('width', 'height', 'aspect', 'view'):
    if getattr(opt, arg) is not None:
        initarg = init_arg_rename.get(arg, arg)
        init_args[initarg] = getattr(opt, arg)
write_arg_rename = {}
write_args = {}
for arg in ('properties', 'real_format', 'int_format', 'property', 'arrows'):
    if getattr(opt, arg) is not None:
        writearg = write_arg_rename.get(arg, arg)
        write_args[writearg] = getattr(opt, arg)

if opt.extra_args is not None:
    init_args.update(eval('dict(%s)' % opt.extra_args))

if opt.write_args is not None:
    write_args.update(eval('dict(%s)' % opt.write_args))

if opt.format is None and outfile is not None:
    opt.format = os.path.splitext(outfile)[1][1:]

if opt.format is None or opt.format == '':
    opt.format = 'xyz'

if opt.aspect is None:
    opt.aspect = 0.75

stdout = False
if outfile is not None:
    stdout = outfile == 'stdout'
    try:
        outfile = AtomsWriter(outfile, format=opt.format, **init_args)
    except RuntimeError, re:
        p.error(str(re))

read_args = {}
if opt.atoms_ref is not None:
    read_args['atoms_ref'] = opt.atoms_ref
if opt.atom_range is not None:
    read_args['range'] = opt.atom_range
if opt.extract_params and 'stdin' not in infiles:
    read_args['range'] = 'empty'
if opt.read_args is not None:
    read_args.update(eval("dict(%s)" % opt.read_args))

if opt.time_ordered_series:
    infiles = time_ordered_series(infiles)

if opt.load_all:
    all_configs = AtomsList(infiles,
                            format=opt.in_format,
                            **read_args)[opt.range.start:opt.range.stop:opt.range.step]
else:
    all_configs = AtomsReader(infiles,
                              format=opt.in_format,
                              start=opt.range.start,
                              stop=opt.range.stop,
                              step=opt.range.step,
                              **read_args)

try:
    show_progress = not opt.extract_params and not stdout and not opt.no_print_at and len(all_configs) > 1 and sys.stderr.isatty() and not opt.load_all
except AttributeError:
    show_progress = False

if show_progress:
    from quippy.progbar import ProgressBar
    pb = ProgressBar(0,len(all_configs),80,showValue=True)

if opt.exec_before is not None:
    exec(opt.exec_before)
for i, at in enumerate(all_configs):
    try:
        process(at, i)
    except (IndexError, ValueError, RuntimeError), re:
        p.error(str(re))
    if show_progress: pb(i)

if show_progress:
    print

if opt.exec_after is not None:
    exec(opt.exec_after)

if outfile is not None:
    try:
        outfile.close()
    except AttributeError:
        pass

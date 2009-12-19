#!/usr/bin/env python

from quippy import *
import sys, optparse, itertools

p = optparse.OptionParser(usage='%prog [options] <input file> [<output file>]')

p.add_option('-r', '--range', action='store', help="""Range of frames to include. Should be either a single frame
number or a slice [start]:[stop][:step]. If -r is omitted,
default is all frames. We use Fortran indexing, so slices include the
"stop" frame. Negative indices count backwards from the end of the file,
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
specified in this command.""")
p.add_option('-P', '--params', action='store', help="""Comma or colon separated list of parameters (per-frame variables) to
include in output file. For example:
  -P Energy,MaxForce
The default is to output all parameters.
For XYZ output these parameters are in addition to the special
parameters "Lattice" and "Properties" which are always written.
The original order of parameters in the input file is preserved.""")
p.add_option('-f', '--format', action='store', help="""Explicitly specify output format, e.g. --format=xyz
Supported formats: cell, nc, pos, pov, xyz.""")
p.add_option('-m', '--merge', action='store', help="""Merge two input files. An auxilary input file name should be given.""")
p.add_option('-M', '--merge-properties', action='store', help="""List of properties to overwrite from MERGE file. Default is all properties.""")
p.add_option('-e', '--exec_code', action='store', help="""Python code to execute on each frame before writing it to output file. Atoms object is
available as `at`.""")

opt, args = p.parse_args()

if len(args) != 1 and len(args) != 2:
   p.error('One input file must be specified. Output file is optional.')

try:
   infile, outfile = args
except ValueError:
   infile, = args
   outfile = None

if infile == '-':  outfile = 'stdin'
if outfile == '-': outfile = 'stdout'

class SliceParser(object):
   def __getitem__(self, idx):
      return idx

def parse_comma_colon_list(L):
   if ':' in L:
      L = L.split(':')
   elif ',' in L:
      L = L.split(',')
   else:
      L = [L]

   return [k.lower() for k in L]


if opt.range is not None:
   try:
      opt.range = eval('SliceParser()[%s]' % opt.range)
   except:
      p.error('Cannot parse slice "%s" - should be in format [start]:[stop][:step]')
else:
   # Default is all frames
   opt.range = slice(1, None, None)

if opt.lattice is not None:
   opt.lattice = [ float(x) for x in opt.lattice.split() ]
   if len(opt.lattice) == 9:
      opt.lattice = farray(opt.lattice).reshape((3,3),order='F')
   elif len(opt.lattice) == 3:
      opt.lattice = farray(diag(opt.lattice))
   elif len(opt.lattice) == 1:
      opt.lattice = opt.lattice*fidentity(3)
   else:
      p.error('LATTICE should consist of 1, 3 or 9 numbers -- got %r' % opt.lattice)


if opt.properties is not None:
   opt.properties = parse_comma_colon_list(opt.properties)

if opt.params is not None:
   opt.params = parse_comma_colon_list(opt.params)

if opt.merge is not None:
   merge_config = Atoms(opt.merge)

   if opt.merge_properties is not None:
      opt.merge_properties = parse_comma_colon_list(opt.merge_properties)
   else:
      opt.merge_properties = merge_config.properties.keys()

def process(at):
   # Override lattice
   if opt.lattice is not None:
      at.set_lattice(opt.lattice)

   # Filter parameters
   if opt.params is not None:
      for k in at.params.keys():
         k = k.lower()
         if not k in opt.params:
            del at.params[k]

   # Merge from merge_config
   if opt.merge:
      for k in opt.merge_properties:
         at.add_property(k, getattr(merge_config, k))

   # Execute user code
   if opt.exec_code is not None:
      exec(opt.exec_code)
      
   # Do the writing
   if outfile is not None:
      if opt.properties is None:
         outfile.write(at)
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
            outfile.write(at, properties=opt.properties)
         except TypeError:
            p.error('Cannot specify property filtering when writing to file "%s"' % outfile)


all_configs = AtomsList(infile)
if outfile is not None:
   outfile = AtomsWriter(outfile, format=opt.format)

if len(all_configs) == 1:
   opt.range = 1

if isinstance(opt.range, slice):
   # multiple frames

   from quippy.progbar import ProgressBar
   pb = ProgressBar(0,len(frange(*opt.range.indices(len(all_configs)))),80,showValue=True)
   
   for i, at in fenumerate(itertools.islice(all_configs, opt.range.start+1, opt.range.stop, opt.range.step)):
      process(at)
      pb(i)

   print
                    
else:
   # single frame
   process(AtomsList(infile)[opt.range])

if outfile is not None:
   try:
      outfile.close()
   except AttributeError:
      pass


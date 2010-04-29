#!/usr/bin/env python
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HND X
# HND X   libAtoms+QUIP: atomistic simulation library
# HND X
# HND X   Portions of this code were written by
# HND X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# HND X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# HND X
# HND X   Copyright 2006-2010.
# HND X
# HND X   Not for distribution
# HND X
# HND X   Portions of this code were written by Noam Bernstein as part of
# HND X   his employment for the U.S. Government, and are not subject
# HND X   to copyright in the USA.
# HND X
# HND X   When using this software, please cite the following reference:
# HND X
# HND X   http://www.libatoms.org
# HND X
# HND X  Additional contributions by
# HND X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# HND X
# HND XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

# This script can be used to interface QUIP to CASTEP. Within QUIP, you
# should create a FilePot potential like this:
# 
# type(Potential) :: pot
# ...
# call initialise(pot, "FilePot command=./castep_driver.py")
# 
# To use the CASTEP driver, you must have Python 2.3 or later
# (http://www.python.org) installed, along with the NumPy numerical
# extension module (http://www.numpy.org). You also need the 'PyAtoms'
# python package somewhere on your PYTHONPATH. This package can be
# found in the 'pyatoms' SVN module in the repository.
# 
# This FilePot script invokes CASTEP to get forces and energy of a
# cluster. This script converts XYZ input file to castep cell file,
# then uses template .cell and .param files for the rest of the CASTEP
# parameters. The atom number mangling necessary to convert between
# simple atom number 1..N and the (species, number within species)
# numbering that CASTEP uses is dealt with by this script. 
#
# There are some parameters at the top of the file castep_driver.py
# which you should edit. The most important are TEMPLATE, which
# specifies the template .cell and .param files to be used (see the
# examples in example.cell and example.param in this directory) and
# CASTEP which is the path to the CASTEP executable to be used. Since
# the calculation takes place in a subdirectory of where the QUIP
# program is running, this will typically begin with '..'. 
# 
# Both these parmeters can also be set by environment variables,
# named CASTEP_TEMPLATE and CASTEP respectively.
# 
# The CASTEP parameter should contain a %s where the seed name should
# go, e.g.
# 
# export CASTEP="../castep %s"
# 
# or
# 
# NUMNODES=8
# export CASTEP="mpirun -np $NUMNODES ../run_castep %s"
# 
# Any questions, email James Kermode <jrk33@cam.ac.uk>


#----------------------------------------------------------------
# Dependancies 
#----------------------------------------------------------------

# Python standard library modules
import sys, string, os, os.path, shutil, glob, operator, xml.dom.minidom

# NumPy <http://www.numpy.org>
from numpy import *

# PyAtoms package
from pyatoms import *

#----------------------------------------------------------------
# Parameters
#----------------------------------------------------------------

# Template used for .cell and .param files
if os.environ.has_key('CASTEP_TEMPLATE'):
   CASTEP_TEMPLATE = os.environ['CASTEP_TEMPLATE']
else:
   CASTEP_TEMPLATE = 'castep'

# Command used to execute castep, with a %s where seed name should go
if os.environ.has_key('CASTEP'):
   CASTEP = os.environ['CASTEP']
else:
   CASTEP = 'mpiexec -np 64 ../run_castep %s'

# If there's no %s, put seed name at end of string
if CASTEP.find('%s') == -1:
   CASTEP = CASTEP + ' %s'

# Should we try to use .check files to speed up calculations?
USE_CHECK_FILES = True
if 'CASTEP_USE_CHECK_FILES' in os.environ:
   USE_CHECK_FILES = True

# Keep a copy of all check files, for debugging purposes
SAVE_ALL_CHECK_FILES = False
if 'CASTEP_SAVE_ALL_CHECK_FILES' in os.environ:
   SAVE_ALL_CHECK_FILES = True

# Keep a copy of all cell and param files, for debugging purposes
SAVE_ALL_INPUT_FILES = False
if 'CASTEP_SAVE_ALL_INPUT_FILES' in os.environ:
   SAVE_ALL_INPUT_FILES = True

# If set to True, don't actually run CASTEP
TEST_MODE = False

# If set to True, create a queue of calculations for later execution
BATCH_QUEUE = False

# If set to True, read from previously calculated queue
BATCH_READ  = False

# If any force is larger than this threshold, repeat the
# calculation without checkfile (units: eV/A)
MAX_FORCE_THRESHOLD = 150.0

# Working directory for CASTEP. Set this to a local scratch
# directory if network file performance is poor.
WORKING_DIR = '.'

# Should we hash cell file and name directory accordingly?
DO_HASH = False

# Number of decimal places to include in hash of cell
HASH_NDIGITS = 2

#----------------------------------------------------------------
# Subroutines 
#----------------------------------------------------------------

def die(message):
   "Print error message and abort"
   sys.stdout.write('castep_driver ABORT   : '+message+'\n')
   os.chdir(orig_dir)
   sys.exit(1)

def warn(message):
   "Print message to STDOUT"
   sys.stdout.write('castep_driver WARNING : '+message+'\n')

def error(message):
   "Print message to STDOUT"
   sys.stdout.write('castep_driver ERROR   : '+message+'\n')
   

def info(message):
   "Print message to STDOUT"
   sys.stdout.write('castep_driver INFO    : '+message+'\n')


def castep_check_pspots(cluster, cell, param):
   """Check pseudopotential files are present, and that we have one for
   each element present in cluster. Also tries to check spin polarisation
   of system matches that specified in parameters."""

   # Get list of unique elements in cluster
   pspot_dict = {}
   for el in cluster.species:
      pspot_dict[el] = None
   elements = pspot_dict.keys()

   if 'SPECIES_POT' in cell:
      for line in cell['SPECIES_POT']:
         element, pspot = map(string.strip, line.split())
         pspot_dict[element] = pspot

   valence_charges = {}
   missing_pspots = False
   
   for el in elements:
      pspot = pspot_dict[el]
      if pspot is None:
         # Not fatal since on-the-fly pseudo-potential will be generated
         # however, we can't estimate then n_electrons here
         warn('Pseudopotential missing for element %s' % el)
         valence_charges[el] = '?'
         missing_pspots = True
      elif not pspot.startswith('/"'):
         shutil.copy(orig_dir+'/'+pspot, '.')
         pspot_lines = open(pspot,'r').readlines()

         # First check for old-style pspot report
         zlines = filter(lambda s: s.startswith('    |  z ='), pspot_lines)
         if len(zlines) == 1:
            # Extract number of valence electrons
            fields = zlines[0].split()
            zv = float(fields[fields.index('zv')+2])
         else:
            # Now try newer style for OTF generated pspots
            zlines = filter(lambda s: s.startswith('   | Element: '), pspot_lines)
            if len(zlines) != 1:
               raise ValueError('Malformed pseudopotential file: does not contain pseudopotential report')
            fields = zlines[0].split()
            zv = float(fields[fields.index('charge:')+1])

         valence_charges[el] = zv
      else:
         # It's an on-the-fly pseudopotential
         fields = pspot[2:-2].split('|')
         valence_charges[el] = float(fields[0])

   info('unique elements and valence electrons: %s' % valence_charges)

   if not missing_pspots:
      n_electrons = reduce(operator.add, \
                           map(lambda el: valence_charges[el]*count(cluster.species == el), \
                               elements))

      info('total electrons %.1f' % n_electrons)
      if (param.has_key('spin_polarised') and param['spin_polarised'].lower() == 'false' and \
          int(n_electrons) % 2 != 0):
         raise ValueError('spin polarised set to false but got odd number of electrons!')
      if (param.has_key('spin_polarised') and param['spin_polarised'].lower() == 'true' and \
          int(n_electrons) % 2 == 0):
         raise ValueError('spin polarised set to true but got even number of electrons!')


def castep_run(cell, param,  stem, castep, log=None):
   "Invoke castep and return True if it completed successfully"
   # Write parameter file ...
   param.write(stem+'.param')
   # ... and cell file
   cell.write(stem+'.cell')

   if TEST_MODE:
      info('test mode: not running castep')
      
      # Simulate making check file
      check_file = open(stem+'.check','w')
      check_file.close()

   else:
      # Remove old output file and error files
      try: 
         os.remove(stem+'.castep')
      except:
         pass
      for f in glob.glob('%s*.err' % stem):
         os.remove(f)

      os.system(castep % stem)

   if SAVE_ALL_CHECK_FILES:
      if os.path.exists('%s.check' % stem):
         n = 0
         while os.path.exists('%s.check.%d' % (stem, n)):
            n += 1
         shutil.copyfile('%s.check' % stem, '%s.check.%d' % (stem, n))

   if SAVE_ALL_INPUT_FILES:
      if os.path.exists('%s.cell' % stem):
         n = 0
         while os.path.exists('%s.cell.%d' % (stem, n)):
            n += 1
         shutil.copyfile('%s.cell' % stem, '%s.cell.%d' % (stem, n))
      if os.path.exists('%s.param' % stem):
         n = 0
         while os.path.exists('%s.param.%d' % (stem, n)):
            n += 1
         shutil.copyfile('%s.param' % stem, '%s.param.%d' % (stem, n))

   err_files = glob.glob('%s*.err' % stem)
   got_error = False
   if (len(err_files) > 0):
      for f in err_files:
         error_text = open(f).read().strip()
         if error_text != '':
            got_error = True
            error(error_text)

   # Write log file here so that it will always be written
   if log is not None and os.path.exists('%s.castep' % stem):
      logf = open(log, 'a')
      castep_output = open('%s.castep' % stem, 'r').readlines()
      logf.writelines(castep_output)
      logf.close()

   return not got_error


import md5

def hash_atoms(at, ndigits=HASH_NDIGITS):
   """Hash an atoms object with a precision of ndigits decimal digits.
   Species, lattice and fractional positions are fed to MD5 to form the hash."""

   doround = lambda x: round(x, ndigits)

   # Compute fractional positions, round them to ndigits, then sort them
   # for hash stability
   flat_frac_pos = array([ dot(at.pos[i,:],at.g) for i in range(at.n) ]).flatten()
   flat_frac_pos = map(doround, flat_frac_pos)
   flat_frac_pos.sort()

   m = md5.new()
   m.update(str(map(doround, at.lattice.flat)))
   m.update(str(at.species))
   m.update(str(flat_frac_pos))

   return m.hexdigest()


class ParamError(Exception):
   pass


#----------------------------------------------------------------
# Main code
#----------------------------------------------------------------

# Save starting directory
orig_dir = os.getcwd()

if len(sys.argv) < 3:
   die('Usage: [-t|-q|-r] %s <xyzfile> <outputfile>' % sys.argv[0])

# If first command line option is '-t' set TEST_MODE to true
if sys.argv[1] == '-t':
   TEST_MODE = True
   sys.argv.pop(0)

# If first command line option  is '-q' set BATCH_QUEUE to true
if sys.argv[1] == '-q':
   BATCH_QUEUE = True
   sys.argv.pop(0)

if sys.argv[1] == '-r':
   BATCH_READ = True
   sys.argv.pop(0)

xyzfile = sys.argv[1]
outfile = sys.argv[2]


stem = os.path.basename(xyzfile)
if stem[-4:] == '.xyz': # Remove extension
   stem = stem[:-4]

params = {}

if CASTEP_TEMPLATE[-4:] != '.xml':
   # Read template parameter file
   try:
      params[0] = ('default', castep.CastepParam(CASTEP_TEMPLATE+'.param'))
   except IOError:
      die("Can't open parameter file %s" % CASTEP_TEMPLATE+'.param')
   except ValueError, message:
      die(str(message))

   # Now read the template cell file
   try:
      cell = castep.CastepCell(CASTEP_TEMPLATE+'.cell')
   except IOError:
      die("Can't open cell file %s" % CASTEP_TEMPLATE+'.cell')
   except ValueError, message:
      die(str(message))
else:
   # Read <castep_cell> and <castep_param> XML stanzas
   xml = xml.dom.minidom.parse(CASTEP_TEMPLATE)
   try:
      elements = xml.documentElement.getElementsByTagName('castep_param')
      if len(elements) == 0:
         raise ValueError('No <castep_param> element in XML file')

      for element in elements:
         if len(element.childNodes) != 1 or element.firstChild.nodeType != xml.TEXT_NODE:
            raise ValueError('<castep_param> element should have exactly one Text node')
      
         if not 'label' in element.attributes.keys():
            label = 'default'
         else:
            label = element.attributes['label'].value
         
         if not 'order' in element.attributes.keys():
            order = len(params)
         else:
            order = element.attributes['order'].value
            
         info('got parameter set %s, order %s' % (label, order))
         params[int(order)] = (label, castep.CastepParam(element.firstChild.data.split('\n')))

   except ValueError, message:
      die(str(message))

   # Now read the template cell file
   try:
      cell = castep.CastepCell(xml=xml)
   except ValueError, message:
      die(str(message))


# Read extended XYZ input file containing cluster
cluster = Atoms(xyzfile)

# remove old output file, if it's there
if os.path.exists(outfile):
   os.remove(outfile)

if DO_HASH:
   path = WORKING_DIR+'/'+stem+'_'+str(hash_atoms(cluster))
elif BATCH_QUEUE or BATCH_READ:
   try:
      batch_id = int(open('castep_driver_batch_id').read())
   except IOError:
      batch_id = 0
   path = '%s/%s-%08d' % (WORKING_DIR, stem, batch_id)
   batch_id += 1
   open('castep_driver_batch_id', 'w').write(str(batch_id))
else:
   path = WORKING_DIR+'/'+stem

# Make working directory if necessary
if not os.path.isdir(path):
   os.mkdir(path)
os.chdir(path)


if not BATCH_READ:
   # Load up old cluster, if it's there
   if os.path.exists(stem+'.xyz.old'):
      info('found old cluster in file %s' % stem+'.xyz.old')
      try:
         old_cluster = Atoms(stem+'.xyz.old')
      except IOError:
         die('error opening old cluster file %s' % stem+'.xyz.old')

      if (all(cluster.lattice == old_cluster.lattice)):
         info('lattice matches that of previous cluster')
      else:
         warn('lattice mismatch with prevous cluster')

      if (cluster.n == old_cluster.n):
         info('RMS position difference is %.3f A' % rms_diff(cluster, old_cluster))
      else:
         warn('number mismatch with previous cluster')

   else:
      old_cluster = cluster

   # Write lattice and atomic positions into cell object
   cell.update_from_atoms(cluster)

   # Check pseudopotentials are present and
   # copy pseudopotential files into working directory
   try:
      castep_check_pspots(cluster, cell, params[0][1])
   except IOError, (errno, msg):
      die('IO Error reading pseudopotentials (%d,%s)' % (errno,msg))
   except ValueError, message:
      die(str(message))

   # If we only got one parameter set use old mechanism to fork it into
   # multiple parameter sets - first DM, then EDFT, them DM without reuse,
   # then EDFT without reuse
   if len(params) == 1 and params[0][0] == 'default':
      info('Only default paramater set found: forking it to produce')
      info(' "standard", "without DM", "without reuse" and "without DM or reuse".')
      standard_params = params[0][1].copy()

      # Fall back to EDFT if DM fails to converge
      if standard_params.has_key('metals_method') and \
             standard_params['metals_method'].lower() == 'dm(edft)':
         params[0][1]['metals_method'] = 'dm'
         if 'max_scf_cycles_dm' in standard_params:
            params[0][1]['max_scf_cycles'] = standard_params['max_scf_cycles_dm']
         params[1] = ('without DM', standard_params.copy())
         params[1][1]['metals_method'] = 'edft'
         if 'max_scf_cycles_edft' in standard_params:
            params[1][1]['max_scf_cycles'] = standard_params['max_scf_cycles_edft']

      # Don't reuse checkpoint file
      params[2] = ('without reuse', standard_params.copy())
      params[2][1]['reuse'] = 'NULL'

      # With neither DM nor checkfile reuse
      if (standard_params.has_key('metals_method') and
          standard_params['metals_method'].lower() == 'dm(edft)'):
         params[3] = ('without DM or reuse', standard_params.copy())
         params[3][1]['reuse'] = 'NULL'
         params[3][1]['metals_method'] = 'edft'
         if 'max_scf_cycles_edft' in standard_params:
            params[3][1]['max_scf_cycles'] = standard_params['max_scf_cycles_edft']


# Main CASTEP invocation loop.
# Loop goes round as many times as there are parameter sets
sorted_param_keys = params.keys()[:]
sorted_param_keys.sort()
sorted_param_keys.reverse()
try:
   while sorted_param_keys:

      order = sorted_param_keys.pop()
      name, param = params[order]
      info('running CASTEP with parameter set <%s>' % name)

      # Ensure iprint is >= 2 since CASTEP doesn't output parameters properly for iprint=1
      if 'iprint' in param:
         if int(param['iprint']) < 2:
            param['iprint'] = '2'
      else:
         param['iprint'] = 2

      # Try to reuse check file if it's there and lattice matches
      # that of previous cluster
      if USE_CHECK_FILES and (('reuse' in param and param['reuse'].upper() != 'NULL') or not 'reuse' in param):
         if (os.path.exists(stem+'.check') and
             cluster.n == old_cluster.n):
#             all(cluster.lattice == old_cluster.lattice)):
            info('check file found: trying to reuse it')
            param['reuse'] = 'default'
         else:
            if cluster.n != old_cluster.n:
               info('cannot reuse check file - atom number mismatch')
            elif not all(cluster.lattice == old_cluster.lattice):
               info('cannot reuse check file - lattice mismatch')
            else:
               info('check file not found')
               
            # Discard parameter sets in order until we find one where
            # we don't try to reuse the electronic density.
            while True:
               try:
                  order = sorted_param_keys.pop()
               except IndexError:
                  raise ParamError
               name, param = params[order]
               if param.get('reuse', 'default').upper() != 'NULL': continue
               sorted_param_keys.append(order) # found one with no reuse, put it back
               break
            continue
      else:
         info('not reusing check file')
         param['reuse'] = 'NULL'

      if 'max_scf_cycles_dm'   in param: del param['max_scf_cycles_dm']
      if 'max_scf_cycles_edft' in param: del param['max_scf_cycles_edft']

      # Run castep
      if not BATCH_READ and not BATCH_QUEUE:
         if not castep_run(cell, param, stem, CASTEP, log=stem+'.castep_log'):
            error('castep run failed')
            continue

      if BATCH_QUEUE:
         param.write(stem+'.param')
         cell.write(stem+'.cell')
         
         info('batch queue mode: not running castep')
         cluster.params['energy'] = 0.0
         cluster.params['virial'] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
         cluster.add_property('force', 0.0, ncols=3)
         break
      else:
         try:
            cluster = castep.read_castep_output(stem+'.castep', cluster)
         except IOError, message:
            error('error parsing .castep file: %s' % message)
            continue
         except ValueError, message:
            error('error parsing .castep file: %s' % message)
            continue

         info('castep completed in %.1f s' % cluster.params['castep_run_time'])

         norm_f = norm(cluster.force)
         max_force_atom = norm_f.argmax()
         max_force = norm_f[max_force_atom]

         if hasattr(cluster,'hybrid'):
            info('max force is %.2f eV/A on atom %d hybrid=%r' % \
                 (max_force,max_force_atom,cluster.hybrid[max_force_atom]))
         else:
            info('max force is %.2f eV/A on atom %d' % (max_force,max_force_atom))

         if max_force > MAX_FORCE_THRESHOLD:
            error('max force is very large - repeating calculation %.6f %.6f' % (max_force, MAX_FORCE_THRESHOLD) )
            continue
         else:
            break
   else:
      raise ParamError

except ParamError:
   die('all parameter sets exhausted. Terminal error.')


# Save cluster for comparison with next time
old_cluster_file = open(stem+'.xyz.old', 'w')
cluster.write_xyz(old_cluster_file)
old_cluster_file.close()

# Finally change back to original working directory and write output file
os.chdir(orig_dir)
output_file = open(outfile, 'w')
cluster.write_xyz(output_file)
output_file.close()

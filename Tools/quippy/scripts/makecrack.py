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
from quippy.crack import *
from numpy import *
import sys
         
if __name__ == '__main__':

   try:
      params = CrackParams()

      if len(sys.argv[1:]) < 1:
         print('Usage: makecrack [-V] <stem>')
         print('Reads parameter file <stem>.xml and writes NetCDF output file <stem>.nc')
         print('If -V option is given then XML file is validated against DTD.')
         print('')
         print('Available parameters and their default values are:')
         params.print_()

      validate = False
      if sys.argv[1] == '-V':
         validate = True
         del sys.argv[1]

      stem = sys.argv[1]
      xmlfilename = stem+'.xml'

      print('Reading parameters from file %s with XML validation %s.' %
            (xmlfilename, {True:'enabled', False: 'disabled'}[validate]))

      xmlfile = InOutput(xmlfilename,INPUT)
      params.read_xml(xmlfile, validate=validate)
      xmlfile.close()

      crack_slab = makecrack(params, stem)
      if params.io_netcdf:
         crack_slab.write(stem+'.nc')
      else:
         crack_slab.write(stem+'.xyz')

   except RuntimeError, re:
      if is_interactive_shell():
         raise
      else:
         sys.stderr.write('error: %s\n' % str(re))
         sys.exit(1)

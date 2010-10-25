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
import optparse

p = optparse.OptionParser(usage='%prog [options] <init args>')

p.add_option('-l', '--label', action='store', help='Label name')
p.add_option('-o', '--outfile', action='store', help='Output file name (default stdout)', default='stdout')

opt, args = p.parse_args()

if len(args) == 0:
   p.error('No Potential init args given')

xml = """<params>

%s
</params>""" % quip_xml_parameters(' '.join(args), opt.label)

if opt.outfile == 'stdout':
    print xml
else:
    out = open(opt.outfile, 'w')
    out.write(xml)
    out.close()

#!/usr/bin/env python
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# H0 X
# H0 X   libAtoms+QUIP: atomistic simulation library
# H0 X
# H0 X   Portions of this code were written by
# H0 X     Albert Bartok-Partay, Silvia Cereda, Gabor Csanyi, James Kermode,
# H0 X     Ivan Solt, Wojciech Szlachta, Csilla Varnai, Steven Winfield.
# H0 X
# H0 X   Copyright 2006-2010.
# H0 X
# H0 X   These portions of the source code are released under the GNU General
# H0 X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# H0 X
# H0 X   If you would like to license the source code under different terms,
# H0 X   please contact Gabor Csanyi, gabor@csanyi.net
# H0 X
# H0 X   Portions of this code were written by Noam Bernstein as part of
# H0 X   his employment for the U.S. Government, and are not subject
# H0 X   to copyright in the USA.
# H0 X
# H0 X
# H0 X   When using this software, please cite the following reference:
# H0 X
# H0 X   http://www.libatoms.org
# H0 X
# H0 X  Additional contributions by
# H0 X    Alessio Comisso, Chiara Gattinoni, and Gianpietro Moras
# H0 X
# H0 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from xml.parsers.xmlproc import xmlproc
from xml.parsers.xmlproc import xmlval
from xml.parsers.xmlproc import xmldtd

def validate_xml(xml_filename, dtd_filename):
    """Validate a given XML file with a given external DTD.

    If the XML file is not valid, an error message will be printed
    to sys.stderr, and the program will be terminated with a non-zero
    exit code.  If the XML file is valid, nothing will be printed.
    """
    dtd = xmldtd.load_dtd(dtd_filename)
    parser = xmlproc.XMLProcessor()
    parser.set_application(xmlval.ValidatingApp(dtd, parser))
    parser.dtd = dtd
    parser.ent = dtd
    # If you want to override error handling, subclass
    # xml.parsers.xmlproc.xmlapp.ErrorHandler and do
    #   parser.set_error_handler(MyErrorHandler(parser))
    parser.parse_resource(xml_filename)
    # If you have xml data only as a string, use instead
    #   parser.feed(xml_data)
    #   parser.close()

if __name__ == '__main__':
    import sys
    if len(sys.argv) != 3:
        print 'Usage: validate.py <xml file> <dtd file>'
        sys.exit(0)
    validate_xml(sys.argv[1], sys.argv[2])
    print 'No errors found.'

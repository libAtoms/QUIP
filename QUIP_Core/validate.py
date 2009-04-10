#!/usr/bin/env python

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

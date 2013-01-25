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

"""Utility functions which will be imported into top-level quippy namespace"""

import sys, os, xml.dom.minidom
import glob
from dictmixin import PuPyDictionary

__all__ = ['infer_format', 'args_str', 'dict_to_args_str', 'parse_slice',
           'parse_comma_colon_list', 'loadstring',
           'quip_xml_parameters', 'is_interactive_shell',
           'parse_params', 'args_str_to_dict',
           'most_recent_file', 'time_ordered_glob', 'most_recent_files']

def infer_format(file, format, lookup):
    """Infer the correct format to read from or write to `file`

       If `format` is not None, there is nothing to do. Otherwise we
       look for file extension in `lookup`."""

    filename = None
    if format is None:
        if isinstance(file, basestring):
            filename = file
            if file in lookup:
                format = file
            else:
                file = os.path.expanduser(file)
                base = os.path.basename(file)
                base, ext = os.path.splitext(base)
                if ext == '':
                    format = base
                else:
                    format = ext[1:]
        else:
            format = file.__class__

    return filename, file, format

def args_str(D):
    """Construct args string from file, string or mapping object"""
    return str(PuPyDictionary(D))

dict_to_args_str = args_str # synonym for args_str() function

def parse_slice(S):
    """Parse string containing slice in form [start]:[stop]:[range] and return slice instance."""

    class SliceParser(object):
        def __getitem__(self, idx):
            return idx

    return eval('SliceParser()[%s]' % S)

def parse_params(s):
   """Read key=value pairs from a string or list of string and return a standard Python dictionary"""
   p = PuPyDictionary(s)
   return dict(p)

args_str_to_dict = parse_params # synonym for parse_params() function
   

def parse_comma_colon_list(L):
    """Parse a comma or colon seperated string into a list, converting each entry to lower-case."""
    if ':' in L:
        L = L.split(':')
    elif ',' in L:
        L = L.split(',')
    else:
        L = [L]

    return [k.lower() for k in L]

def loadstring(s):
    import StringIO
    from numpy import loadtxt
    return loadtxt(StringIO.StringIO(s.replace('[','').replace(']','')))

def quip_xml_parameters(name, label=None):
    """Search for an QUIP XML parameter set matching `name' and, optionally `label`.

    Result is the XML stanza as a string. Looks under
    ${QUIP_ROOT}/QUIP_Core/parameters for a file matching `name`,
    `name`.xml, ip.parms.`name`.xml or tightbind.params.`name`.xml. If
    `label` is given, it should be a text string matching the label
    attribute of a QUIP <params> XML stanza. If a matching parameter
    set can't be found, we raise an :exc:`IOError`.
    """

    from quippy import QUIP_ROOT
    from xml.parsers.expat import ExpatError
    from xml.dom.minidom import Element

    param_dir = os.path.join(QUIP_ROOT, 'QUIP_Core/parameters')

    if name.lower().startswith('filepot') or name.lower().startswith('callbackpot'):
        # FilePot and CallbackPot don't need XML params
        return ''

    if name.upper().startswith('IP ') or name.upper().startswith('TB '):
        name = name[3:]

    xml_files = [ os.path.join(param_dir, name),
                  os.path.join(param_dir, name) + '.xml',
                  os.path.join(param_dir, 'ip.parms.'+name) + '.xml',
                  os.path.join(param_dir, 'tightbind.parms.'+name) + '.xml' ]

    for xml_file in xml_files:
        if os.path.exists(xml_file): break
    else:
        raise IOError('Cannot find XML parameter file matching "%s"' % name)

    xml_string = open(xml_file).read()

    if label is not None:
        # Read XML and extract tag with matching label attribute
        try:
            d = xml.dom.minidom.parseString(xml_string)
        except ExpatError:
            # Add surrounding document element
            xml_string = r'<params>%s</params>' % xml_string
            d = xml.dom.minidom.parseString(xml_string).firstChild

        for child in d.childNodes:
            if not isinstance(child, Element): continue
            if child.hasAttribute('label'):
                if child.getAttribute('label') == label:
                    break
        else:
            raise IOError('Cannot find XML stanza in file "%s" matching label "%s"' % (xml_file, label))

        return str(child.toxml())
    else:
        # Return entire contents of file
        return str(xml_string)

def is_interactive_shell():
    """Return True if we're running in an interactive shell (including IPython)"""

    if hasattr(sys, 'ps1'):
        return True
    try:
        __IPYTHON__
    except NameError:
        return False
    else:
        return True


def read_text_file(fh):
    """
    Read all lines of text from a file

    fh can be a file, file-like object or filename. If `fh` is a
    filename, file is automatically opened and closed.

    Returns a tuple (filename, lines).
    """

    opened = False
    filename = fh
    if isinstance(fh, basestring):
        fh = open(fh, 'r')
        opened = True
    lines = [line for line in fh]
    if opened:
        fh.close()
    return filename, lines

def time_ordered_glob(pattern):
   """
   Return a list of files matching `pattern` sorted by modification time
   """
   
   return sorted(glob.glob(pattern), key=lambda f: os.stat(f).st_mtime)

def most_recent_file(pattern):
    """
    Find the most recent file matching glob `pattern`
    """

    return time_ordered_glob(pattern)[-1]


def most_recent_files(dir_pattern, file_pattern, suffix=''):
    """
    Find the most recent file matching file_pattern in each of the directories matching dir_pattern
    """
    return [most_recent_file(os.path.join(dir, file_pattern))+suffix for dir in glob.glob(dir_pattern)]


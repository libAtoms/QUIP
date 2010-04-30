# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# HP X
# HP X   pyatoms: atomistic simulations tools
# HP X
# HP X   Copyright James Kermode 2010
# HP X
# HP X   These portions of the source code are released under the GNU General
# HP X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
# HP X
# HP X   If you would like to license the source code under different terms,
# HP X   please contact James Kermode, james.kermode@gmail.com
# HP X
# HP X   When using this software, please cite the following reference:
# HP X
# HP X   http://www.jrkermode.co.uk/PyAtoms
# HP X
# HP XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

from ordereddict import OrderedDict

import string, re

class ParamReader(OrderedDict):
   """Subclass of OrderedDictionary for reading key/value pairs from strings or files.
      The original order of items is maintained. Values that looks like floats or ints
      or lists of floats or ints are automatically converted on reading."""
   
   def __init__(self, s=None, D=None):
      OrderedDict.__init__(self)
      if s is not None:
         self.read(s)
      if D is not None:
         self.update(D)

   def copy(self):
      return ParamReader(D=OrderedDict.copy(self))

   def __copy__(self):
      return self.copy()
         
   def parse(self, s):
        key_quoted_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*"([^"]+)"\s*')
        key_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*([-0-9A-Za-z_.:\[\]()]+)\s*')

        s = s.strip()

        while 1:
           # Match quoted string first, then fall through to plain key=value
           m = key_quoted_value.match(s)
           if m is None:
              m = key_value.match(s)
              if m is not None:
                 s = key_value.sub('', s, 1)
           else:
              s = key_quoted_value.sub('', s, 1)

           if m is None: break # No more matches

           key = m.group(1)
           value = m.group(2)

           # Try to convert to (list of) floats or ints
           try:
              numvalue = []
              for x in string.split(value):
                 if x.find('.') == -1:
                    numvalue.append(int(float(x)))
                 else:
                    numvalue.append(float(x))
              if (len(numvalue) == 1): 
                 numvalue = numvalue[0] # Only one number
              value = numvalue
           except ValueError:
              pass

           self[key] = value

   def read(self, f):
      if type(f) == type(''):
         self.parse(f)
      else:
         try:
            for line in f:
               self.parse(line)
         except TypeError:
            raise TypeError("Don't know how to read from object - "+str(f))


   def __repr__(self):
      return "ParamReader('"+str(self)+"')"


   def __str__(self):

      if len(self) == 0: return ''
      
      s = ''
      for key in self.keys():
         val = self[key]
         if type(val) == type([]):
            val = ' '.join(map(str, val))

         if type(val) == type('') and ' ' in val:
            s = s + '%s="%s" ' % (key, val)
         else:
            s = s + '%s=%s ' % (key, str(val))

      return s.strip()


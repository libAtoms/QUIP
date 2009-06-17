from ordereddict import OrderedDict

import string, re

class ParamReader(OrderedDict):
   """Subclass of OrderedDictionary for reading key/value pairs from strings or files.
      The original order of items is maintained. Values that looks like floats or ints
      or lists of floats or ints are automatically converted on reading."""
   
   def __init__(self, source=None):
      OrderedDict.__init__(self)
      if source is not None:
         self.read(source)

   def copy(self):
      return ParamReader(OrderedDict.copy(self))

   def __copy__(self):
      return self.copy()
         
   def parse(self, s):
        key_quoted_value = re.compile(r'([A-Za-z_]+[A-Za-z0-9_]*)\s*=\s*["\{\}]([^"\{\}]+)["\{\}]\s*')
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
      if isinstance(f, str):
         self.parse(f)
      elif isinstance(f, dict):
         self.update(f)
      else:
         try:
            for line in f:
               self.parse(line)
         except TypeError:
            raise TypeError("Don't know how to read from object - "+str(f))


   def __repr__(self):
      return "%s('%s')" % (self.__class__.__name__, str(self))


   def __str__(self):
      return self.asstring()

   def asstring(self, sep=' '):
      if len(self) == 0: return ''

      type_val_map = {(bool,True): 'T', (bool,False): 'F'}
      
      s = ''
      for key in self.keys():
         val = self[key]

         if hasattr(val, '__iter__'):
            val = ' '.join([ str(type_val_map.get((type(x),x), x)) for x in val ])
         else:
            val = type_val_map.get((type(val),val), val)

         if type(val) == type('') and ' ' in val:
            s = s + '%s="%s"%s' % (key, val, sep)
         else:
            s = s + '%s=%s%s' % (key, str(val), sep)

      return s.strip()


   def write(self, f):
      f.write(self.asstring(sep='\n'))


def args_str(D):
   """Construct args string from file, string or mapping object"""
   return str(ParamReader(D))

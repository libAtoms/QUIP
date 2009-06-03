def args_str(d):
   """Construct args string from dictionary-like object"""

   val_map = {True: 'T', False: 'F'}

   if len(d) == 0: return ''

   s = ''
   for key in d.keys():
      val = d[key]

      val = val_map.get(val, val)
      
      if hasattr(val, '__iter__'):
         val = ' '.join(map(str, val))

      if type(val) == type('') and ' ' in val:
         s = s + '%s="%s" ' % (key, val)
      else:
         s = s + '%s=%s ' % (key, str(val))

   return s.strip()


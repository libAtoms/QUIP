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

from quippy import FortranPotential

class Potential(FortranPotential):

   __doc__ = FortranPotential.__doc__
   callback_map = {}
   
   def __init__(self, args_str, *args, **kwargs):
      if args_str.lower().startswith('callbackpot') and not 'label' in args_str:
         args_str = args_str + ' label=%d' % id(self)
      FortranPotential.__init__(self, args_str, *args, **kwargs)

      if args_str.lower().startswith('callbackpot'):
         FortranPotential.set_callback(self, Potential.callback)

   __init__.__doc__ = FortranPotential.__init__.__doc__

   @staticmethod
   def callback(at_ptr):
      from quippy import Atoms
      at = Atoms(fpointer=at_ptr, finalise=False)
      if at.params['label'] not in Potential.callback_map:
         raise ValueError('Unknown Callback potential label %s' % at.params['label'])
      Potential.callback_map[at.params['label']](at)

   def set_callback(self, callback):
      Potential.callback_map[str(id(self))] = callback



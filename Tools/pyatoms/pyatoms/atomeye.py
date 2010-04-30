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

import sys, re, os, os.path, pexpect, atexit, string, time, thread
from threading import Thread

from numpy import *
from atoms import *

# Location of atomeye executable
atomeye_path = 'A'

# Thread to munge the atomeye output pipe
class AtomEyeReader(Thread):
   def __init__(self, myatomeye):
      Thread.__init__(self)
      self.setDaemon(True) # Daemonise so we don't stop python closing
      self.atomeye = myatomeye
      self.lines = []
      
   def run(self):

      # Reg exps for things we want to echo to console
      atom_clicked = re.compile(r'^([0-9]*) atom ')
      atom_distance = re.compile(r'^x\[[0-9]*\]')
      file_change = re.compile(r'^Loading configuration from "(.*)"')
      load_finished = re.compile('^All bin-related allocations freed')

      while 1:
         
         # Check child is still there
         if not self.atomeye.child.isalive():
            break

         try:
            line = self.atomeye.readline()
         except pexpect.EOF:
            break
         except pexpect.TIMEOUT:
            continue

         # Save all the output
         self.lines.append(line.strip())

         m = atom_clicked.search(line)
         if not m is None:
            self.atomeye.on_click_atom(int(m.group(1)))

         if atom_distance.search(line):
            print line.strip()

         # Frame we're displaying has changed. Update and call hook function
         m = file_change.search(line)
         if not m is None:
            self.atomeye.on_change_frame(m.group(1))

         if load_finished.search(line) and self.atomeye.busy_loading.locked():
            self.atomeye.busy_loading.release()


   # Print everything out
   def printlog(self):
      for line in self.lines:
         print line

         

# Class to launch and communicate with an AtomEye process
class AtomEye(object):
   def __init__(self, atoms=None, *showargs, **kwshowargs):
      self.child = None
      self.wipe()
      self.start()
      atexit.register(self.close)
      if atoms is not None:
         self.show(atoms, *showargs, **kwshowargs)
      

   def wipe(self):
      self._cleanup_files()
      
      self._atoms = None

      self.verbose = False

      self.selecting = False
      self.selection = []
      
      self.atoms_cache = []
      self.atoms_refs = []
      self.cfg_files = []
      
      self.aux_props = None
      self.saved_view = None
      self.cfgbasename = os.tmpnam()
      self.viewfile = '%s.view' % self.cfgbasename


   def start(self):
      self.close()
      self.child = pexpect.spawn(atomeye_path, timeout=5)
      self.child.delaybeforesend = 0.0
      self.busy_loading = thread.allocate_lock()
      self.busy_loading.acquire()
      self.readline = self.child.readline
      self.reader = AtomEyeReader(self)
      self.printlog = self.reader.printlog
      self.reader.start()
      self('toggle_parallel_projection')


   def __call__(self, s):
      self.child.sendline(s)


   def __del__(self):
      self.close()

   def close(self):
      if self.child is not None and self.child.isalive():
         self('quit')
         self.reader.join()

      self._cleanup_files()


   def _cleanup_files(self):
      if not hasattr(self,'cfg_files'): return
      for f in self.cfg_files:
         if os.path.exists(f):
            os.unlink(f)

      if os.path.exists(self.viewfile):
         os.unlink(self.viewfile)


   def on_change_frame(self, cfgfile):
      "Print frame specific parameters"

      try:
         self._atoms = self.atoms_refs[self.cfg_files.index(cfgfile)]
      except IndexError:
         raise ValueError('Unknown cfg file %s' % cfgfile)

      if self.verbose:
         print
         self.status()
         print


   def on_click_atom(self, i):
      if self._atoms is None: return

      # Take care of AtomEye cell doubling
      if i >= self._atoms.n: i = mod(i,self._atoms.n)

      if self.selecting:
         self.selection.append(i)

      # Print atom information and properties
      if self.verbose:
         print
         print 'Atom %d (%s)' % (i, self._atoms.species[i])
         for key in self._atoms.properties:
            print '%-15s = ' % (key) + str(getattr(self._atoms,key)[i])
         print
      

   def select(self):
      "Return a list of atoms, selected by clicking on them in AtomEye view window"
      self.selection = []
      self.selecting = True
      raw_input('Right click on some atoms. Press enter to finish...')
      self.selecting = False
      return self.selection
   

   def highlight(self, atomlist):
      if not 'highlight' in self._atoms.properties:
         self._atoms.add_property('highlight',0)

      self._atoms.highlight[:] = 0
      self._atoms.highlight[atomlist] = 1

      self.show(self._atoms, aux_prop='highlight')


   def animate(self, atomseq, *showargs, **kwshowargs):
      for atoms in atomseq:
         self.show(atoms, save_view=False, *showargs, **kwshowargs)


   def status(self):
      if self._atoms is None:
         print 'No atoms loaded'
      else:
         print '%d Atoms' % self._atoms.n
         for key in self._atoms.params:
            print '%-15s = %s' % (key, str(self._atoms.params[key]))


   def show(self, atoms=None, aux_prop=None, geometry=None,
            width=None,height=None,revert_view=False,save_view=True,
            *alignargs,**alignkwargs):

      if revert_view: self.revert_view()
      if geometry is not None: self.set_geometry(geometry)
      if atoms is not None: self.set_atoms(atoms)
      if aux_prop is not None: self.set_aux_prop(aux_prop)

      if width is not None and height is not None:
         raise ValueError("Can't set both width and height at same time")
      if width  is not None:  self.set_width(width)
      if height is not None:  self.set_height(height)
      if save_view: self.save_view()

      # Pass remaining args onto align
      if len(alignargs) != 0 or len(alignkwargs) != 0:
         self.align(*alignargs, **alignkwargs)


   def set_atoms(self, atoms):
      restarted = False
      if self.child is None or not self.child.isalive():
         self.start()
         restarted = True

      # Wait if AtomEye is currently loading a file
      while self.busy_loading.locked(): pass

      new_atoms = False
      if atoms is not None:
         new_atoms = not self._atoms is atoms
         self._atoms = atoms

      # Make a new cache entry if this an atoms object we haven't seen before
      id_list = map(id, self.atoms_refs)
      if not id(self._atoms) in id_list:
         self.atoms_cache.append(None)
         self.atoms_refs.append(None)
         self.cfg_files.append('%s_%08d.cfg' % (self.cfgbasename, len(self.cfg_files)))
         id_list.append(id(self._atoms))
         
      # Find index of self.atoms in atoms_ref list. cfg_files and atoms_cache
      # have same ordering, so index is valid for them too.
      index = id_list.index(id(self.atoms))
      
      # Has atoms changed since we last wrote a .cfg for this atom id?
      rewritten = False
      if self._atoms != self.atoms_cache[index] or not os.path.exists(self.cfg_files[index]):
         print 'writing cfg...'
         self.aux_props = self._atoms.write_cfg(self.cfg_files[index])[5:]
         self.atoms_refs[index] = self._atoms
         self.atoms_cache[index] = self._atoms.copy()
         rewritten = True

      if rewritten or restarted or new_atoms:
         self.busy_loading.acquire()
         self('load_config %s' % self.cfg_files[index])

   def get_atoms(self):
      return self._atoms

   atoms = property(get_atoms, set_atoms)

   def set_aux_prop(self, aux_prop):
      try:
         self('aux_property_coloring %d' % (self.aux_props.index(aux_prop)))
      except ValueError:
         raise ValueError("Can't find aux_prop %s in Atoms object" % aux_prop)

   def get_aux_prop(self):
      if self.aux_props is None: return None
      view = self.get_view()
      proplines = filter(lambda x: x.startswith('set n->auxiliary_idx'), view)
      if len(proplines) != 1:
         raise ValueError('Bad view - missing auiliary_idx line %r' % proplines)

      dum1, dum2, aux_idx = proplines[0].split()
      return self.aux_props[int(aux_idx)]

   aux_prop = property(get_aux_prop, set_aux_prop)


   def align(self, align=None, catom=None, cpos=None, cfrac_pos=None, cbond=None):
      if catom is not None:
         cpos = self._atoms.frac_pos[catom,:] + array([0.5,0.5,0.5])
      elif cpos is not None:
         cpos = dot(cpos, self._atoms.g) + array([0.5,0.5,0.5])
      elif cfrac_pos is not None:
         cpos = cfrac_pos + array([0.5,0.5,0.5])
      elif align is not None:
         if align == 'CoM':
            cpos = (array([dot(self._atoms.frac_pos[:,0],self._atoms.mass),
                           dot(self._atoms.frac_pos[:,1],self._atoms.mass),
                           dot(self._atoms.frac_pos[:,2],self._atoms.mass)])
                    /sum(self._atoms.mass)) + array([0.5,0.5,0.5])
         elif align == 'CrackPos':
            if not 'CrackPos' in self._atoms.params:
               raise ValueError("No CrackPos found in Atoms object")
            cpos = array([self._atoms.params['CrackPos'],0.0,0.0])
            cpos = dot(cpos, self._atoms.g) + array([0.5,0.5,0.5])
         elif align == 'center':
            cpos = array([0.5,0.5,0.5])
         else:
            raise ValueError("Unknown alignment option %s" % align)
      elif cbond is not None:
         cpos = 0.5*(self._atoms.frac_pos[cbond[0]]+ self._atoms.frac_pos[cbond[1]])
      else:
         cpos = array([0.5,0.5,0.5]) # default is centre of cell

      self('xtal_origin_goto %f %f %f' % tuple(cpos))
      

   def get_view(self):
      if os.path.exists(self.viewfile):
         os.unlink(self.viewfile)
      self('save %s' % self.viewfile)
      while not os.path.exists(self.viewfile):
         time.sleep(0.1)
      view = open(self.viewfile).readlines()
      return view


   def set_view(self, view):
      open(self.viewfile,'w').writelines(view)
      self('load_script %s' % self.viewfile)
      self('shift_xtal 0 0')

   view = property(get_view, set_view)


   def save_view(self):
      "Save this view so that self.revert_view() will bring us back to it"
      self.saved_view = self.get_view()


   def revert_view(self):
      "Revert viewport to last saved view"
      if self.saved_view is not None:
         self.set_view(self.saved_view)


   def _get_size_and_geometry(self, view=None):
      
      if view is None:
         view = self.get_view()

      geometryline = filter(lambda x: x.startswith('resize'),view)
      if len(geometryline) != 1:
         raise ValueError('Bad view - resize is mssing %r' % geometryline)
      dummy, wx, wy = geometryline[0].split()
      geometry = map(float, (wx,wy))
      
      xline = filter(lambda x: x.startswith('set AX_3D->x'),view)
      if len(xline) != 1:
         raise ValueError('Bad view - AX3D-> x is missing %r' % xline)
      dum1, dum2, x, y, z = xline[0].split()
      x, y, z = map(float, (x, y, z))
      
      width = (2.49915-z)/746.942 # semi-empirical fit
      height = width/geometry[0]*geometry[1]

      return width, height, geometry

   def set_width(self, width):
      cur_width = self.get_width()
      if width > cur_width:
         delta = 1.0 - float(width)/float(cur_width)
      else:
         delta = float(cur_width)/float(width) - 1.0
      self('advance %f' % delta)

   def get_width(self):
      width, height, geometry = self._get_size_and_geometry()
      return width

   width = property(get_width, set_width)

   def set_height(self, height):
      cur_height = self.get_height()
      if height > cur_height:
         delta = 1.0 - float(height)/float(cur_height)
      else:
         delta = float(cur_height)/float(height) - 1.0
      self('advance %f' % delta)

   def get_height(self):
      width, height, geometry = self._get_size_and_geometry()
      return height      

   height = property(get_height, set_height)

   def set_geometry(self, geometry):
      self('resize %d %d' % geometry)

   def get_geometry(self):
      width, height, geometry = self._get_size_and_geometry()
      return geometry

   geometry = property(get_geometry, set_geometry)

   

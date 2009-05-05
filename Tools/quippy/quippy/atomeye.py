"""High-level interface to atomeye. Low-level interface is in _atomeye extension module."""

import threading, _atomeye, time, sys
from math import ceil, log10

from quippy import CInOutput
__cio = CInOutput()

__window_id = 0 # for now, we always have only one window
__thread = None

atoms = None
__olddata = None

frame = 0
delta = 1

paint_property = None
paint_value = 1

def on_atom_click(idx):
    if atoms is None: return
    theat = atoms
    if hasattr(atoms, '__iter__'):
        theat = atoms[frame]
    
    if idx > theat.n:
        idx = idx % theat.n
    print "frame %d, atom %d clicked" % (frame, idx)
    
    if paint_property is not None and theat.has_property(paint_property):
        getattr(theat, paint_property)[idx] = paint_value
        redraw()

def on_advance(mode):
    global atoms, frame, delta
    
    if not hasattr(atoms,'__iter__'): return
    if mode == 'forward':
        frame += delta
        if frame >= len(atoms):
            frame = 0
    elif mode == 'backward':
        frame -= delta
        if frame < 0:
            frame = len(atoms)-1
    elif mode == 'first':
        frame = 0
    elif mode == 'last':
        frame = len(atoms)-1
    redraw()
    

def paint(property='selection',value=1,fill=0):
    global paint_property, paint_value
    if atoms is None: return
    if not atoms.has_property(property):
        atoms.add_property(property, fill)
    paint_property = property
    paint_value = value
    _atomeye.load_libatoms(__window_id, __cio.update(atoms), 'paint')
    aux_property_coloring(paint_property)

def start():
    global __thread
    if isAlive(): return
    __thread = threading.Thread(target=_atomeye.start, args=(on_atom_click,on_advance))
    __thread.setDaemon(True)
    __thread.start()

    # wait for AtomEye to be initialised succesfully
    while not _atomeye.isAlive():
        time.sleep(0.1)

    for funcname in dir(sys.modules[__name__]):
        h = help(funcname)
        if not 'unknown command' in h:
            getattr(sys.modules[__name__],funcname).__doc__ = h


def show(obj, property=None, frame=None, animate=False):
    global atoms

    if not isAlive(): start()

    atoms = obj
    if hasattr(obj,'__iter__'):
        if frame is not None:
            if frame < 0: frame = len(atoms)-frame
            if frame >= len(atoms): frame=len(atoms)-1
            setattr(sys.modules[__name__],'frame',frame)
        else:
            frame = getattr(sys.modules[__name__],'frame')
        if animate:
            while True:
                on_advance('forward')
        else:
            redraw(property=property)
    else:
        redraw(property=property)

    
def isAlive():
    return _atomeye.isAlive()

def redraw(property=None):
    global __olddata

    if not isAlive(): 
        raise RuntimeError('AtomEye not running')
    
    if atoms is None:
        raise RuntimeError('Nothing to view -- set atomeye.atoms to Atoms or sequence of Atoms')

    theat = atoms
    if hasattr(atoms, '__iter__'):
        theat = atoms[frame]
        fmt = "%%0%dd" % ceil(log10(len(atoms)+1))
        title = 'atoms[%s/%s]' % (fmt % frame, fmt % len(atoms))
    else:
        title = 'atoms'

    #if __olddata is None or __oldproperty != property or not theat.data.equal(__olddata):
    #    __olddata = theat.data.copy()
    _atomeye.load_libatoms(__window_id, __cio.update(theat), title)
    if property is not None:
        aux_property_coloring(property)
    _atomeye.redraw(__window_id)

def run_command(command, expect_output=False):
    if not isAlive(): 
        raise RuntimeError('AtomEye not running')
    res = _atomeye.run_command(__window_id, command)
    if not expect_output:
        if res is not None:
            raise RuntimeError(res)
    else:
        return res
    

def help(command):
    return _atomeye.help(__window_id, command)

def close():
    if not isAlive(): 
        raise RuntimeError('AtomEye not running')
    _atomeye.close(__window_id)

def set(key, value):
    res = _atomeye.run_command(__window_id, "set %s %s" % (str(key), str(value)))
    if res is not None:
        raise ValueError(res)

def save(filename):
    run_command("save %s" % str(filename))

def load_script(arg):
    run_command("load_script %s" % str(arg))

def key(key):
    run_command("key %s" % key)

def toggle_coordination_coloring():
    run_command("toggle_coordination_coloring")

def translate(axis, delta):
    run_command("translate %d %f " % (axis, delta))

def shift_xtal(axis, delta):
    run_command("shift_xtal %d %f" % (axis, delta))

def rotate(axis, theta):
    run_command("rotate %d %f" % (axis, theta))

def advance(delta):
    run_command("advance %f" % delta)

def shift_cutting_plane(delta):
    run_command("shift_cutting_plane %f" % delta)

def change_bgcolor(color):
    run_command("change_bgcolor %f %f %f" % (color[0], color[1], color[2]))

def change_atom_r_ratio(delta):
    run_command("change_atom_r_ratio %f" % delta)

def change_bond_radius(delta):
    run_command("change_bond_radius %f" % delta)

def change_view_angle_amplification(delta):
    run_command("change_view_angle_amplification %f" % delta)

def toggle_parallel_projection():
    run_command("toggle_parallel_projection")

def toggle_bond_mode():
    run_command("toggle_bond_mode" )

def normal_coloring():
    run_command("normal_coloring")

def aux_property_coloring(auxprop):
    if isinstance(auxprop,int): auxprop = str(auxprop)
    run_command("aux_property_coloring %s" % auxprop)

def central_symmetry_coloring():
    run_command("central_symmetry_coloring")

def change_aux_property_threshold(lower_upper, delta):
    if isinstance(lower_upper, int): lower_upper = str(lower_upper)
    run_command("change_aux_property_threshold %s %f" % (lower_upper, delta))

def reset_aux_property_thresholds():
    run_command("reset_aux_property_thresholds")

def toggle_aux_property_thresholds_saturation():
    run_command("toggle_aux_property_thresholds_saturation")

def toggle_aux_property_thresholds_rigid():
    run_command("toggle_aux_property_thresholds_rigid")

def rcut_patch(sym1, sym2, inc_dec, delta=None):
    run_command("rcut_patch start %s %s" % (sym1,sym2))
    if delta is None:
        run_command("rcut_patch %s" % inc_dec)
    else:
        run_command("rcut_patch %s %f" % (inc_dec, delta))
    run_command("rcut_patch finish")

def select_gear(gear):
    run_command("select_gear %d" % gear)

def cutting_plane(n, d, s):
    run_command("cutting_plane %d %f %f %f %f %f %f" % \
                             (n, d[0], d[1], d[2], s[0], s[1], s[2]))

def shift_cutting_plane_to_anchor(n):
    run_command("shift_cutting_plane_to_anchor %d" % n)

def delete_cutting_plane(n):
    run_command("delete_cutting_plane %d" % n)

def flip_cutting_plane(n):
    run_command("flip_cutting_plane %d" % n)

def capture(filename, resolution=None):
    if resolution is None: resolution = ""
    format = filename[filename.rindex('.')+1:]
    run_command("capture %s %s %s" % (format, filename, resolution))

def change_wireframe_mode():
    run_command("change_wireframe_mode")

def change_cutting_plane_wireframe_mode():
    run_command("change_cutting_plane_wireframe_mode")

def load_config(filename):
    run_command("load_config %s" % filename)

def load_config_advance(command):
    run_command("load_config_advance %s" % command)

def script_animate(filename):
    run_command("script_animate %s" % filename)

def load_atom_color(filename):
    run_command("load_atom_color %s" % filename)

def load_aux(filename):
    run_command("load_aux %s" % filename)

def look_at_the_anchor():
    run_command("look_at_the_anchor")

def observer_goto():
    run_command("observer_goto")

def xtal_origin_goto(s):
    run_command("xtal_origin_goto %f %f %f" % (s[0], s[1], s[2]))

def find_atom(i):
    run_command("find_atom %d" % i)

def resize(width, height):
    run_command("resize %d %d" % (width, height))

def change_aux_colormap(n):
    run_command("change_aux_colormap %d" % n)

def print_atom_info(i):
    run_command("print_atom_info %d" % i)

def save_atom_indices():
    run_command("save_atom_indices")

def change_central_symm_neighbormax():
    run_command("change_central_symm_neighbormax")

def timer(label):
    run_command("timer %s" % label)

def isoatomic_reference_imprint():
    run_command("isoatomic_reference_imprint")

def toggle_shell_viewer_mode():
    run_command("toggle_shell_viewer_mode")

def toggle_xtal_mode():
    run_command("toggle_xtal_mode")

def change_shear_strain_subtract_mean():
    run_command("change_shear_strain_subtract_mean")

def zoom_to_fit():
    pass

// HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
// HQ X
// HQ X   quippy: Python interface to QUIP atomistic simulation library
// HQ X
// HQ X   Copyright James Kermode 2010
// HQ X
// HQ X   These portions of the source code are released under the GNU General
// HQ X   Public License, version 2, http://www.gnu.org/copyleft/gpl.html
// HQ X
// HQ X   If you would like to license the source code under different terms,
// HQ X   please contact James Kermode, james.kermode@gmail.com
// HQ X
// HQ X   When using this software, please cite the following reference:
// HQ X
// HQ X   http://www.jrkermode.co.uk/quippy
// HQ X
// HQ XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

#include <Python.h>
#include <pthread.h>
#include <signal.h>

#include "numpy/arrayobject.h"

#include <xyz_netcdf.h>
#include <atomeyelib.h>

static char atomeye_doc[] = 
"This module interfaces to AtomEye.";

static PyObject *on_click_atom_pyfunc = NULL;
static PyObject *on_advance_pyfunc = NULL;
static PyObject *on_close_pyfunc = NULL;
static PyObject *on_new_pyfunc = NULL;

static Atoms *atomeye_atoms = NULL;

static void on_click_atom(int iw, int atom)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_click_atom_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,i)", iw, atom);       
    PyEval_CallObject(on_click_atom_pyfunc, arglist); 
    Py_DECREF(arglist);                               
  }
  PyGILState_Release(state);
}

static void on_close(int iw)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_close_pyfunc != NULL) {
    arglist = Py_BuildValue("(i)", iw);
    PyEval_CallObject(on_close_pyfunc, arglist);
    Py_DECREF(arglist);
  }
  PyGILState_Release(state);
}

static void on_advance(int iw, char *instr)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_advance_pyfunc != NULL) {
    arglist = Py_BuildValue("(i,s)", iw, instr);
    PyEval_CallObject(on_advance_pyfunc, arglist);
    Py_DECREF(arglist);                           
  }
  PyGILState_Release(state);
}

static void on_new(int iw)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  if (on_new_pyfunc != NULL) {
    arglist = Py_BuildValue("(i)", iw);       
    PyEval_CallObject(on_new_pyfunc, arglist);
    Py_DECREF(arglist);                       
  }
  PyGILState_Release(state);
}

static int update_atoms_structure(PyObject *pyat)
{
  int i,j,lookup[3];
  PyObject *n = NULL, *data = NULL, *lattice = NULL, *properties = NULL, *items = NULL, 
    *item = NULL, *key = NULL, *value = NULL,
    *intsize = NULL, *realsize = NULL, *strsize = NULL, *logicalsize = NULL,
    *p = NULL, *valuetuple = NULL,
    *intarray = NULL, *realarray = NULL, *strarray = NULL, *logicalarray = NULL;

  atomeye_atoms->n_param = 0;

  /* atoms.data - quippy.Table or compatible type */
  if ((data = PyObject_GetAttrString(pyat, "data")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data");
    goto fail;    
  }

  /* atoms.n - int */
  if ((n = PyObject_GetAttrString(pyat, "n")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.n");
    goto fail;
  }
  if (!PyInt_Check(n)) {
    PyErr_SetString(PyExc_TypeError, "atoms.n must be an integer");
    goto fail;
  }
  atomeye_atoms->n_atom = (size_t)PyInt_AsLong(n);

  /* atoms.data.intsize - int */
  if ((intsize = PyObject_GetAttrString(data, "intsize")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.intsize");
    goto fail;
  }
  atomeye_atoms->n_int = (int)PyInt_AsLong(intsize);
  
  /* atoms.data.realsize - int */
  if ((realsize = PyObject_GetAttrString(data, "realsize")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.realsize");
    goto fail;
  }
  atomeye_atoms->n_real = (int)PyInt_AsLong(realsize);

  /* atoms.data.strsize - int */
  if ((strsize = PyObject_GetAttrString(data, "strsize")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.strsize");
    goto fail;
  }
  atomeye_atoms->n_str = (int)PyInt_AsLong(strsize);

  /* atoms.data.logicalsize - int */
  if ((logicalsize = PyObject_GetAttrString(data, "logicalsize")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.logicalsize");
    goto fail;
  }
  atomeye_atoms->n_logical = (int)PyInt_AsLong(logicalsize);

  /* atoms.properties  - quippy.Dictionary or compatible type */
  if ((properties = PyObject_GetAttrString(pyat, "properties")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.properties");
    goto fail;
  }
  if (!PyMapping_Check(properties)) {
    PyErr_SetString(PyExc_TypeError, "atoms.properties must be a mapping type");
    goto fail;
  }
  atomeye_atoms->n_property = (int)PyMapping_Size(properties);

  if ((items = PyMapping_Items(properties)) == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot get items from atoms.properties");
    goto fail;
  }

  for (i=0; i< atomeye_atoms->n_property; i++) {
    if ((item = PySequence_GetItem(items, i)) == NULL) {
      PyErr_SetString(PyExc_ValueError, "Error getting item from atoms.properties.items()");
      goto fail;
    }
    key = PySequence_GetItem(item, 0);
    value = PySequence_GetItem(item, 1);
    valuetuple = PySequence_Tuple(value);
    
    if(!PyString_Check(key)) {
      PyErr_SetString(PyExc_TypeError, "Keys in atoms.proeprtiees.keys() must be of string type");
      goto fail;
    }

    if (!PySequence_Check(value)) {
      PyErr_SetString(PyExc_TypeError, "Values in atoms.properties.values() must be of sequence type");
      goto fail;
    }

    if (PySequence_Size(value) != 3) {
      PyErr_SetString(PyExc_ValueError, "Values in atoms.properties.values() must be of length 3.");
      goto fail;
    }

    for (j=0; j<3; j++) {
      p = PySequence_GetItem(valuetuple, j);
      lookup[j] = PyInt_AsLong(p);
      if (lookup[j] == -1) goto fail;
      Py_DECREF(p);
    }
    strcpy(atomeye_atoms->property_name[i], PyString_AsString(key));
    atomeye_atoms->property_type[i] = lookup[0];
    atomeye_atoms->property_start[i] = lookup[1]-1;
    atomeye_atoms->property_ncols[i] = lookup[2]-lookup[1]+1;

    Py_DECREF(key);
    Py_DECREF(value);
    Py_DECREF(valuetuple);
    Py_DECREF(item);
  }

  /* atoms.lattice */
  if ((lattice = PyObject_GetAttrString(pyat, "lattice")) == NULL) {
    PyErr_SetString(PyExc_AttributeError, "Missing atoms.lattice");
    goto fail;
  }
  if (PyArray_NDIM(lattice) != 2) {
    PyErr_SetString(PyExc_ValueError, "atoms.lattice must have 2 dimensions");
    goto fail;
  }
  if (PyArray_DIM(lattice,0) != 3 || PyArray_DIM(lattice,1) != 3) {
    PyErr_SetString(PyExc_ValueError, "atoms.lattice must have shape (3,3)");
    goto fail;
  }
  if (PyArray_TYPE(lattice) != NPY_DOUBLE) {
    PyErr_SetString(PyExc_ValueError, "atoms.lattice must have type double");
    goto fail;
  }
  for (i=0; i<3; i++)
    for (j=0; j<3; j++)
      atomeye_atoms->lattice[i][j] = *(double *)PyArray_GETPTR2(lattice, i, j);

  /* atoms.data.int */
  if (atomeye_atoms->n_int > 0) {
    if ((intarray = PyObject_GetAttrString(data, "int")) == NULL) {
      PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.int");
      goto fail;
    }
    if (PyArray_NDIM(intarray) != 2) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.int must have 2 dimensions");
      goto fail;
    }
    if (PyArray_DIM(intarray, 0) != atomeye_atoms->n_int || PyArray_DIM(intarray, 1) != atomeye_atoms->n_atom) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.int must have shape (atoms.data.intsize,atoms.n)");
      goto fail;
    }
    if (PyArray_TYPE(intarray) != NPY_INT) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.int must have type int");
      goto fail;
    }
    if (!PyArray_ISFORTRAN(intarray)) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.int must be Fortran-contigous");
      goto fail;
    }
    atomeye_atoms->int_data = PyArray_DATA(intarray);
  } else
    atomeye_atoms->int_data = NULL;

  /* atoms.data.real */
  if (atomeye_atoms->n_real > 0) {
    if ((realarray = PyObject_GetAttrString(data, "real")) == NULL) {
      PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.real");
      goto fail;
    }
    if (PyArray_NDIM(realarray) != 2) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.real must have 2 dimensions");
      goto fail;
    }
    if (PyArray_DIM(realarray, 0) != atomeye_atoms->n_real || PyArray_DIM(realarray, 1) != atomeye_atoms->n_atom) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.real must have shape (atoms.data.realsize,atoms.n)");
      goto fail;
    }
    if (PyArray_TYPE(realarray) != NPY_DOUBLE) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.real must have type real");
      goto fail;
    }
    if (!PyArray_ISFORTRAN(realarray)) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.real must be Fortran-contigous");
      goto fail;
    }
    atomeye_atoms->real_data = PyArray_DATA(realarray);
  } else
    atomeye_atoms->real_data = NULL;

  /* atoms.data.str */
  if (atomeye_atoms->n_str > 0) {
    if ((strarray = PyObject_GetAttrString(data, "str")) == NULL) {
      PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.str");
      goto fail;
    }
    if (PyArray_NDIM(strarray) != 3) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.str must have 3 dimensions");
      goto fail;
    }
    if (PyArray_DIM(strarray, 0) != PROPERTY_STRING_LENGTH ||
	PyArray_DIM(strarray, 1) != atomeye_atoms->n_str || 
	PyArray_DIM(strarray, 2) != atomeye_atoms->n_atom) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.str must have shape (atoms.data.strsize,atoms.n)");
      goto fail;
    }
    if (PyArray_TYPE(strarray) != NPY_STRING) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.str must have type str");
      goto fail;
    }
    if (!PyArray_ISFORTRAN(strarray)) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.str must be Fortran-contigous");
      goto fail;
    }
    atomeye_atoms->str_data = PyArray_DATA(strarray);
  } else
    atomeye_atoms->str_data = NULL;

  /* atoms.data.logical */
  if (atomeye_atoms->n_logical > 0) {
    if ((logicalarray = PyObject_GetAttrString(data, "logical")) == NULL) {
      PyErr_SetString(PyExc_AttributeError, "Missing attribute atoms.data.logical");
      goto fail;
    }
    if (PyArray_NDIM(logicalarray) != 2) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.logical must have 2 dimensions");
      goto fail;
    }
    if (PyArray_DIM(logicalarray, 0) != atomeye_atoms->n_logical || PyArray_DIM(logicalarray, 1) != atomeye_atoms->n_atom) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.logical must have shape (atoms.data.logicalsize,atoms.n)");
      goto fail;
    }
    if (PyArray_TYPE(logicalarray) != NPY_BOOL && PyArray_TYPE(logicalarray) != NPY_INT) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.logical must have type logical or type int");
      goto fail;
    }
    if (!PyArray_ISFORTRAN(logicalarray)) {
      PyErr_SetString(PyExc_ValueError, "atoms.data.logical must be Fortran-contigous");
      goto fail;
    }
    atomeye_atoms->logical_data = PyArray_DATA(logicalarray);
  } else
    atomeye_atoms->logical_data = NULL;

  Py_DECREF(n);                       
  Py_DECREF(data);
  Py_DECREF(intsize);
  Py_DECREF(realsize);
  Py_DECREF(strsize);
  Py_DECREF(logicalsize);
  Py_DECREF(properties);
  Py_DECREF(items);
  Py_DECREF(lattice);
  Py_XDECREF(intarray);
  Py_XDECREF(realarray);
  Py_XDECREF(strarray);
  Py_XDECREF(logicalarray);
  return 1;

 fail:
  Py_XDECREF(n);                       
  Py_XDECREF(data);
  Py_XDECREF(intsize);
  Py_XDECREF(realsize);
  Py_XDECREF(strsize);
  Py_XDECREF(logicalsize);
  Py_XDECREF(properties);
  Py_XDECREF(items);
  Py_XDECREF(intarray);
  Py_XDECREF(realarray);
  Py_XDECREF(strarray);
  Py_XDECREF(logicalarray);
  Py_XDECREF(lattice);
  Py_XDECREF(key);
  Py_XDECREF(value);
  Py_XDECREF(valuetuple);
  Py_XDECREF(item);
  return 0;
}

static char atomeye_open_window_doc[]=
  "iw = _atomeye.open_window(copy=-1,atoms=None,nowindow=) -- open a new AtomEye window";

static PyObject*
atomeye_open_window(PyObject *self, PyObject *args)
{
  int icopy = -1, iw, argc;
  char outstr[255];
  char *argv[3];
  PyObject *pyat = NULL;
  static int atomeye_initialised = 0;
  int nowindow = 0;

  if (!PyArg_ParseTuple(args, "|iOi", &icopy, &pyat, &nowindow))
    return NULL;

  if (!atomeye_initialised) {
    atomeye_atoms = malloc(sizeof(Atoms));  // Who will free this structure?
    if (atomeye_atoms == NULL) return 0;
    atoms_init(atomeye_atoms);

    argv[0] = (char *)malloc(20);
    argv[1] = (char *)malloc(20);
    strcpy(argv[0], "A");
    strcpy(argv[1], "-nostdin");
    argc = 2;
    if (nowindow) {
      strcpy(argv[2], "-nowindow");
      argc = 3;
    }
  
    if (pyat != NULL && pyat != Py_None) {
      if (!update_atoms_structure(pyat)) return NULL;
      atomeyelib_init(argc, argv, (void *)atomeye_atoms);
    } else
      atomeyelib_init(argc, argv, NULL);

    atomeyelib_set_handlers(&on_click_atom, &on_close, &on_advance, &on_new);

    free(argv[0]);
    free(argv[1]);

    atomeye_initialised = 1;
  }

  iw = atomeyelib_open_window(icopy);

  if (iw == -1) {
    sprintf(outstr, "Bad copy window id %d", icopy);
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  return PyInt_FromLong((long)iw);
}

static char atomeye_set_handlers_doc[]=
  "_atomeye.set_handlers(on_click, on_close, on_advance)";

static PyObject*
atomeye_set_handlers(PyObject *self, PyObject *args)
{
  if (!PyArg_ParseTuple(args, "OOOO", &on_click_atom_pyfunc, &on_close_pyfunc, &on_advance_pyfunc, &on_new_pyfunc))
    return NULL;

  if (!PyCallable_Check(on_click_atom_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  if (!PyCallable_Check(on_close_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  if (!PyCallable_Check(on_advance_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  if (!PyCallable_Check(on_new_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  Py_INCREF(on_click_atom_pyfunc);
  Py_INCREF(on_advance_pyfunc);
  Py_INCREF(on_close_pyfunc);
  Py_INCREF(on_new_pyfunc);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_redraw_doc[]=
  "_atomeye.redraw(iw) -- redraw window";

static PyObject*
atomeye_redraw(PyObject *self, PyObject *args)
{
  int iw;
  char outstr[255];

  if (!PyArg_ParseTuple(args, "i", &iw))
    return NULL;

  if (!atomeyelib_queueevent(iw, ATOMEYELIB_REDRAW, "", NULL, outstr)) {
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_run_command_doc[]=
  "_atomeye.run_command(iw, instr) -- run an AtomEye command";

static PyObject*
atomeye_run_command(PyObject *self, PyObject *args)
{
  int iw;
  char *instr;
  char outstr[255];
  
  if (!PyArg_ParseTuple(args, "is", &iw, &instr))
    return NULL;

  if (!atomeyelib_queueevent(iw, ATOMEYELIB_RUN_COMMAND, instr, NULL, outstr)) {
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_load_atoms_doc[]=
  "_atomeye.load_atoms(iw, title, atoms) -- load atoms into AtomEye";

static PyObject*
atomeye_load_atoms(PyObject *self, PyObject *args)
{
  int iw;
  char *title;
  PyObject *pyat;
  char outstr[255];
  
  if (!PyArg_ParseTuple(args, "isO", &iw, &title, &pyat))
    return NULL;
  
  if (!update_atoms_structure(pyat)) return NULL;

  if (!atomeyelib_queueevent(iw, ATOMEYELIB_LOAD_ATOMS, title, (void *)atomeye_atoms, outstr)) {
    PyErr_SetString(PyExc_RuntimeError, outstr);    
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_set_title_doc[] =
  "_atomeye.set_title(iw, title) -- set window title";

static PyObject*
atomeye_set_title(PyObject *self, PyObject *args)
{
  int iw;
  char *title;

  if (!PyArg_ParseTuple(args, "is", &iw, &title))
    return NULL;

  atomeyelib_set_title(iw, title);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_wait_doc[] =
  "_atomeye.wait(iw) -- wait for window `iw` to complete all queued events";

static PyObject*
atomeye_wait(PyObject *self, PyObject *args)
{
  int iw;

  if (!PyArg_ParseTuple(args, "i", &iw))
    return NULL;
  
  atomeyelib_wait(iw);
  return Py_None;
}


static PyMethodDef atomeye_methods[] = {
  {"open_window", atomeye_open_window, METH_VARARGS, atomeye_open_window_doc},
  {"set_handlers", atomeye_set_handlers, METH_VARARGS, atomeye_set_handlers_doc},
  {"redraw", atomeye_redraw, METH_VARARGS, atomeye_redraw_doc},
  {"run_command", atomeye_run_command, METH_VARARGS, atomeye_run_command_doc},
  {"load_atoms", atomeye_load_atoms, METH_VARARGS, atomeye_load_atoms_doc},
  {"set_title", atomeye_set_title, METH_VARARGS, atomeye_set_title_doc},
  {"wait", atomeye_wait, METH_VARARGS, atomeye_wait_doc},
  {NULL, NULL}
};

PyMODINIT_FUNC
init_atomeye(void)
{
  Py_InitModule3("_atomeye", atomeye_methods, atomeye_doc);
  import_array();
}

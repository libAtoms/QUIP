#include <Python.h>
#include <pthread.h>
#include <signal.h>

#include <xyz_netcdf.h>
#include <atomeyelib.h>

static char atomeye_doc[] = 
"This module interfaces to AtomEye.";

static PyObject *on_click_atom_pyfunc = NULL;
static PyObject *on_advance_pyfunc = NULL;

static int atomeye_initialised = 0;

static void on_click_atom(int atom)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  arglist = Py_BuildValue("(i)", atom);               // Build argument list
  PyEval_CallObject(on_click_atom_pyfunc, arglist);   // Call Python
  Py_DECREF(arglist);                                 // Trash arglist
  PyGILState_Release(state);
}

static void on_advance(char *instr)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  arglist = Py_BuildValue("(s)", instr);               // Build argument list
  PyEval_CallObject(on_advance_pyfunc, arglist);       // Call Python
  Py_DECREF(arglist);                                 // Trash arglist
  PyGILState_Release(state);
}


static void on_close(void)
{
  fprintf(stderr, "atomeyemodule on_close handler called.\n");
  atomeye_initialised = 0;
}


static char atomeye_start_doc[] = 
  "start(on_click_handler, on_advance_handler) -- start AtomEye.";

static PyObject*
atomeye_start(PyObject *self, PyObject *args)
{
  char *argv[2];

  if (!PyArg_ParseTuple(args, "OO", &on_click_atom_pyfunc, &on_advance_pyfunc))
    return NULL;

  if (!PyCallable_Check(on_click_atom_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  if (!PyCallable_Check(on_advance_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  Py_INCREF(on_click_atom_pyfunc);
  Py_INCREF(on_advance_pyfunc);

  argv[0] = (char *)malloc(20);
  argv[1] = (char *)malloc(20);
  strcpy(argv[0], "A");
  strcpy(argv[1], "-nostdin");
  
  Py_BEGIN_ALLOW_THREADS;
  atomeyelib_main(2, argv, &on_click_atom, &on_close, &on_advance, &atomeye_initialised);
  Py_END_ALLOW_THREADS;

  free(argv[0]);
  free(argv[1]);
  free(argv[2]);

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_close_doc[] =
  "close(thread_id) -- close AtomEye thread `thread_id`.";

static PyObject*
atomeye_close(PyObject *self, PyObject *args)
{
  int iw;
  if (!PyArg_ParseTuple(args, "i", &iw))
    return NULL;
  
  atomeyelib_close(iw);

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_run_command_doc[] =
  "run_command(window_id, command) -- send `command` to AtomEye window `window_id`.";

static PyObject*
atomeye_run_command(PyObject *self, PyObject *args)
{
  char *command;
  int iw;
  char *outstr = NULL;

  if (!PyArg_ParseTuple(args, "is", &iw, &command))
    return NULL;

  if (!atomeye_initialised) {
    PyErr_SetString(PyExc_RuntimeError, "AtomEye is not initialised");
    return NULL;
  }
  
  atomeyelib_run_command(iw, command, &outstr);

  if (outstr != NULL) {
    return PyString_FromString(outstr);
  }

  Py_INCREF(Py_None);
  return Py_None;  
}

static char atomeye_help_doc[] =
  "help(window_id, command) -- send `command` to AtomEye window `window_id`.";

static PyObject*
atomeye_help(PyObject *self, PyObject *args)
{
  char *command;
  int iw;
  char *outstr = NULL;

  if (!PyArg_ParseTuple(args, "is", &iw, &command))
    return NULL;

  if (!atomeye_initialised) {
    PyErr_SetString(PyExc_RuntimeError, "AtomEye is not initialised");
    return NULL;
  }
  
  atomeyelib_help(iw, command, &outstr);

  if (outstr != NULL) {
    return PyString_FromString(outstr);
  }

  Py_INCREF(Py_None);
  return Py_None;  
}


static char atomeye_redraw_doc[] =
  "redraw(window_id) -- redraw AtomEye window `window_id`.";

static PyObject*
atomeye_redraw(PyObject *self, PyObject *args)
{
  int iw;

  if (!PyArg_ParseTuple(args, "i", &iw))
    return NULL;

  if (!atomeye_initialised) {
    PyErr_SetString(PyExc_RuntimeError, "AtomEye is not initialised");
    return NULL;
  }
  
  atomeyelib_redraw(iw);

  Py_INCREF(Py_None);
  return Py_None;  
}

static char atomeye_load_libatoms_doc[] =
  "load_libatoms(window_id, atoms_ptr, title) -- load from libAtoms Atoms object in memory.";

static PyObject*
atomeye_load_libatoms(PyObject *self, PyObject *args)
{
  int iw;
  Atoms *at_ptr;
  char *title, *outstr = NULL;

  if (!PyArg_ParseTuple(args, "iLs", &iw, &at_ptr, &title))
    return NULL;

  if (!atomeye_initialised) {
    PyErr_SetString(PyExc_RuntimeError, "AtomEye is not initialised");
    return NULL;
  }

  atomeyelib_load_libatoms(iw, at_ptr, title, &outstr);
  if (outstr != NULL) {
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_isAlive_doc[] =
  "isAlive() -- return whether AtomEye has been initialised.";

static PyObject*
atomeye_isAlive(PyObject *self, PyObject *args)
{
  return PyBool_FromLong((long)atomeye_initialised);
}

static PyMethodDef atomeye_methods[] = {
  {"start", atomeye_start, METH_VARARGS, atomeye_start_doc},
  {"close", atomeye_close, METH_VARARGS, atomeye_close_doc},
  {"run_command", atomeye_run_command, METH_VARARGS, atomeye_run_command_doc},
  {"help", atomeye_help, METH_VARARGS, atomeye_help_doc},
  {"redraw", atomeye_redraw, METH_VARARGS, atomeye_redraw_doc},
  {"load_libatoms", atomeye_load_libatoms, METH_VARARGS, atomeye_load_libatoms_doc},
  {"isAlive", atomeye_isAlive, METH_VARARGS, atomeye_isAlive_doc},
  {NULL, NULL}
};

PyMODINIT_FUNC
init_atomeye(void)
{
  Py_InitModule3("_atomeye", atomeye_methods, atomeye_doc);
}

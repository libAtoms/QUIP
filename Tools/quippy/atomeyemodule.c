#include <Python.h>
#include <pthread.h>
#include <signal.h>

#include <xyz_netcdf.h>
#include <atomeyelib.h>

static char atomeye_doc[] = 
"This module interfaces to AtomEye.";

static char atomeye_start_doc[] = 
  "start(filename) -- start AtomEye and load from `filename`.";

static PyObject *on_click_atom_pyfunc = NULL;

static void on_click_atom(int atom)
{
  PyObject *arglist;
  PyGILState_STATE state;

  state = PyGILState_Ensure();
  arglist = Py_BuildValue("(i)", atom+1);               // Build argument list
  PyEval_CallObject(on_click_atom_pyfunc, arglist);     // Call Python
  Py_DECREF(arglist);                                   // Trash arglist
  PyGILState_Release(state);
}

static PyObject*
atomeye_start(PyObject *self, PyObject *args)
{
  char *argv[3];
  const char *filename;

  if (!PyArg_ParseTuple(args, "sO", &filename, &on_click_atom_pyfunc))
    return NULL;

  if (!PyCallable_Check(on_click_atom_pyfunc)) {
    PyErr_SetString(PyExc_TypeError, "Need a callable object!");
    return NULL;
  }

  Py_INCREF(on_click_atom_pyfunc);

  argv[0] = (char *)malloc(20);
  argv[1] = (char *)malloc(20);
  argv[2] = (char *)malloc(strlen(filename)+1);
  strcpy(argv[0], "A");
  strcpy(argv[1], "-nostdin");
  strcpy(argv[2], filename);
  
  Py_BEGIN_ALLOW_THREADS;
  atomeyelib_main(3, argv, &on_click_atom);
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
  
  atomeyelib_run_command(iw, command, &outstr);

  if (outstr != NULL) {
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
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

  if (!PyArg_ParseTuple(args, "ils", &iw, &at_ptr, &title))
    return NULL;

  atomeyelib_load_libatoms(iw, at_ptr, title, &outstr);
  if (outstr != NULL) {
    PyErr_SetString(PyExc_RuntimeError, outstr);
    return NULL;
  }

  Py_INCREF(Py_None);
  return Py_None;
}

static char atomeye_set_output_doc[] =
  "set_output(on_off) -- enable or disable output to stdout and stderr.";

static PyObject*
atomeye_set_output(PyObject *self, PyObject *args)
{
  int on_off;
  
  if (!PyArg_ParseTuple(args, "i", &on_off))
    return NULL;  

  atomeyelib_set_output(on_off);

  Py_INCREF(Py_None);
  return Py_None;  
}

static PyMethodDef atomeye_methods[] = {
  {"start", atomeye_start, METH_VARARGS, atomeye_start_doc},
  {"close", atomeye_close, METH_VARARGS, atomeye_close_doc},
  {"run_command", atomeye_run_command, METH_VARARGS, atomeye_run_command_doc},
  {"redraw", atomeye_redraw, METH_VARARGS, atomeye_redraw_doc},
  {"load_libatoms", atomeye_load_libatoms, METH_VARARGS, atomeye_load_libatoms_doc},
  {"set_output", atomeye_set_output, METH_VARARGS, atomeye_set_output_doc},
  {NULL, NULL}
};

PyMODINIT_FUNC
init_atomeye(void)
{
  Py_InitModule3("_atomeye", atomeye_methods, atomeye_doc);
}

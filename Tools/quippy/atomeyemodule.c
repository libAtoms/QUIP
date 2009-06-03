#include <Python.h>
#include <pthread.h>
#include <signal.h>

#include <xyz_netcdf.h>
#include <atomeyelib.h>

static char atomeye_doc[] = 
"This module interfaces to AtomEye.";

static PyObject *on_click_atom_pyfunc = NULL;
static PyObject *on_advance_pyfunc = NULL;
static PyObject *on_close_pyfunc = NULL;
static PyObject *on_new_pyfunc = NULL;

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

static char atomeye_open_window_doc[]=
  "iw = _atomeye.open_window(copy=-1,atoms=NULL) -- open a new AtomEye window";

static PyObject*
atomeye_open_window(PyObject *self, PyObject *args)
{
  int icopy = -1, iw;
  char outstr[255];
  char *argv[2];
  void *data = NULL;
  static int atomeye_initialised = 0;

  if (!PyArg_ParseTuple(args, "|iL", &icopy, &data))
    return NULL;

  if (!atomeye_initialised) {
    argv[0] = (char *)malloc(20);
    argv[1] = (char *)malloc(20);
    strcpy(argv[0], "A");
    strcpy(argv[1], "-nostdin");
  
    atomeyelib_init(2, argv, data);
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
  "_atomeye.run_command(iw, instr, data) -- run an AtomEye command";

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
  "_atomeye.load_atoms(iw, title, data) -- load atoms into AtomEye";

static PyObject*
atomeye_load_atoms(PyObject *self, PyObject *args)
{
  int iw;
  char *title;
  void *data;
  char outstr[255];
  
  if (!PyArg_ParseTuple(args, "isL", &iw, &title, &data))
    return NULL;

  if (!atomeyelib_queueevent(iw, ATOMEYELIB_LOAD_ATOMS, title, data, outstr)) {
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


static PyMethodDef atomeye_methods[] = {
  {"open_window", atomeye_open_window, METH_VARARGS, atomeye_open_window_doc},
  {"set_handlers", atomeye_set_handlers, METH_VARARGS, atomeye_set_handlers_doc},
  {"redraw", atomeye_redraw, METH_VARARGS, atomeye_redraw_doc},
  {"run_command", atomeye_run_command, METH_VARARGS, atomeye_run_command_doc},
  {"load_atoms", atomeye_load_atoms, METH_VARARGS, atomeye_load_atoms_doc},
  {"set_title", atomeye_set_title, METH_VARARGS, atomeye_set_title_doc},
  {NULL, NULL}
};

PyMODINIT_FUNC
init_atomeye(void)
{
  Py_InitModule3("_atomeye", atomeye_methods, atomeye_doc);
}

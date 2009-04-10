#include <Python.h>
#include <numpy/arrayobject.h>

static PyObject * 
arraydata(PyObject *self, PyObject *args)
{
  int nd, i;
  npy_intp *dimensions;
  char *data;
  PyArrayObject *array;
  PyObject *shape, *type, *item;
  PyArray_Descr *descr;

  if (!PyArg_ParseTuple(args, "OOl", &shape,&type,&data))
    return NULL;

  if (!PySequence_Check(shape) && !PyInt_Check(shape)) {
    PyErr_SetString(PyExc_ValueError, "shape must be a sequence or an integer");
    return NULL;
  }

  if (PyInt_Check(shape)) {
    // shape is a single integer
    nd = 1;
    dimensions = (npy_intp*)malloc(1*sizeof(npy_intp));
    dimensions[0] = (npy_intp)PyInt_AsLong(shape);
  }
  else {
    // shape is a sequence
    nd = (int)PySequence_Length(shape);

    if (nd == 0) {
      PyErr_SetString(PyExc_ValueError, "len(shape) must be > 0");
      return NULL;
    }
    dimensions = (npy_intp*)malloc(nd*sizeof(int));
    for (i=0; i<nd; i++) {
      item = PySequence_GetItem(shape,i);
      dimensions[i] = (npy_intp)PyInt_AsLong(item);
      Py_DECREF(item);
    }
  }

  if (!PyArray_DescrConverter(type, &descr)) {
    PyErr_SetString(PyExc_ValueError, "dtype must be a numpy data-type");
    return NULL;
  }

//  array = (PyArrayObject*) PyArray_FromDimsAndDataAndDescr(nd, dimensions, descr, data);
  array = (PyArrayObject*) PyArray_NewFromDescr(&PyArray_Type, descr, nd, dimensions, NULL, 
                                                data, NPY_FORTRAN | NPY_WRITEABLE, NULL);
  free(dimensions);
  return (PyObject *)array;
}

static PyMethodDef arraydata_methods[] = {
  {"arraydata", arraydata, METH_VARARGS, 
   "Make an array from shape, dtype and pointer to data.\n\narraydata((d1,...,dn), dtype, dataptr) -> array"},
  {NULL, NULL}
};

static char arraydata_doc[] = 
  "Extension module to create numpy arrays from pointers to memory location";

PyMODINIT_FUNC
initarraydata(void)
{
  Py_InitModule3("arraydata", arraydata_methods, arraydata_doc);
  import_array();
}



#include <Python.h>
#include <numpy/arrayobject.h>

#include <xyz_netcdf.h>

static Atoms *atoms;

static PyObject *
py_xyz_netcdf_open(PyObject *self, PyObject *args)
{
  char *filename;
  int action;
  int append;

  cio_init(&atoms, filename, &action, &append, NULL, NULL, NULL, NULL, NULL, NULL,
	   NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
	   NULL, NULL, NULL, NULL, NULL, NULL);
  
  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
py_xyz_netcdf_close(PyObject *self, PyObject *args)
{
  cio_free(atoms);

  Py_INCREF(Py_None);
  return Py_None;
}

static PyObject*
py_xyz_netcdf_query(PyObject *self, PyObject *args)
{
  int frame, n_frame;
  frame = 0; // default to frame zero.

  if (!PyArg_ParseTuple(args, "|i", &frame))
    return NULL;

  n_frame = cio_query(atoms, &frame);

  return Py_BuildValue("(i)", n_frame);
}

static PyObject*
py_xyz_netcdf_write(PyObject *self, PyObject *args)
{
}

static PyObject * 
py_xyz_netcdf_read(PyObject *self, PyObject *args)
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

static PyMethodDef xyz_netcdf_methods[] = {
  {"xyz_netcdf", xyz_netcdf, METH_VARARGS, 
   "Make an array from shape, dtype and pointer to data.\n\nxyz_netcdf((d1,...,dn), dtype, dataptr) -> array"},
  {NULL, NULL}
};

static char xyz_netcdf_doc[] = 
  "Extension module to read/write to/from XYZ and NetCDF files";

PyMODINIT_FUNC
initxyz_netcdf(void)
{
  Py_InitModule3("xyz_netcdf", xyz_netcdf_methods, xyz_netcdf_doc);
  import_array();
}



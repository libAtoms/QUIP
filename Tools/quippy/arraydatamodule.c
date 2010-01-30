#include <Python.h>
#include <fortranobject.h>

static PyObject*
arraydata(PyObject *self, PyObject *args)
{
  typedef void (*arrayfunc_t)(int*,int*,int*,int*,void*);

  int nd, i, typenum;
  int dim_temp[10];
  npy_intp *dimensions;
  char *data = NULL;
  PyArrayObject *array = NULL;
  PyArray_Descr *descr = NULL;

  int *this = NULL;
  npy_intp this_Dims[1] = {-1};
  const int this_Rank = 1;
  PyArrayObject *capi_this_tmp = NULL;
  int capi_this_intent = 0;
  PyObject *this_capi = NULL;
  PyFortranObject *arrayfunc_capi = NULL;

  if (!PyArg_ParseTuple(args, "OO", &this_capi,&arrayfunc_capi))
    return NULL;

  /* Processing variable this */
  this_Dims[0]=12;
  capi_this_intent |= F2PY_INTENT_IN;
  capi_this_tmp = array_from_pyobj(PyArray_INT,this_Dims,this_Rank,capi_this_intent,this_capi);
  if (capi_this_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(PyExc_TypeError,"failed in converting 1st argument `this' of arraydata to C/Fortran array" );
    goto fail;
  } else {
    this = (int *)(capi_this_tmp->data);
  }

  /* Processing variable arrayfunc */
  if (!PyFortran_Check1(arrayfunc_capi)) {
    PyErr_SetString(PyExc_TypeError, "2nd argument `arrayfunc' is not a fortran object");
    goto fail;
  }
  
  if (arrayfunc_capi->defs[0].rank==-1) {/* is Arrayfunc_Capirtran routine */
    if ((arrayfunc_capi->defs[0].func==NULL)) {
      PyErr_Format(PyExc_RuntimeError, "no function to call");
      goto fail;
    }
    else if (arrayfunc_capi->defs[0].data==NULL) {
      PyErr_Format(PyExc_TypeError, "fortran object is not callable");
      goto fail;
    }
  } else {
    PyErr_Format(PyExc_TypeError, "fortran object is not callable");
    goto fail;
  }

  /* Call arrayfunc_capi routine */
  ((arrayfunc_t)(arrayfunc_capi->defs[0].data))(this, &nd, &typenum, dim_temp, &data);

  if (data == NULL) {
    PyErr_SetString(PyExc_ValueError, "array is NULL");
    goto fail;
  }

  dimensions = (npy_intp*)malloc(nd*sizeof(npy_intp));
  for (i=0; i<nd; i++) {
    dimensions[i] = (npy_intp)(dim_temp[i]);
  }

  /* Construct array */
  descr = PyArray_DescrNewFromType(typenum);
  array = (PyArrayObject*) PyArray_NewFromDescr(&PyArray_Type, descr, nd, dimensions, NULL, 
                                                data, NPY_FORTRAN | NPY_WRITEABLE, NULL);
  free(dimensions);
  if((PyObject *)capi_this_tmp!=this_capi) {
    Py_XDECREF(capi_this_tmp);
  }
  return (PyObject *)array;

 fail:
  Py_XDECREF(descr);
  if(capi_this_tmp != NULL && ((PyObject *)capi_this_tmp!=this_capi)) {
    Py_XDECREF(capi_this_tmp);
  }
  return NULL;
}

static PyMethodDef arraydata_methods[] = {
  {"arraydata", arraydata, METH_VARARGS, 
   "Make an array from integer(12) array containing reference to derived type object,\n and fortran array function.\n\narraydata(fpointer,array_fobj) -> array"},
  {NULL, NULL}
};

static char arraydata_doc[] = 
  "Extension module to create numpy arrays which access existing data at a given memory location";

PyMODINIT_FUNC
initarraydata(void)
{
  Py_InitModule3("arraydata", arraydata_methods, arraydata_doc);
  PyFortran_Type.ob_type = &PyType_Type;
  import_array();
}



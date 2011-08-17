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
#include <fortranobject.h>
#include <libatoms.h>

static PyObject*
get_array(PyObject *self, PyObject *args)
{
  typedef void (*arrayfunc_t)(int*,int*,int*,int*,void*);
  typedef void (*arrayfunc_key_t)(int*,char*,int*,int*,int*,void*,int);

  int nd, i, type, typenum;
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
  char *key = NULL;

  if (!PyArg_ParseTuple(args, "OO|s", &this_capi,&arrayfunc_capi,&key))
    return NULL;

  /* Processing variable this */
  this_Dims[0]=SIZEOF_FORTRAN_T;
  capi_this_intent |= F2PY_INTENT_IN;
  capi_this_tmp = array_from_pyobj(PyArray_INT,this_Dims,this_Rank,capi_this_intent,this_capi);
  if (capi_this_tmp == NULL) {
    if (!PyErr_Occurred())
      PyErr_SetString(PyExc_TypeError,"failed in converting 1st argument `this' of get_array to C/Fortran array" );
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
  if (key == NULL) 
    ((arrayfunc_t)(arrayfunc_capi->defs[0].data))(this, &nd, &type, dim_temp, &data);
  else
    ((arrayfunc_key_t)(arrayfunc_capi->defs[0].data))(this, key, &nd, &type, dim_temp, &data, strlen(key));

  if (data == NULL) {
    PyErr_SetString(PyExc_ValueError, "array is NULL");
    goto fail;
  }

  // Convert from libAtoms type code to numpy typenum
  switch(type) {
  case(T_INTEGER_A):
  case(T_LOGICAL_A):
  case(T_INTEGER_A2):
    typenum = NPY_INT32;
    break;
  case(T_REAL_A):
  case(T_REAL_A2):
    typenum = NPY_DOUBLE;
    break;
  case(T_COMPLEX_A):
    typenum = NPY_COMPLEX128;
    break;
  case(T_CHAR_A):
    typenum = NPY_CHAR;
    break;
  default:
    PyErr_Format(PyExc_TypeError, "Unknown data type %d", type);
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
  {"get_array", get_array, METH_VARARGS, 
   "Make an array from integer(SIZEOF_FORTRAN_T) array containing reference to derived type object,\n and fortran array function.\n\get_array(fpointer,array_fobj[,key]) -> array"},
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



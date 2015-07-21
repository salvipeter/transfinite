#include "Python.h"

#include "surface.hh"

extern "C" {

static PyObject*
transfiniteSurface(PyObject* self, PyObject* args) {
  // As a test of using the geometry library, take a double value x,
  // and compute the length of the vector (x,x,x).
  double x;
  if (!PyArg_ParseTuple(args, "d", &x))
    return NULL;
  Vector3D v(x, x, x);
  return Py_BuildValue("d", v.norm());
}

static PyMethodDef PyTranSurfMethods[] = {
  {"transfiniteSurface", transfiniteSurface, METH_VARARGS, "Create a transfinite surface."},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
initpytransurf(void) {
  (void) Py_InitModule("pytransurf", PyTranSurfMethods);
}

}

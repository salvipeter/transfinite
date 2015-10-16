#include "Python.h"

#include <algorithm>
#include <cstring>
#include <exception>

#include "surface-side-based.hh"
#include "surface-corner-based.hh"
#include "surface-generalized-coons.hh"
#include "surface-composite-ribbon.hh"

#ifdef WIN32
#define DLL_EXPORT __declspec(dllexport)
#else // !WIN32
#define DLL_EXPORT
#endif // WIN32

using Transfinite::DoubleVector;
using Transfinite::Point3D;
using Transfinite::PointVector;
using Transfinite::BSCurve;
using Transfinite::CurveVector;
using Transfinite::BSSurface;

class ParseException : public std::runtime_error {
public:
  ParseException(const std::string &message) : std::runtime_error(message) { }
};

// Conversions

DoubleVector
parseDoubleList(PyObject *py_list) {
  DoubleVector result;
  int n = PyList_Size(py_list);
  if (n < 0)
    throw ParseException("invalid double list");
  result.reserve(n);
  for (int i = 0; i < n; ++i) {
    double x = PyFloat_AsDouble(PyList_GetItem(py_list, i));
    if (PyErr_Occurred())
      throw ParseException("invalid double list");
    result.push_back(x);
  }
  return result;
}

// A point is a list: [x, y, z]
Point3D
parsePoint(PyObject *py_point) {
  DoubleVector p = parseDoubleList(py_point);
  if (p.size() != 3)
    throw ParseException("not a 3D point");
  return Point3D(p[0], p[1], p[2]);
}

PointVector
parsePointList(PyObject *py_list) {
  PointVector result;
  int n = PyList_Size(py_list);
  if (n < 0)
    throw ParseException("invalid point list");
  result.reserve(n);
  for (int i = 0; i < n; ++i)
    result.push_back(parsePoint(PyList_GetItem(py_list, i)));
  return result;
}

// A curve is a tuple: (degree, [knots], [points])
std::shared_ptr<BSCurve>
parseCurve(PyObject *py_curve) {
  PyObject *py_knots, *py_points;
  double degree;
  if (!PyArg_ParseTuple(py_curve, "dOO", &degree, &py_knots, &py_points))
    throw ParseException("invalid curve");
  return std::make_shared<BSCurve>(degree, parseDoubleList(py_knots), parsePointList(py_points));
}

CurveVector
parseCurveList(PyObject *py_list) {
  CurveVector result;
  int n = PyList_Size(py_list);
  if (n < 0)
    throw ParseException("invalid curve list");
  result.reserve(n);
  for (int i = 0; i < n; ++i)
    result.push_back(parseCurve(PyList_GetItem(py_list, i)));
  return result;
}

PyObject *
buildPyList(const std::vector<PyObject *> &list) {
  size_t n = list.size();
  PyObject *py_list = PyList_New(n);
  if (!py_list)
    return NULL;
  for (size_t i = 0; i < n; ++i)
    PyList_SET_ITEM(py_list, i, list[i]);
  return py_list;
}

PyObject *
buildPyDoubleList(const std::vector<double> &list) {
  std::vector<PyObject *> py_list;
  std::transform(list.begin(), list.end(), std::back_inserter(py_list),
                 [](double x) { return Py_BuildValue("d", x); });
  return buildPyList(py_list);
}

PyObject *
buildPyPoint(const Point3D &point) {
  DoubleVector py_point = { point[0], point[1], point[2] };
  return buildPyDoubleList(py_point);
}

PyObject *
buildPyPoint2D(const Point3D &point) {
  DoubleVector py_point = { point[0], point[1] };
  return buildPyDoubleList(py_point);
}

PyObject *
buildPyCurve2D(const std::shared_ptr<BSCurve> &curve) {
  std::vector<PyObject *> py_points;
  size_t n = curve->nrControlPoints();
  py_points.reserve(n);

  for (size_t i = 0; i < n; ++i)
    py_points.push_back(buildPyPoint2D(curve->controlPoint(i)));
  
  return Py_BuildValue("(IOO)", curve->degree(),
                       buildPyDoubleList(curve->knotVector()),
                       buildPyList(py_points));
}

// A surface is a tuple: (deg_u, deg_v, [knots_u], [knots_v], [points])
PyObject *
buildPySurface(const BSSurface &surf) {
  std::vector<PyObject *> py_points;
  size_t n = surf.cnet_.size();
  if (n == 0)
    return NULL;
  size_t m = surf.cnet_[0].size();
  py_points.reserve(n * m);

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j)
      py_points.push_back(buildPyPoint(surf.cnet_[i][j]));
  
  return Py_BuildValue("(IIOOO)",
                       surf.deg_u_, surf.deg_v_,
                       buildPyDoubleList(surf.knots_u_), buildPyDoubleList(surf.knots_v_),
                       buildPyList(py_points));
}

// A trimmed surface is a tuple: (surface, [param_curves])
PyObject *
buildPyTrimmedSurface(const BSSurface &surf) {
  PyObject *py_surf = buildPySurface(surf);

  std::vector<PyObject *> py_param_curves; py_param_curves.reserve(surf.param_curves_.size());
  std::transform(surf.param_curves_.begin(), surf.param_curves_.end(),
                 std::back_inserter(py_param_curves),
                 [](const std::shared_ptr<BSCurve> &c) { return buildPyCurve2D(c); });

  return Py_BuildValue("(OO)", py_surf, buildPyList(py_param_curves));
}

// Python interface

extern "C" {

static PyObject *
transfiniteSurface(PyObject* self, PyObject* args) {
  PyObject *py_cv;
  char *surface_type, *split_or_trim;
  double fit_tol;
  if (!PyArg_ParseTuple(args, "Ossd", &py_cv, &surface_type, &split_or_trim, &fit_tol)) {
    PyErr_SetString(PyExc_RuntimeError, "wrong argument list");
    return NULL;
  }

  CurveVector cv;
  try {
    cv = parseCurveList(py_cv);
  } catch(ParseException e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    return NULL;
  }

  std::shared_ptr<Transfinite::Surface> surf;
  if (std::strcmp(surface_type, "SB") == 0)
    surf = std::make_shared<Transfinite::SurfaceSideBased>();
  else if (std::strcmp(surface_type, "CB") == 0)
    surf = std::make_shared<Transfinite::SurfaceCornerBased>();
  else if (std::strcmp(surface_type, "GC") == 0)
    surf = std::make_shared<Transfinite::SurfaceGeneralizedCoons>();
  else if (std::strcmp(surface_type, "CR") == 0)
    surf = std::make_shared<Transfinite::SurfaceCompositeRibbon>();
  else {
    PyErr_SetString(PyExc_RuntimeError, "invalid surface type");
    return NULL;
  }

  surf->setCurves(cv);
  surf->setupLoop();
  surf->update();

  if (std::strcmp(split_or_trim, "split") == 0) {
    std::vector<BSSurface> result = surf->fitCentralSplit(fit_tol);
    std::vector<PyObject *> py_result; py_result.reserve(result.size());
    std::transform(result.begin(), result.end(), std::back_inserter(py_result),
                   [](const BSSurface &s) { return buildPySurface(s); });
    return buildPyList(py_result);
  } else if (std::strcmp(split_or_trim, "trim") == 0) {
    BSSurface result = surf->fitTrimmed(fit_tol);
    return buildPyTrimmedSurface(result);
  } else {
    PyErr_SetString(PyExc_RuntimeError, "fit mode should be split or trim");
    return NULL;
  }

  return NULL;
}

static PyMethodDef pyTranSurfMethods[] = {
  {"transfiniteSurface", transfiniteSurface, METH_VARARGS, "Create a transfinite surface."},
  {NULL, NULL, 0, NULL}
};

#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT, "pytransurf", NULL, -1, pyTranSurfMethods, NULL, NULL, NULL, NULL
};

DLL_EXPORT PyObject *
PyInit_pytransurf() {
  return PyModule_Create(&moduledef);
}

#else  // PY_MAJOR_VERSION < 3

DLL_EXPORT PyMODINIT_FUNC
initpytransurf() {
  (void) Py_InitModule("pytransurf", pyTranSurfMethods);
}

#endif  // PY_MAJOR_VERSION

}

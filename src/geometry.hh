#pragma once

#include <vector>

#include "Point.hh"
#include "Vector.hh"
#include "BSplineCurve.hh"
#include "CurveWithRMF.hh"
#include "Triangulation.hh"

namespace Transfinite
{

  typedef Point<2, double> Point2D;
  typedef Vector<2, double> Vector2D;

  typedef Point<3, double> Point3D;
  typedef Vector<3, double> Vector3D;

  typedef std::vector<Point2D> Point2DVector;
  typedef std::vector<Point3D> Point3DVector;

  typedef BSplineCurve<Point3D> BSCurve;
  typedef std::vector<BSCurve *> CurveVector;

  typedef CurveWithRMF<Point3D> Fence;
  typedef std::vector<Fence *> FenceVector;

  typedef Triangulation<Point3D> Mesh;

  double const epsilon = 1.0e-8;

} // Transfinite

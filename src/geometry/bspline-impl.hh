#pragma once

#include "geometry.hh"

#include "BSplineCurve.hh"

class BSCurve::BSCurveImpl {
public:
  BSCurveImpl() {}
  BSCurveImpl(size_t degree, DoubleVector knots, PointVector cpts) {
    BSplineCurve<Point<3, double>>::KnotVectorType bsp_knots;
    BSplineCurve<Point<3, double>>::ControlVectorType bsp_cpts;
    bsp_knots.insert(knots.begin(), knots.end());
    for (const auto &cp : cpts) {
      bsp_cpts.push_back(Point<3, double>(cp[0], cp[1], cp[2]));
    }
    c_ = BSplineCurve<Point<3, double>>(degree, bsp_knots, bsp_cpts);
  }
  Point3D eval(double u) const {
    Vector<3, double> p(c_.Eval(u));
    return Point3D(p[0], p[1], p[2]);
  }
  Point3D eval(double u, size_t nr_der, VectorVector &der) const {
    const Vector<3, double> *bsp_der;
    Vector<3, double> p(c_.Eval(u, nr_der, bsp_der));
    der.clear();
    der.reserve(nr_der + 1);
    der.push_back(Vector3D(p[0], p[1], p[2]));
    for (size_t i = 0; i < nr_der; ++i) {
      const Vector<3, double> &v = bsp_der[i];
      der.push_back(Vector3D(v[0], v[1], v[2]));
    }
    return Point3D(p[0], p[1], p[2]);
  }
  DoubleVector knotVector() const {
    return DoubleVector(c_.theKnotVector().begin(), c_.theKnotVector().end());
  }
  void insertKnot(double k) { c_.InsertKnot(k); }
  void reverse() {
    c_.Reverse();
  }
  void normalize() {
    c_.Reparametrize(0.0, 1.0);
  }
  size_t nrControlPoints() const { return c_.NrControlValues(); }
  Point3D controlPoint(size_t i) const {
    const Point<3, double> &p = c_.theControlVector()[i];
    return Point3D(p[0], p[1], p[2]);
  }
  double arcLength(double from, double to) const {
    return c_.EstimateArcLength(from, to);
  }
  void trim(double from, double to) { c_.Trim(from, to); }
private:
  BSplineCurve<Point<3, double>> c_;
};

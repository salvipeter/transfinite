#pragma once

#include "geometry.hh"

namespace Transfinite {

using namespace Geometry;

class Curve {
public:
  virtual ~Curve() { }
  virtual Point3D eval(double u) const = 0;
  virtual Point3D eval(double u, size_t nr_der, VectorVector &der) const = 0;
  virtual void reverse() = 0;
  virtual double arcLength(double from, double to) const = 0;
};

// wrapper
class BSplineCurve : public Curve {
public:
  BSplineCurve(const BSCurve &curve) : bsc(curve) { bsc.normalize(); }
  ~BSplineCurve() { }
  Point3D eval(double u) const override { return bsc.eval(u); }
  Point3D eval(double u, size_t nr_der, VectorVector &der) const override {
    return bsc.eval(u, nr_der, der);
  }
  void reverse() override {
      bsc.reverse();
      bsc.normalize();
  }
  double arcLength(double from, double to) const override {
    return bsc.arcLength(from, to);
  }
private:
  BSCurve bsc;
};

using CurveVector = std::vector<std::shared_ptr<Curve>>;

} // namespace Transfinite

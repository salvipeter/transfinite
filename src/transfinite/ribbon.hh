#pragma once

#include "geometry.hh"
#include "rmf.hh"

namespace Transfinite {

class Ribbon {
public:
  Ribbon();
  virtual ~Ribbon();
  std::shared_ptr<const BSCurve> curve() const;
  std::shared_ptr<BSCurve> curve();
  void setCurve(const std::shared_ptr<BSCurve> &curve);
  void setNeighbors(const std::shared_ptr<Ribbon> &prev, const std::shared_ptr<Ribbon> &next);
  void setMultiplier(double m);
  void setHandler(const Vector3D &h);
  void reset();
  virtual void update();
  virtual Vector3D crossDerivative(double s) const = 0;
  virtual Point3D eval(const Point2D &sd) const;
  virtual Vector3D normal(double s) const;

protected:
  std::shared_ptr<BSCurve> curve_;
  std::weak_ptr<Ribbon> prev_, next_;
  RMF rmf_;
  Vector3D handler_;
  double multiplier_;
  bool handler_initialized_;
};

} // namespace Transfinite

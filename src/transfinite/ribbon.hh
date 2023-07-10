#pragma once

#include "geometry.hh"
#include "rmf.hh"

#include <functional>
#include <optional>

namespace Transfinite {

using NormalFence = std::function<Vector3D(double)>;

class Ribbon {
public:
  Ribbon();
  virtual ~Ribbon();
  std::shared_ptr<const Curve> curve() const;
  std::shared_ptr<Curve> curve();
  void setCurve(const std::shared_ptr<Curve> &curve);
  void setNeighbors(const std::shared_ptr<Ribbon> &prev, const std::shared_ptr<Ribbon> &next);
  double multiplier() const;
  void setMultiplier(double m);
  std::optional<Vector3D> handler() const;
  void setHandler(const Vector3D &h);
  void overrideNormalFence(const std::shared_ptr<NormalFence> &fence);
  void reset();
  virtual void update();
  virtual Vector3D crossDerivative(double s) const = 0;
  virtual Point3D eval(const Point2D &sd) const;
  Vector3D normal(double s) const;

protected:
  std::shared_ptr<Curve> curve_;
  std::weak_ptr<Ribbon> prev_, next_;
  RMF rmf_;
  std::shared_ptr<NormalFence> normal_fence_;
  Vector3D handler_;
  double multiplier_;
  bool handler_initialized_;
};

} // namespace Transfinite

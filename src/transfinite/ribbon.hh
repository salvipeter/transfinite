#pragma once

#include "geometry.hh"
#include "rmf.hh"

class Ribbon {
public:
  virtual ~Ribbon();
  std::shared_ptr<const BSCurve> curve() const;
  std::shared_ptr<BSCurve> curve();
  void setCurve(const std::shared_ptr<BSCurve> &curve);
  void setNeighbors(const std::shared_ptr<Ribbon> &prev, const std::shared_ptr<Ribbon> &next);
  virtual void update();
  virtual Point3D eval(const Point2D &sd) const = 0;
  Vector3D normal(double s) const;

protected:
  std::shared_ptr<BSCurve> curve_;
  std::shared_ptr<Ribbon> prev_, next_;
  RMF rmf_;
};

#pragma once

#include "geometry.hh"
#include "rmf.hh"

class Ribbon {
public:
  virtual ~Ribbon();
  std::shared_ptr<BSCurve> curve() const;
  void setCurve(const std::shared_ptr<BSCurve> &curve);
  void setNeighbors(const std::shared_ptr<Ribbon> &prev, const std::shared_ptr<Ribbon> &next);
  virtual void update();
  virtual Point3D eval(const Point2D &sd) const = 0;

protected:
  std::shared_ptr<BSCurve> curve_;
  std::shared_ptr<Ribbon> prev_, next_;
  RMF rmf_;
};

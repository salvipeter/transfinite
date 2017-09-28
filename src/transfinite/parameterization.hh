#pragma once

#include <map>

#include "geometry.hh"

namespace Transfinite {

using namespace Geometry;

class Domain;

class Parameterization {
public:
  virtual ~Parameterization();
  void setDomain(const std::shared_ptr<Domain> &new_domain);
  virtual void update();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const = 0;
  virtual Point2DVector mapToRibbons(const Point2D &uv) const;
  virtual Point2D inverse(size_t i, const Point2D &pd) const;

protected:
  struct PointComparator {
    bool operator()(const Point2D &p, const Point2D &q) const {
      return p[0] < q[0] || (p[0] == q[0] && p[1] < q[1]);
    }
  };

  size_t next(size_t i, size_t j = 1) const { return (i + j) % n_; }
  size_t prev(size_t i, size_t j = 1) const { return (i + n_ - j) % n_; }

  size_t n_;
  std::shared_ptr<Domain> domain_;

private:
  using ParameterCache = std::map<Point2D, Point2DVector, PointComparator>;
  mutable ParameterCache cache_;
};

} // namespace Transfinite

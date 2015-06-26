#pragma once

#include <map>
#include <memory>

#include "geometry.hh"

class Domain;

class Parameterization {
public:
  virtual ~Parameterization();
  void setDomain(const std::shared_ptr<Domain> &new_domain);
  virtual void invalidate();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const = 0;
  Point2DVector mapToRibbons(const Point2D &uv) const;

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

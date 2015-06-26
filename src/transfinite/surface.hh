#pragma once

#include "geometry.hh"

class Domain;
class Parameterization;
class Ribbon;

// Invalidation rules:
// - boundary curve changes => domain || ribbon
//     (if the domain remains the same, the invalidate routine should invalidate the ribbon)
// - domain changes => parameterization
// - ribbon handler changes => ribbon
// - resolution changes => parameterization

// Notes:
// - set... functions clear the appropriate cache data

class Surface {
public:
  enum class Invalidation {
    InvalidateNothing = 0,
    InvalidateDomain = 1,              // => InvalidateParameterization
    InvalidateParameterization = 2,    // => InvalidateBlendFunctions, InvalidateAllRibbons
    InvalidateRibbon = 3,
    InvalidateBlendFunctions = 4,
    InvalidateAllRibbons = 5
  };

  Surface(size_t sides);
  virtual ~Surface();
  void setSide(size_t i, const BSCurve *curve);
  void setSides(const CurveVector &curves);
  void setFence(size_t i, const Fence *fence);
  void setFences(const FenceVector &fences);
  void invalidate(Invalidation what, size_t i = 0); // i only for InvalidateRibbon
  const Ribbon *ribbon(size_t i) const;
  Point3D evaluate(const Point2D &uv) const;
  Mesh evaluate(size_t resolution) const;
protected:
  Domain *domain_;
  std::vector<Parameterization *> param_;
  std::vector<Ribbon *> ribbons_;
};

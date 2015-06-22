#pragma once

#include "geometry.hh"

namespace Transfinite
{

  class Domain;
  class Parameterization;
  class Ribbon;

  // Notes:
  // - boundary curve changes => domain || ribbon
  //     (if the domain remains the same, the invalidate routine should invalidate the ribbon)
  // - fence changes => ribbon
  // - ribbon handler changes => ribbon
  // - resolution changes => parameterization

  enum Invalidation
  {
    InvalidateNothing = 0,
    InvalidateDomain = 1,              // => InvalidateParameterization
    InvalidateParameterization = 2,    // => InvalidateBlendFunctions, InvalidateAllRibbons
    InvalidateRibbon = 3,
    InvalidateBlendFunctions = 4,
    InvalidateAllRibbons = 5
  };

  // Notes:
  // - set... functions clear the appropriate cache data

  class Surface
  {
  public:
    Surface(size_t sides);
    virtual ~Surface();
    void setSide(size_t i, BSCurve const *curve);
    void setSides(CurveVector const &curves);
    void setFence(size_t i, Fence const *fence);
    void setFences(FenceVector const &fences);
    void invalidate(Invalidation what, size_t i = 0); // i only for InvalidateRibbon
    Ribbon const *ribbon(size_t i) const;
    Point3D evaluate(Point2D const &uv) const;
    Mesh evaluate(size_t resolution) const;
  protected:
    Domain *domain;
    std::vector<Parameterization *> param;
    std::vector<Ribbon *> ribbons;
  };

} // Transfinite

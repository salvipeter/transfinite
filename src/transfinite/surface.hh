#pragma once

#include "geometry.hh"

class Domain;
class Parameterization;
class Ribbon;

class Surface {
public:
  Surface();
  virtual ~Surface();
  void setGamma(bool use_gamma);
  void setCurve(size_t i, const std::shared_ptr<BSCurve> &curve);
  void setCurves(const CurveVector &curves);
  void setupLoop();
  virtual void update(size_t i);
  virtual void update();
  virtual Point3D eval(const Point2D &uv) const = 0;
  TriMesh eval(size_t resolution) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const = 0;
  DoubleVector blendSideSingular(const Point2DVector &sds) const;
  double gamma(double d) const;

  size_t next(size_t i, size_t j = 1) const { return (i + j) % n_; }
  size_t prev(size_t i, size_t j = 1) const { return (i + n_ - j) % n_; }

  size_t n_;
  std::shared_ptr<Domain> domain_;
  std::shared_ptr<Parameterization> param_;
  std::vector<std::shared_ptr<Ribbon>> ribbons_;
  bool use_gamma_;
};

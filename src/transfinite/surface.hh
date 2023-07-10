#pragma once

#include "curves.hh"

#include <functional>
#include <optional>

namespace Transfinite {

using namespace Geometry;

using NormalFence = std::function<Vector3D(double)>;
  
class Domain;
class Parameterization;
class Ribbon;

class Surface {
public:
  Surface();
  Surface(const Surface &) = default;
  virtual ~Surface();
  Surface &operator=(const Surface &) = default;
  void setGamma(bool use_gamma);
  void setCurve(size_t i, const std::shared_ptr<Curve> &curve);
  void setCurves(const CurveVector &curves);
  virtual void setupLoop();
  double ribbonMultiplier(size_t i) const;
  void setRibbonMultiplier(size_t i, double m);
  std::optional<Vector3D> ribbonHandler(size_t i) const;
  void setRibbonHandler(size_t i, const Vector3D &h);
  void overrideNormalFence(size_t i, const std::shared_ptr<NormalFence> &fence);
  void resetRibbon(size_t i);
  virtual void update(size_t i);
  virtual void update();
  std::shared_ptr<const Domain> domain() const;
  std::shared_ptr<const Parameterization> parameterization() const;
  std::shared_ptr<const Ribbon> ribbon(size_t i) const;
  virtual Point3D eval(const Point2D &uv) const = 0;
  virtual TriMesh eval(size_t resolution) const;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const = 0;
  Point3D cornerCorrection(size_t i, double s1, double s2) const;
  Point3D sideInterpolant(size_t i, double si, double di) const;
  Point3D cornerInterpolant(size_t i, const Point2DVector &sds) const;
  Point3D cornerInterpolantD(size_t i, const Point2DVector &sds) const;
  DoubleVector blendCorner(const Point2DVector &sds) const;
  DoubleVector blendSideSingular(const Point2DVector &sds) const;
  DoubleVector blendCornerDeficient(const Point2DVector &sds) const;

  size_t next(size_t i, size_t j = 1) const { return (i + j) % n_; }
  size_t prev(size_t i, size_t j = 1) const { return (i + n_ - j) % n_; }

  size_t n_;
  std::shared_ptr<Domain> domain_;
  std::shared_ptr<Parameterization> param_;
  std::vector<std::shared_ptr<Ribbon>> ribbons_;

private:
  struct CornerData {
    Point3D point;
    Vector3D tangent1, tangent2, twist1, twist2;
  };

  void updateCorner(size_t i);
  void updateCorners();
  double gamma(double d) const;
  static Vector3D rationalTwist(double u, double v, const Vector3D &f, const Vector3D &g);

  std::vector<CornerData> corner_data_;
  bool use_gamma_;
};

} // namespace Transfinite

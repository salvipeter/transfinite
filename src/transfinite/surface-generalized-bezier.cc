#include "domain-regular.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-dummy.hh"
#include "surface-generalized-bezier.hh"

namespace Transfinite {

using DomainType = DomainRegular;
using ParamType = ParameterizationBarycentric;
using RibbonType = RibbonDummy;

SurfaceGeneralizedBezier::SurfaceGeneralizedBezier() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

SurfaceGeneralizedBezier::~SurfaceGeneralizedBezier() {
}

/*
  This is a trivial implementation of the surface evaluator.
  It is much slower than the one given below,
  as it delegates weight computation to the `weight` function,
  but it is easier to understand. Also see the comments in `weight`.

Point3D
SurfaceGeneralizedBezier::eval(const Point2D &uv) const {
  Point3D surface_point(0,0,0);
  double weight_sum = 0.0;
  for (size_t i = 0; i < n_; ++i) {
    for (size_t k = 0; k < layers_; ++k) {
      for (size_t j = 0; j <= degree_; ++j) {
        double blend = weight(i, j, k, uv);
        surface_point += nets_[i][j][k] * blend;
        weight_sum += blend;
      }
    }
  }
  surface_point += central_cp_ * (1.0 - weight_sum);
  return surface_point;
}
*/

Point3D
SurfaceGeneralizedBezier::eval(const Point2D &uv) const {
  Point3D surface_point(0,0,0);
  Point2DVector sds = param_->mapToRibbons(uv);
  double weight_sum = 0.0;
  for (size_t i = 0; i < n_; ++i) {
    const double &si   = sds[i][0];
    const double &di_1 = sds[prev(i)][1];
    const double &di   = sds[i][1];
    const double &di1  = sds[next(i)][1];
    if (di + di1 < epsilon || di_1 + di < epsilon) {
      if (di + di1 < epsilon)
        return nets_[i][degree_][0];
      else
        return nets_[i][0][0];
    }
    double alpha = di_1 / (di_1 + di);
    double beta  = di1  / (di1  + di);
    DoubleVector bl_s, bl_d;
    bernstein(degree_, si, bl_s);
    bernstein(degree_, di, bl_d);
    for (size_t k = 0; k < layers_; ++k) {
      for (size_t j = 0; j <= degree_; ++j) {
        double blend = bl_s[j] * bl_d[k];
        if (k < 2 && (j < 2 || j > degree_ - 2)) {
          if (j < 2)
            blend *= alpha;
          else
            blend *= beta;
        } else if (j == k || j == degree_ - k)
          blend *= 0.5;
        else if (j < k || j > degree_ - k)
          blend = 0.0;
        surface_point += nets_[i][j][k] * blend;
        weight_sum += blend;
      }
    }
  }
  surface_point += central_cp_ * (1.0 - weight_sum);
  return surface_point;
}

size_t
SurfaceGeneralizedBezier::degree() const {
  return degree_;
}

size_t
SurfaceGeneralizedBezier::layers() const {
  return layers_;
}

void
SurfaceGeneralizedBezier::initNetwork(size_t n, size_t degree) {
  if (n == n_ && degree == degree_)
    return;

  n_ = n;
  degree_ = degree;
  layers_ = (degree_ + 1) / 2;

  // Just allocate a 3-dimensional matrix for control points
  nets_.clear(); nets_.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    ControlNet cn; cn.reserve(degree_ + 1);
    for (size_t j = 0; j <= degree_; ++j) {
      PointVector cs(layers_);
      cn.push_back(cs);
    }
    nets_.push_back(cn);
  }
}

void
SurfaceGeneralizedBezier::setupLoop() {
  CurveVector curves;
  DoubleVector knots;
  PointVector cpts(degree_ + 1);
  knots.insert(knots.end(), degree_ + 1, 0.0);
  knots.insert(knots.end(), degree_ + 1, 1.0);
  for (size_t i = 0; i < n_; ++i) {
    for (size_t j = 0; j <= degree_; ++j)
      cpts[j] = nets_[i][j][0];
    curves.push_back(std::make_shared<BSCurve>(degree_, knots, cpts));
  }
  setCurves(curves);

  Surface::setupLoop();

  if (domain_->update())
    param_->update();
}

Point3D
SurfaceGeneralizedBezier::centralControlPoint() const {
  return central_cp_;
}

void
SurfaceGeneralizedBezier::setCentralControlPoint(const Point3D &p) {
  central_cp_ = p;
}

Point3D
SurfaceGeneralizedBezier::controlPoint(size_t i, size_t j, size_t k) const {
  return nets_[i][j][k];
}

void
SurfaceGeneralizedBezier::setControlPoint(size_t i, size_t j, size_t k, const Point3D &p) {
  nets_[i][j][k] = p;
  if (j < layers_)
    nets_[prev(i)][degree_-k][j] = p;
  else if (degree_ - j < layers_)
    nets_[next(i)][k][degree_-j] = p;
}

double
SurfaceGeneralizedBezier::weight(size_t i, size_t j, size_t k, const Point2D &uv) const {
  // Points outside the diagonals have 0 weight (except the (0,1) point)
  if (k >= 2 && (j < k || j > degree_ - k))
    return 0.0;

  // Otherwise we need the local parameters
  Point2DVector sds = param_->mapToRibbons(uv);
  const double &si   = sds[i][0];
  const double &di_1 = sds[prev(i)][1];
  const double &di   = sds[i][1];
  const double &di1  = sds[next(i)][1];

  // At the corners give back 0.5 for the corner point, 0.0 otherwise
  // (the rational weights used below have a singularity at these points)
  if (di + di1 < epsilon)
    return j == degree_ && k == 0 ? 0.5 : 0.0;
  if (di_1 + di < epsilon)
    return j == 0 && k == 0 ? 0.5 : 0.0;

  // Basic blend is computed by Bernstein functions
  double blend = bernstein(j, degree_, si) * bernstein(k, degree_, di);

  // At the 2x2 corners use rational weights
  if (k < 2 && (j < 2 || j > degree_ - 2)) {
    if (j < 2)
      return blend * di_1 / (di_1 + di);
    else
      return blend * di1 / (di1  + di);
  }

  // On the diagonals use half of the weight
  if (j == k || j == degree_ - k)
    return blend * 0.5;

  // Between the diagonals use the blend as it is
  return blend;
}

std::shared_ptr<Ribbon>
SurfaceGeneralizedBezier::newRibbon() const {
  return std::make_shared<RibbonType>();
}

void
SurfaceGeneralizedBezier::bernstein(size_t n, double u, DoubleVector &coeff) {
  coeff.clear(); coeff.reserve(n + 1);
  coeff.push_back(1.0);
  double u1 = 1.0 - u;
  for (size_t j = 1; j <= n; ++j) {
    double saved = 0.0;
    for (size_t k = 0; k < j; ++k) {
      double  tmp = coeff[k];
      coeff[k] = saved + tmp * u1;
      saved = tmp * u;
    }
    coeff.push_back(saved);
  }
}

double
SurfaceGeneralizedBezier::bernstein(size_t i, size_t n, double u) {
  DoubleVector tmp(n + 1, 0.0);
  tmp[n-i] = 1.0;
  const double u1 = 1.0 - u;
  for (size_t k = 1; k <= n; ++k)
    for (size_t j = n; j >= k; --j)
      tmp[j] = tmp[j-1] * u + tmp[j] * u1;
  return tmp[n];
}

} // namespace Transfinite

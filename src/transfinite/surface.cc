#include <algorithm>

#include "domain.hh"
#include "parameterization.hh"
#include "ribbon.hh"
#include "surface.hh"
#include "utilities.hh"

// For the 4-sided central split fit
#include "parameterization-bilinear.hh"

namespace Transfinite {

Surface::Surface()
  : n_(0), use_gamma_(true) {
}

Surface::~Surface() {
}

void
Surface::setGamma(bool use_gamma) {
  use_gamma_ = use_gamma;
}

void
Surface::setCurve(size_t i, const std::shared_ptr<BSCurve> &curve) {
  if (n_ <= i) {
    ribbons_.resize(i + 1);
    n_ = i + 1;
  }
  ribbons_[i] = newRibbon();
  ribbons_[i]->setCurve(curve);
  domain_->setSide(i, curve);
}

void
Surface::setCurves(const CurveVector &curves) {
  ribbons_.clear();
  ribbons_.reserve(curves.size());
  for (const auto &c : curves) {
    ribbons_.push_back(newRibbon());
    ribbons_.back()->setCurve(c);
  }
  domain_->setSides(curves);
  n_ = curves.size();
}

void
Surface::setupLoop() {
  // Tasks:
  // - propagate adjacency information
  // - normalize curves
  // - reverse curves when needed (and normalize once again, for safety)
  for (size_t i = 0; i < n_; ++i)
    ribbons_[i]->curve()->normalize();
  for (size_t i = 0; i < n_; ++i) {
    std::shared_ptr<Ribbon> rp = ribbons_[prev(i)], rn = ribbons_[next(i)];
    ribbons_[i]->setNeighbors(rp, rn);
    if (i == 0) {
      Point3D r_start = ribbons_[i]->curve()->eval(0.0);
      Point3D r_end = ribbons_[i]->curve()->eval(1.0);
      Point3D rn_start = rn->curve()->eval(0.0);
      Point3D rn_end = rn->curve()->eval(1.0);
      double end_to_start = (r_end - rn_start).norm();
      double end_to_end = (r_end - rn_end).norm();
      double start_to_start = (r_start - rn_start).norm();
      double start_to_end = (r_start - rn_end).norm();
      if (std::min(start_to_start, start_to_end) < std::min(end_to_start, end_to_end)) {
        ribbons_[i]->curve()->reverse();
        ribbons_[i]->curve()->normalize();
      }
    } else {
      Point3D r_start = ribbons_[i]->curve()->eval(0.0);
      Point3D r_end = ribbons_[i]->curve()->eval(1.0);
      Point3D rp_end = rp->curve()->eval(1.0);
      if ((r_end - rp_end).norm() < (r_start - rp_end).norm()) {
        ribbons_[i]->curve()->reverse();
        ribbons_[i]->curve()->normalize();
      }
    }
  }
}

std::optional<Vector3D>
Surface::ribbonHandler(size_t i) const {
  return ribbons_[i]->handler();
}

void
Surface::setRibbonHandler(size_t i, const Vector3D &h) {
  ribbons_[i]->setHandler(h);
}

double
Surface::ribbonMultiplier(size_t i) const {
  return ribbons_[i]->multiplier();
}

void
Surface::setRibbonMultiplier(size_t i, double m) {
  ribbons_[i]->setMultiplier(m);
}

void
Surface::resetRibbon(size_t i) {
  ribbons_[i]->reset();
}

void
Surface::update(size_t i) {
  if (domain_->update())
    param_->update();
  ribbons_[i]->update();
  updateCorner(prev(i));
  updateCorner(i);
}

void
Surface::update() {
  if (domain_->update())
    param_->update();
  for (auto &r : ribbons_)
    r->update();
  updateCorners();
}

std::shared_ptr<const Domain>
Surface::domain() const {
  return domain_;
}

std::shared_ptr<const Parameterization>
Surface::parameterization() const {
  return param_;
}

std::shared_ptr<const Ribbon>
Surface::ribbon(size_t i) const {
  return ribbons_[i];
}

TriMesh
Surface::eval(size_t resolution) const {
  TriMesh mesh = domain_->meshTopology(resolution);
  Point2DVector uvs = domain_->parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [this](const Point2D &uv) { return eval(uv); });
  mesh.setPoints(points);
  return mesh;
}

Point3D
Surface::cornerCorrection(size_t i, double s1, double s2) const {
  // Assumes that both s1 and s2 are 0 at the corner,
  // s1 increases towards corner (i-1), and s2 towards corner (i+1).
  s1 = inrange(0, gamma(s1), 1);
  s2 = inrange(0, gamma(s2), 1);
  return corner_data_[i].point
    + corner_data_[i].tangent1 * s1
    + corner_data_[i].tangent2 * s2
    + rationalTwist(s1, s2, corner_data_[i].twist2, corner_data_[i].twist1) * s1 * s2;
}

Point3D
Surface::sideInterpolant(size_t i, double si, double di) const {
  si = inrange(0, si, 1);
  di = std::max(gamma(di), 0.0);
  return ribbons_[i]->eval(Point2D(si, di));
}

Point3D
Surface::cornerInterpolant(size_t i, const Point2DVector &sds) const {
  double si = sds[i][0], si1 = sds[next(i)][0];
  return sideInterpolant(i, si, si1) + sideInterpolant(next(i), si1, 1.0 - si)
    - cornerCorrection(i, 1.0 - si, si1);
}

Point3D
Surface::cornerInterpolantD(size_t i, const Point2DVector &sds) const {
  double di = sds[i][1], di1 = sds[next(i)][1];
  return sideInterpolant(i, 1.0 - di1, di) + sideInterpolant(next(i), di, di1)
    - cornerCorrection(i, di1, di);
}

DoubleVector
Surface::blendCorner(const Point2DVector &sds) const {
  DoubleVector blf; blf.reserve(n_);

  size_t close_to_boundary = 0;
  for (const auto &sd : sds) {
    if (sd[1] < epsilon)
      ++close_to_boundary;
  }

  if (close_to_boundary > 0) {
    for (size_t i = 0; i < n_; ++i) {
      size_t ip = next(i);
      if (close_to_boundary > 1)
        blf.push_back(sds[i][1] < epsilon && sds[ip][1] < epsilon ? 1.0 : 0.0);
      else if (sds[i][1] < epsilon) {
        double tmp = std::pow(sds[ip][1], -2);
        blf.push_back(tmp / (tmp + std::pow(sds[prev(i)][1], -2)));
      } else if (sds[ip][1] < epsilon) {
        double tmp = std::pow(sds[i][1], -2);
        blf.push_back(tmp / (tmp + std::pow(sds[next(ip)][1], -2)));
      } else
        blf.push_back(0.0);
    }
  } else {
    double denominator = 0.0;
    for (size_t i = 0; i < n_; ++i) {
      blf.push_back(std::pow(sds[i][1] * sds[next(i)][1], -2));
      denominator += blf.back();
    }
    std::transform(blf.begin(), blf.end(), blf.begin(),
                   [denominator](double x) { return x / denominator; });
  }

  return blf;
}

DoubleVector
Surface::blendSideSingular(const Point2DVector &sds) const {
  DoubleVector blf; blf.reserve(n_);

  size_t close_to_boundary = 0;
  for (const auto &sd : sds) {
    if (sd[1] < epsilon)
      ++close_to_boundary;
  }

  if (close_to_boundary > 0) {
    double blend_val = 1.0 / close_to_boundary;
    for (const auto &sd : sds)
      blf.push_back(sd[1] < epsilon ? blend_val : 0.0);
  } else {
    double denominator = 0.0;
    for (const auto &sd : sds) {
      blf.push_back(std::pow(sd[1], -2));
      denominator += blf.back();
    }
    std::transform(blf.begin(), blf.end(), blf.begin(),
                   [denominator](double x) { return x / denominator; });
  }

  return blf;
}

DoubleVector
Surface::blendCornerDeficient(const Point2DVector &sds) const {
  DoubleVector blf; blf.reserve(n_);
  for (size_t i = 0; i < n_; ++i) {
    size_t ip = next(i);
    if (sds[i][1] < epsilon && sds[ip][1] < epsilon) {
      blf.push_back(1.0);
      continue;
    }
    blf.push_back((sds[ip][1] * hermite(0, 1.0 - sds[i][0]) * hermite(0, sds[i][1] ) +
                   sds[i][1]  * hermite(0,    sds[ip][0]  ) * hermite(0, sds[ip][1])) /
                  (sds[i][1] + sds[ip][1]));
  }
  return blf;
}

void
Surface::updateCorner(size_t i) {
  static const double step = 1.0e-4;
  size_t ip = next(i);

  VectorVector der;
  Vector3D d1, d2;
  corner_data_[i].point = ribbons_[i]->curve()->eval(1.0, 1, der);
  corner_data_[i].tangent1 = -der[1];
  ribbons_[ip]->curve()->eval(0.0, 1, der);
  corner_data_[i].tangent2 = der[1];
  d1 = ribbons_[i]->crossDerivative(1.0);
  d2 = ribbons_[i]->crossDerivative(1.0 - step);
  corner_data_[i].twist1 = (d2 - d1) / step;
  d1 = ribbons_[ip]->crossDerivative(0.0);
  d2 = ribbons_[ip]->crossDerivative(step);
  corner_data_[i].twist2 = (d2 - d1) / step;
}

void
Surface::updateCorners() {
  corner_data_.resize(n_);
  for (size_t i = 0; i < n_; ++i)
    updateCorner(i);
}

double
Surface::gamma(double d) const {
  if (use_gamma_)
    return d / (2.0 * d + 1.0);
  return d;
}

Vector3D
Surface::rationalTwist(double u, double v, const Vector3D &f, const Vector3D &g) {
  if (std::abs(u + v) < epsilon)
    return Vector3D(0,0,0);
  return (f * u + g * v) / (u + v);
}

#ifndef NO_SURFACE_FIT

namespace {
  double knot_snap_tolerance;

  double snapToKnot(double x, const BSCurve &c) {
    DoubleVector knots = c.knotVector();
    for (double i : knots)
      if (fabs(i - x) < knot_snap_tolerance)
        return i;
    return x;
  }

  void unifyKnots(BSCurve &c1, BSCurve &c2) {
    bool changed;
    DoubleVector k1 = c1.knotVector(), k2 = c2.knotVector();
    do {
      changed = false;
      for (DoubleVector::const_iterator i = k1.begin(), j = k2.begin(),
             ie = k1.end(), je = k2.end(); i != ie && j != je; ++i, ++j) {
        if (fabs(*i - *j) < knot_snap_tolerance && *i != 1.0 && *j != 1.0) {
          // Treat these knots as if they were the same.
          // Better than inserting very small intervals in both.
          // Note, that we do not have this option at the ends of the knot vectors. 
        } else if (*i < *j) {
          c2.insertKnot(*i);
          k2.insert(j, *i);
          changed = true;
          break;
        } else if (*i > *j) {
          c1.insertKnot(*j);
          k1.insert(i, *j);
          changed = true;
          break;
        }
      }
    } while (changed);
  }
}

std::vector<BSSurface>
Surface::fitCentralSplit(double fit_tol, double knot_snapping_tol, size_t sampling_density) const {
  // Plan
  // ----
  // when n == 4:
  // 1. create 4 edge curves: edge_curves[0..3]
  // 2. unify the opposite knot vectors
  // 3. fit
  // when n != 4:
  // 1. create n dividing curves: dividing_curves[0..n-1]
  // 2. for all subsurfaces:
  //    3. create 2 half-edge curves: edge_curves[0..1]
  //    4. determine the appropriate dividing curves
  //    5. unify the opposite knot vectors
  //    6. fit

  knot_snap_tolerance = knot_snapping_tol;
  std::vector<BSSurface> surfaces;

  if (n_ == 4) {
    // 1. Generate edge curves
    BSCurve edge_curves[4];
    for (size_t i = 0; i < 4; ++i) {
      std::shared_ptr<Ribbon> r = ribbons_[i];
      edge_curves[i] = *r->curve();
      if (i > 1) {
        edge_curves[i].reverse();
        edge_curves[i].normalize();
      }
    }

    // 2. Unify the opposite knot vectors
    unifyKnots(edge_curves[0], edge_curves[2]);
    unifyKnots(edge_curves[1], edge_curves[3]);

    // 3. Fit
    SurfaceFitter fitter;
    fitter.setDegreeU(3);
    fitter.setDegreeV(3);

    ParameterizationBilinear param_bilinear;
    param_bilinear.setDomain(domain_);
    param_bilinear.update();

    // Problem: we want to pass parameter values on a [0,1]x[0,1] square,
    //          but the actual domain can be anything, so we do a bilinear mapping,
    //          based on the first ribbon
    size_t &u_step = sampling_density;
    Point2D center = domain_->center();
    Point3D point = eval(center);
    Point2D bilinear = param_bilinear.mapToRibbon(0, center);
    fitter.newPointGroup(fit_tol);
    fitter.addParamPoint(bilinear, point);

    for (size_t j = 1; j <= u_step; ++j) {
      double j_coeff = (double)j / u_step;
      for (size_t k = 0; k < 4; k++) {
        for (size_t i = 0; i < j; i++) {
          Point2D p = center * (1.0 - j_coeff) + domain_->edgePoint(k, (double)i/j) * j_coeff;
          point = eval(p);
          bilinear = param_bilinear.mapToRibbon(0, p);
          fitter.addParamPoint(bilinear, point);
        }
      }
    }
        
    fitter.setKnotVectorU(edge_curves[0].knotVector());
    fitter.setKnotVectorV(edge_curves[1].knotVector());

    const size_t nr_ctrl[] = { edge_curves[0].nrControlPoints(),
                               edge_curves[1].nrControlPoints() };

    // Setup control point constraints
    for (size_t i = 0; i < 2; ++i) {
      size_t end = nr_ctrl[1-i] - 1;
      const BSCurve &e1 = edge_curves[i], &e2 = edge_curves[i+2];
      for (size_t j = 0; j < nr_ctrl[i]; ++j) {
        if (i == 0) {
          fitter.addControlPoint(j,   0,   e1.controlPoint(j));
          fitter.addControlPoint(j,   end, e2.controlPoint(j));
        } else {
          fitter.addControlPoint(end, j,   e1.controlPoint(j));
          fitter.addControlPoint(0,   j,   e2.controlPoint(j));
        }
      }
    }

    fitter.fit();
    surfaces.push_back(fitter.surface());
  } else {
    surfaces.reserve(n_);
    Point2D center = domain_->center();
    size_t &u_step = sampling_density;

    // 0. Generate dividing curves (from the side to the center)
    std::vector<BSCurve> dividing_curves; dividing_curves.reserve(n_);
    DoubleVector edge_midpoints; edge_midpoints.resize(n_);
    Point2DVector edge_mid_u; edge_mid_u.resize(n_);
    for (size_t i = 0; i < n_; ++i) {
      CurveFitter bf;
      bf.setDegree(3);
      bf.newPointGroup(fit_tol);
      double edge_half;
      {
        std::shared_ptr<Ribbon> r = ribbons_[i];
        BSCurve c = *r->curve();
        double mid = snapToKnot(0.5, c);
        edge_midpoints[i] = mid;
        edge_half = mid;
      }
      Point2D u0 = domain_->edgePoint(i, edge_half);
      edge_mid_u[i] = u0;
      for (size_t j = 0; j <= u_step; j++) {
        double u = (double)j / u_step;
        Point2D param = u0 * (1.0 - u) + center * u;
        Point3D point = eval(param);
        bf.addParamPoint(u, point);
      }
      size_t nr_cpts = 5;       // kutykurutty
      bf.setNrControlPoints(nr_cpts);
      bf.addControlPoint(0, eval(u0));
      bf.addControlPoint(nr_cpts - 1, eval(center));
      bf.fit();
      dividing_curves.push_back(bf.curve());
    }

    for (size_t i = 0; i < n_; ++i) {
      // 1. Generate edge curves
      BSCurve edge_curves[2];
      for (size_t j = 0; j < 2; ++j) {
        size_t ip = next(i, j);
        std::shared_ptr<Ribbon> r = ribbons_[ip];
        edge_curves[j] = *r->curve();
        double mid = edge_midpoints[ip];
        if (j == 0) {
          double end = snapToKnot(1.0, edge_curves[j]);
          edge_curves[j].trim(mid, end);
        } else {
          double start = snapToKnot(0.0, edge_curves[j]);
          edge_curves[j].trim(start, mid);
        }
        edge_curves[j].normalize();
      }

      // 2. Unify the opposite knot vectors
      BSCurve div_curves[2];
      div_curves[0] = dividing_curves[next(i)];
      div_curves[0].reverse();
      div_curves[0].normalize();
      div_curves[1] = dividing_curves[i];
      div_curves[1].normalize();
      unifyKnots(edge_curves[0], div_curves[0]);
      unifyKnots(edge_curves[1], div_curves[1]);

      // 3. Fit
      SurfaceFitter fitter;
      fitter.setDegreeU(3);
      fitter.setDegreeV(3);
      fitter.newPointGroup(fit_tol);

      Point2D u1 = edge_mid_u[i];
      Point2D v1 = edge_mid_u[next(i)];
      Point2D u_end = domain_->vertices()[i];

      size_t v_step = u_step;

      for (size_t j = 0; j <= u_step; j++) {
        double u = (double)j / u_step;
        Point2D u_v0 = u1 * (1.0 - u) + u_end * u;
        Point2D u_v1 = center * (1.0 - u) + v1 * u;
        for (size_t k = 0; k <= v_step; k++) {
          double v = (double)k / v_step;
          Point2D param = u_v0 * (1.0 - v) + u_v1 * v;
          Point3D point = eval(param);
          fitter.addParamPoint(Point2D(u, v), point);
        }
      }

      fitter.setKnotVectorU(edge_curves[0].knotVector());
      fitter.setKnotVectorV(edge_curves[1].knotVector());

      size_t const nr_ctrl[] = { edge_curves[0].nrControlPoints(),
                                 edge_curves[1].nrControlPoints() };

      // Setup control point constraints
      for (size_t k = 0; k < 2; ++k) {
        size_t end = nr_ctrl[1-k] - 1;
        const BSCurve &e = edge_curves[k], &d = div_curves[k];
        for (size_t j = 0; j < nr_ctrl[k]; ++j) {
          if (k == 0) {
            fitter.addControlPoint(j,   0,   e.controlPoint(j));
            fitter.addControlPoint(j,   end, d.controlPoint(j));
          } else {
            fitter.addControlPoint(end, j,   e.controlPoint(j));
            fitter.addControlPoint(0,   j,   d.controlPoint(j));
          }
        }
      }

      fitter.fit();
      surfaces.push_back(fitter.surface());
    }
  }
  return surfaces;
}

BSSurface
Surface::fitTrimmed(double fit_tol, size_t resolution, size_t max_cpts_u, size_t max_cpts_v,
                    double curvature_weight, double oscillation_weight) const {
  Point2DVector uvs = domain_->parameters(resolution);

  SurfaceFitter fitter;

  // Add all points
  fitter.newPointGroup(fit_tol);
  for (const auto &uv : uvs)
    fitter.addParamPoint(uv, eval(uv));

  // Add boundary points as separate point groups
  for (size_t i = 0; i < n_; ++i) {
    fitter.newPointGroup(fit_tol);
    for (size_t j = 0; j <= resolution; ++j) {
      double t = (double)j / resolution;
      Point2D uv = domain_->edgePoint(i, t);
      fitter.addParamPoint(uv, eval(uv));
    }
  }

  fitter.setDegreeU(3);
  fitter.setDegreeV(3);
  fitter.setNrControlPointsU(6);
  fitter.setNrControlPointsV(6);
  fitter.setMaxNrControlPointsU(max_cpts_u);
  fitter.setMaxNrControlPointsV(max_cpts_v);
  fitter.setCurvatureWeight(curvature_weight);
  fitter.setOscillationWeight(oscillation_weight);

  fitter.fitWithCarrierSurface();

  // Save trimming information with the surface
  BSSurface s = fitter.surface(); s.curves_.reserve(n_);
  std::transform(ribbons_.begin(), ribbons_.end(), std::back_inserter(s.curves_),
                 [](const std::shared_ptr<Ribbon> &r) { return r->curve(); });

  // Fit curves on the trim parameters
  for (size_t i = 1; i <= n_; ++i) {
    Point2DVector params = fitter.parameters(i);
    PointVector params3d; params3d.reserve(params.size());
    std::transform(params.begin(), params.end(), std::back_inserter(params3d),
                   [](const Point2D &p) { return Point3D(p[0], p[1], 0.0); });
    CurveFitter bf;
    bf.setDegree(3);
    bf.newPointGroup(fit_tol);
    for (size_t j = 0; j <= resolution; ++j) {
      double t = (double)j / resolution;
      bf.addParamPoint(t, params3d[j]);
    }
    size_t nr_cpts = 5;       // kutykurutty
    bf.setNrControlPoints(nr_cpts);
    bf.addControlPoint(0, params3d[0]);
    bf.addControlPoint(nr_cpts - 1, params3d[resolution]);
    bf.fit();
    s.param_curves_.push_back(std::make_shared<BSCurve>(bf.curve()));
  }

  return s;
}

#endif  // NO_SURFACE_FIT

} // namespace Transfinite

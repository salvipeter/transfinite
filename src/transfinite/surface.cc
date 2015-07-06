#include <algorithm>

#include "domain.hh"
#include "parameterization.hh"
#include "ribbon.hh"
#include "surface.hh"

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
  s1 = std::min(std::max(gamma(s1), 0.0), 1.0);
  s2 = std::min(std::max(gamma(s2), 0.0), 1.0);
  return corner_data_[i].point
    + corner_data_[i].tangent1 * s1
    + corner_data_[i].tangent2 * s2
    + rationalTwist(s1, s2, corner_data_[i].twist2, corner_data_[i].twist1) * s1 * s2;
}

Point3D
Surface::sideInterpolant(size_t i, double si, double di) const {
  si = std::min(std::max(si, 0.0), 1.0);
  di = std::max(gamma(di), 0.0);
  return ribbons_[i]->eval(Point2D(si, di));
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

double
Surface::blendHermite(double x) {
  double x2 = x * x;
  return 2.0 * x * x2 - 3.0 * x2 + 1.0;
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
    fitter.setTolerance(fit_tol);
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
      bf.setTolerance(fit_tol);
      bf.setDegree(3);

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
      fitter.setTolerance(fit_tol);
      fitter.setDegreeU(3);
      fitter.setDegreeV(3);

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
Surface::fitTrimmed(...) const {
  bool is_4_sided = n_ == 4;
  bool smart_parametrization =
    O("Smart Parametrization for All") || O("Smart Parametrization for 4 sided");
  bool use_simple_rotation = O("Use Only Rotation");
  double extra_rotation = O("Extra Rotation Angle") * M_PI / 180.0;

  int max_ribbon_id = -1;
  double ribbon_3D_ratio = 1.0;
  if (smart_parametrization) {
    // Find longest ribbon
    double max_ribbon_length = -1.0;
    for (size_t i = 0; i < n_; ++i) {
      double ribbon_length = ribbons_[i]->curve()->arcLength(0.0, 1.0);
      if (ribbon_length > max_ribbon_length) {
        max_ribbon_length = ribbon_length;
        max_ribbon_id = i;
      }
    }
    if (max_ribbon_id == -1 || max_ribbon_length < epsilon) {
      // Error - longest ribbon not found
      smart_parametrization = false;
    } else {
      double ribbon_width =
        (ribbons_[max_ribbon_id]->crossDerivative(0.0).norm()
         + ribbons_[max_ribbon_id]->crossDerivative(1.0).norm())
        / 2.0;
      if (ribbon_width > epsilon)
        ribbon_3D_ratio = ribbon_width / max_ribbon_length;
    }
  }

  SurfaceFitter fitter;
  fitter.setTolerance(tol);

  Point2DVector uvs = domain_->parameters(resolution);
  int size = uvs.size();
  for (int i = 0; i < size; ++i) {
    Point2D param = uvs[i];
    if (use_simple_rotation)
      param = param_.mapToRibbon(max_ribbon_id, param, extra_rotation);
    else if (smart_parametrization) {
      param = param_.mapToRibbon(max_ribbon_id, param);
      if (O("Extra Smart Parametrization for All"))
        param[1] *= ribbon_3D_ratio;
    }
    fitter.addParamPoint(param, eval(uvs[i]));
  }

  // TODO: add second point group (tol) "border point group"

  int border_curve_sample_density = O("Border Curve Sample Density");
  int mesh_resolution = O("Surface Triangles");
  int num_inner_brd_points = O("Ratio of Boundary/Mesh Resolution");
  border_curve_sample_density = mesh_resolution * num_inner_brd_points;
  for (size_t i = 0; i < n_; ++i) {
    for (int j = 0; j <= border_curve_sample_density; ++j) {
      double u = (double)j / border_curve_sample_density;
      Point2D ribbon_param(u, 0.0);
      Point2D param = param_.mapToRibbon(i, ribbon_param);

      param[0] = std::round(1.0e5 * param[0]) * 1.0e-5;
      param[1] = std::round(1.0e5 * param[1]) * 1.0e-5;

      Point3D point = ribbons_[i]->curve()->eval(ribbon_param[0]);
      if (use_simple_rotation)
        param = param_.mapToRibbon(max_ribbon_id, param, extra_rotation);
      else if (smart_parametrization) {
        param = param_.mapToRibbon(max_ribbon_id, param);
        if (O("Extra Smart Parametrization for All")) {
          param[1] *= ribbon_3D_ratio;
        }
      }
      fitter.addParamPoint(param, point); // border point group
    }
  }

  bool extend_ribbons = !smart_parametrization && O("Extend Ribbons");
  if (extend_ribbons) {
    // add third point group (tol * 5.0) "exterior point group"
    double ribbon_extension = O("Ribbon Extension Ratio");
    double u_step = 1.0 / O("Ribbon Sample Density");
    double v_step = u_step;
    double domain_size = 1.0 + O("Domain Enlargement");

    for (size_t i = 0; i < n_; ++i) {
      for (double u = 0.0; u <= 1.0; u += u_step) {
        for (double v = 0.0; v >= -ribbon_extension; v -= v_step) {
          Point2D ribbon_param(u,v);
          Point2D param = param_.mapToRibbon(i, ribbon_param);
          if (std::abs(param[0]) <= domain_size && std::abs(param[1]) <= domain_size) {
            Point3D point = ribbons_[i]->eval(ribbon_param);
            fitter.addParamPoint(param, point); // exterior point group
          }
        }
      }
    }

    // points from triangles between ribbons
    bool extend_triangles = O("Extend Triangles");
    if (extend_triangles) {
      int triangle_sample_density = O("Ribbon Triangle Sample Density");

      for (size_t i = 0; i < n_; ++i) {
        const Point2D &a_param = domain_.vertices()[i];
        Point3D a = ribbons[i]->curve()->eval(0.0);

        Point2D b_ribbon_param(0.0, -ribbon_extension);
        Point2D b_param = param_.mapToRibbon(i, b_ribbon_param);
        Point3D b = ribbons_[i]->eval(b_ribbon_param);

        Point2D c_ribbon_param(1.0, -ribbon_extension);
        Point2D c_param = param_.mapToRibbon(prev(i), c_ribbon_param);
        Point3D c = ribbons_[prev(i)]->eval(c_ribbon_param);

        for (int j = 1; j <= triangle_sample_density; ++j) {
          double j_coeff = (double)j / triangle_sample_density;
          Point2D ab_param = a_param * (1.0 - j_coeff) + b_param * j_coeff;
          Point3D ab = a * (1.0 - j_coeff) + b * j_coeff;
          Point2D ac_param = a_param * (1.0 - j_coeff) + c_param * j_coeff;
          Point3D ac = a * (1.0 - j_coeff) + c * j_coeff;

          for (int k = 0; k <= j; ++k) {
            double k_coeff = (double)k / j;

            Point2D abc_param = ab_param * (1.0 - k_coeff) + ac_param * k_coeff;

            if (std::abs(abc_param[0]) <= domain_size && std::abs(abc_param[1]) <= domain_size) {
              S::PointType abc = ab * (1.0 - k_coeff) + ac * k_coeff;
              fitter.addParamPoint(abc_param, abc); // exterior point group
            }
          }
        }
      }
    }
  }

  fitter.setDegreeU(3);
  fitter.setDegreeV(3);
  fitter.setNrControlPointsU(6);
  fitter.setNrControlPointsV(6);
  fitter.setOutlierPercentage(O("Fitting Outlier Percentage"));
  fitter.setLocalOutlierPercentages(O("Fitting Local Outliers"));

  int minNrControlPointsU = O("Minimum in U");
  int minNrControlPointsV = O("Minimum in V");
  int maxNrControlPointsU = O("Maximum in U");
  int maxNrControlPointsV = O("Maximum in V");
  if (minNrControlPointsU)
    fitter.setNrControlPointsU(minNrControlPointsU);
  if (minNrControlPointsV)
    fitter.setNrControlPointsV(minNrControlPointsV);
  if (maxNrControlPointsU)
    fitter.setMaxNrControlPointsU(maxNrControlPointsU);
  if (maxNrControlPointsV)
    fitter.setMaxNrControlPointsV(maxNrControlPointsVb);
  if (maxNrControlPointsU && maxNrControlPointsV)
    fitter.setMaxNrControlPoints(maxNrControlPointsU * maxNrControlPointsV);

  bool lsq_dist_method = O("Minimize Dist");
  bool lsq_ctrl_dist_method = false;
  if (lsq_dist_method || !lsq_ctrl_dist_method) {
    double weight = O("Weight for Dist");
    BSSF_FN_LsqDistance<S::PointType>* functional = new BSSF_FN_LsqDistance<S::PointType>();
    functional->Weight() = weight;
    fitter.AddFunctional(functional);
  }
  if (O("Minimize Curvature")) {
    Real weight = O("Weight for Curvature");
    fitter.AddFunctional(new BSSF_FN_Curvature<S::PointType>(1.0, weight));
  }
  if (O("Minimize Oscillation")) {
    Real weight = O("Weight for Oscillation");
    fitter.AddFunctional(new BSSF_FN_Oscillation<S::PointType>(1.0, weight));
  }
  if (O("Minimize Area")) {
    Real weight = O("Weight for Area");
    fitter.AddFunctional(new BSSF_FN_Area<S::PointType>(1.0, weight));
  }

  BSSF_SmoothnessFunctional<S::PointType>::DirectionType smoothness_direction =
    (BSSF_SmoothnessFunctional<S::PointType>::DirectionType) (O("Smoothness Direction") - 1);
  int smoothness_method = O("Smoothness Method");
  switch (smoothness_method) {
  case 2:
    fitter.AddFunctional(new BSSF_FN_Curvature<S::PointType>(smoothness_direction));
    break;
  case 3:
    fitter.AddFunctional(new BSSF_FN_Oscillation<S::PointType>(smoothness_direction));
    break;
  case 4:
    fitter.AddFunctional(new BSSF_FN_Area<S::PointType>(smoothness_direction));
    break;
  default:
    break;
  }

  int knot_insertion = O("Knot Insertion Strategy");
  switch (knot_insertion) {
  case 2:
    fitter.SetKnotInserter(new BSSF_KI_LargestSummedDev<S::PointType>());
    break;
  case 3:
    fitter.SetKnotInserter(new BSSF_KI_LargestSummedDevDir<S::PointType>());
    break;
  case 4:
    fitter.SetKnotInserter(new BSSF_KI_Triangulation<S::PointType>());
    break;
  case 5:
    fitter.SetKnotInserter(new BSSF_KI_FlippingTriangles<S::PointType>());
    break;
  default:
    break;
  }

  fitter.OptimizeParameters() = O("Optimize Parameters");
  fitter.fit();

  return fitter.surface();
}

#endif  // NO_SURFACE_FIT

} // namespace Transfinite

#include <algorithm>

#include "domain.hh"
#include "parameterization.hh"
#include "ribbon.hh"
#include "surface.hh"

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

  double snapToKnot(double x, S::CurveType const &c) {
    S::CurveType::KnotVectorType const &knots = c.theKnotVector();
    double tolerance = knot_snap_tolerance * fabs(knots.back() - knots.front()); 
    for (S::CurveType::KnotVectorType::const_iterator i = knots.begin(), ie = knots.end();
         i != ie; ++i)
      if (fabs(*i - x) < tolerance)
        return *i;
    return x;
  }

  void unifyKnots(S::CurveType &c1, S::CurveType &c2) {
    // Assumes that both curves are parameterized in [0,1]
    bool changed;
    S::CurveType::KnotVectorType const &k1 = c1.theKnotVector();
    S::CurveType::KnotVectorType const &k2 = c2.theKnotVector();
    do {
      changed = false;
      for (S::CurveType::KnotVectorType::const_iterator i = k1.begin(), j = k2.begin(),
             ie = k1.end(), je = k2.end(); i != ie && j != je; ++i, ++j) {
        if (fabs(*i - *j) < knot_snap_tolerance && *i != 1.0 && *j != 1.0) {
          // Treat these knots as if they were the same.
          // Better than inserting very small intervals in both.
          // Note, that we do not have this option at the ends of the knot vectors. 
        } else if (*i < *j) {
          c2.InsertKnot(*i);
          changed = true;
          break;
        } else if (*i > *j) {
          c1.InsertKnot(*j);
          changed = true;
          break;
        }
      }
    } while (changed);
  }
}

void
Surface::fitCentralSplit(double knot_snapping_tol, size_t sampling_density) const {
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

  if (n == 4) {
    // 1. Generate edge curves
    S::CurveType edge_curves[4];
    for (size_t i = 0; i < 4; ++i) {
      SRibbon *r = Ribbons()[i];
      edge_curves[i] = *r->FreeCurve()->Curve();
      double start = r->Low();
      double end = r->Up();
      if (start > end) {
        Real const lower = edge_curves[i].LowerParam();
        Real const upper = edge_curves[i].UpperParam();
        edge_curves[i].Reparametrize((end-lower)+start, (start-upper)+end);
        std::swap(start, end);
      }
      start = snapToKnot(start, edge_curves[i]);
      end = snapToKnot(end, edge_curves[i]);

      edge_curves[i].Trim(start, end);
      if (i > 1)
        edge_curves[i].Reverse();
      edge_curves[i].Reparametrize(0.0, 1.0);
    }

    // 2. Unify the opposite knot vectors
    unifyKnots(edge_curves[0], edge_curves[2]);
    unifyKnots(edge_curves[1], edge_curves[3]);

    // 3. Fit
    SurfaceFitter::FitterType& fitter = sf->AddFitter(object_id, 0);

    BSS_Fitter<S::PointType>::ParamPointSetReference point_set = fitter.ParamPointSet();
    point_set.push_back(BSS_Fitter<S::PointType>::ParamPointGroupType(tol));
    BSS_Fitter<S::PointType>::ParamPointGroupReference pg(point_set.front());

    // Problem: we want to pass parameter values on a [0,1]x[0,1] square,
    //          but the actual domain can be anything, so we do a bilinear mapping,
    //          based on the first ribbon
    S::PointType point;
    int u_step = sampling_density;
    S::Point2D const &centre = domain.Center();
    if (GetSurfacePoint(centre, point)) {
      S::Point2D const bilinear =
        domain.MapNSidedDomainTo4SidedDomain(centre, 0, SDomain::BilinearSweep);
      pg.push_back(BSS_Fitter<S::PointType>::ParamPointType(point, bilinear));
      if (!parametrization_check)
        fitted_points.push_back(point);
    }

    for (int j = 1; j <= u_step; ++j) {
      double const j_coeff = (double)j / u_step;
      for (size_t k = 0; k < n; k++) {
        for (int i = 0; i < j; i++) {
          S::Point2D const p = AffinComb(centre, j_coeff, domain.edge_point(k, (double)i/j));
          if (GetSurfacePoint(p, point)) {
            S::Point2D const bilinear =
              domain.MapNSidedDomainTo4SidedDomain(p, 0, SDomain::BilinearSweep);
            pg.push_back(BSS_Fitter<S::PointType>::ParamPointType(point, bilinear));
            if (!parametrization_check)
              fitted_points.push_back(point);
          }
        }
      }
    }
        
    fitter.DegreeU() = 3;
    fitter.DegreeV() = 3;

    fitter.KnotVectorU() = edge_curves[0].theKnotVector();
    fitter.KnotVectorV() = edge_curves[1].theKnotVector();

    size_t const nr_ctrl[] = { edge_curves[0].NrControlValues(),
                               edge_curves[1].NrControlValues() };

    // Setup control points constraints
    BSS_Fitter<S::PointType>::ConstraintsType cst(nr_ctrl);
    for (size_t i = 0; i < 2; ++i) {
      size_t end = nr_ctrl[1-i] - 1;
      S::CurveType::ControlVectorType const &ctrl1 = edge_curves[i].theControlVector();
      S::CurveType::ControlVectorType const &ctrl2 = edge_curves[i+2].theControlVector();
      for (size_t j = 0; j < nr_ctrl[i]; ++j) {
        if (i == 0) {
          cst[j][ 0 ] = ctrl1[j];
          cst[j][end] = ctrl2[j];
        } else {
          cst[end][j] = ctrl1[j];
          cst[ 0 ][j] = ctrl2[j];
        }
      }
    }
    std::swap(fitter.CtrlConstraints(), cst);

    fitter.Enlargement() = 0.0;

    fitter.AddFunctional(new BSSF_FN_LsqDistance<S::PointType>());
    fitter.AddFunctional(new BSSF_FN_Curvature<S::PointType>());
    fitter.OptimizeParameters() = false;

    QString dump_file = getDumpFilePath(sf->GetFilename());
    if (!dump_file.isEmpty())
      fitter.DumpOut( dump_file.toAscii().data());

    if (progress)
      SketchesWindow::sketches_window->GetProgressIndicator()->setProgress(++(*progress));
  } else {
    S::Point2D const &center = domain.Center();
    int u_step = OptionHandler::instance().GetOption("Patch Sample Density")->IntValue();
    // 0. Generate dividing curves (from the side to the center)
    std::vector<S::CurveType> dividing_curves; dividing_curves.reserve(n);
    std::vector<double> edge_midpoints; edge_midpoints.reserve(n); edge_midpoints.resize(n);
    std::vector<S::Point2D> edge_mid_u; edge_mid_u.reserve(n); edge_mid_u.resize(n);
    for (size_t i = 0; i < n; ++i) {
      BSC_Fitter<S::PointType> bf;
      bf.Degree() = 3;
      BSC_Fitter<S::PointType>::ParamPointGroupType ppoints(tol);

      double edge_half;
      {
        SRibbon *r = Ribbons()[i];
        S::CurveType c = *r->FreeCurve()->Curve();
        double mid = snapToKnot((r->Low() + r->Up()) / 2.0, c);
        edge_midpoints[i] = mid;
        edge_half = (mid - r->Low()) / (r->Up() - r->Low());
      }
      S::Point2D const u0 = domain.edge_point(i, edge_half);
      edge_mid_u[i] = u0;
      for (int j = 0; j <= u_step; j++) {
        double const u = (double)j / u_step;
        S::Point2D const param = AffinComb(u0, u, center);
        S::PointType point;
        if (GetSurfacePoint(param, point))
          ppoints.push_back(BSC_Fitter<S::PointType>::ParamPointType(point, u));
      }
      BSC_Fitter<S::PointType>::ConstraintsType cst(ppoints.size());
      cst.front() = ppoints.front().Point();
      cst.back() = ppoints.back().Point();
      std::swap(bf.CtrlConstraints(), cst);
      bf.ParamPointSet().push_back(ppoints);
      bf.AddFunctional(new BSCF_FN_LsqDistance<S::PointType>());
      bf.AddFunctional(new BSCF_FN_Curvature<S::PointType>());
      bf.Enlargement() = 0.0;
      bf.OptimizeParameters() = false;
      bf.Fit();
      dividing_curves.push_back(bf.Curve());
      // TODO: this shouldn't be needed - why does it not interpolate the constraints?
      size_t end = dividing_curves.back().NrControlValues() - 1;
      dividing_curves.back().theControlVector()[0] = ppoints.front().Point();
      dividing_curves.back().theControlVector()[end] = ppoints.back().Point();
    }
    for (size_t i = 0; i < n; ++i) {
      // 1. Generate edge curves
      S::CurveType edge_curves[2];
      for (size_t j = 0; j < 2; ++j) {
        size_t ip = (i + j) % n;
        SRibbon *r = Ribbons()[ip];
        edge_curves[j] = *r->FreeCurve()->Curve();
        double mid = edge_midpoints[ip];
        if (j == 0) {
          double end = snapToKnot(r->Up(), edge_curves[j]);
          if (mid > end) {
            Real const lower = edge_curves[j].LowerParam();
            Real const upper = edge_curves[j].UpperParam();
            edge_curves[j].Reparametrize((end-lower)+mid, (mid-upper)+end);
            std::swap(mid, end);
            mid = snapToKnot(mid, edge_curves[j]);
            end = snapToKnot(end, edge_curves[j]);
          }
          edge_curves[j].Trim(mid, end);
        } else {
          double start = snapToKnot(r->Low(), edge_curves[j]);
          if (start > mid) {
            Real const lower = edge_curves[j].LowerParam();
            Real const upper = edge_curves[j].UpperParam();
            edge_curves[j].Reparametrize((mid-lower)+start, (start-upper)+mid);
            std::swap(mid, start);
            mid = snapToKnot(mid, edge_curves[j]);
            start = snapToKnot(start, edge_curves[j]);
          }
          edge_curves[j].Trim(start, mid);
        }
        edge_curves[j].Reparametrize(0.0, 1.0);
      }

      // 2. Unify the opposite knot vectors
      S::CurveType div_curves[2];
      div_curves[0] = dividing_curves[(i+1)%n];
      div_curves[0].Reverse();
      div_curves[0].Reparametrize(0.0, 1.0);
      div_curves[1] = dividing_curves[i];
      div_curves[1].Reparametrize(0.0, 1.0);
      unifyKnots(edge_curves[0], div_curves[0]);
      unifyKnots(edge_curves[1], div_curves[1]);

      // 3. Fit
      SurfaceFitter::FitterType& fitter = sf->AddFitter(object_id, i);

      BSS_Fitter<S::PointType>::ParamPointSetReference point_set = fitter.ParamPointSet();
      point_set.push_back( BSS_Fitter<S::PointType>::ParamPointGroupType(tol));
      BSS_Fitter<S::PointType>::ParamPointGroupReference pg(point_set.front());

      S::Point2D const u1 = edge_mid_u[i];
      S::Point2D const v1 = edge_mid_u[(i+1)%n];
      S::Point2D const u_end = domain.edge_point(i, 1.0);

      int v_step = u_step;

      for (int j = 0; j <= u_step; j++) {
        double const u = (double)j / u_step;
        S::Point2D const u_v0 = AffinComb(u1, u, u_end);
        S::Point2D const u_v1 = AffinComb(center, u, v1);
        for (int k = 0; k <= v_step; k++) {
          double const v = (double)k / v_step;
          S::Point2D const param(AffinComb(u_v0, v, u_v1));
          S::PointType point;
          if (GetSurfacePoint(param, point)) {
            pg.push_back(BSS_Fitter<S::PointType>::ParamPointType(point, S::Point2D(u, v)));
            if (!parametrization_check)
              fitted_points.push_back(point);
          }
        }
      }

      fitter.DegreeU() = 3;
      fitter.DegreeV() = 3;

      fitter.KnotVectorU() = edge_curves[0].theKnotVector();
      fitter.KnotVectorV() = edge_curves[1].theKnotVector();

      size_t const nr_ctrl[] = { edge_curves[0].NrControlValues(),
                                 edge_curves[1].NrControlValues() };
#if 1
      size_t const nr_div_ctrl[] = { div_curves[0].NrControlValues(),
                                     div_curves[1].NrControlValues() };
      if (nr_ctrl[0] != nr_div_ctrl[0] || nr_ctrl[1] != nr_div_ctrl[1])
        {
          S_ERROR(QString("QuadFit: different number of ctrlpoints! edge_crv(%1,%2) div_crv(%3,%4)")
                  .arg(nr_ctrl[0]).arg(nr_ctrl[1]).arg(nr_div_ctrl[0]).arg(nr_div_ctrl[1]));

          QString knots_str;
          for (int i1=0; i1<2; i1++)
            {
              for (int i2=0; i2<2; i2++)
                {
                  knots_str = QString("QuadFit: knots on %1[%2]: ").arg(i1 ? "div_curves" : "edge_curves").arg(i2);
                  S::CurveType* crv = (i1 ? &div_curves[i2] : &edge_curves[i2]);
                  for (UInt i3=0, i3e=crv->theKnotVector().size(); i3<i3e; i3++)
                    {
                      knots_str += QString("%1%2").arg(i3 ? "," : "").arg(crv->theKnotVector()[i3]);
                    }
                  S_ERROR(knots_str);
                }
            }
        }
#endif
      // Setup control points constraints
      BSS_Fitter<S::PointType>::ConstraintsType cst(nr_ctrl);
      for (size_t k = 0; k < 2; ++k) {
        size_t end = nr_ctrl[1-k] - 1;
        S::CurveType::ControlVectorType const &ctrl1 = edge_curves[k].theControlVector();
        S::CurveType::ControlVectorType const &ctrl2 = div_curves[k].theControlVector();
        for (size_t j = 0; j < nr_ctrl[k]; ++j) {
          if (k == 0) {
            cst[j][ 0 ] = ctrl1[j];
            cst[j][end] = ctrl2[j];
          } else {
            cst[end][j] = ctrl1[j];
            cst[ 0 ][j] = ctrl2[j];
          }
        }
      }
      std::swap(fitter.CtrlConstraints(), cst);

      fitter.Enlargement() = 0.0;

      fitter.AddFunctional(new BSSF_FN_LsqDistance<S::PointType>());
      fitter.AddFunctional(new BSSF_FN_Curvature<S::PointType>());
      fitter.OptimizeParameters() = false;

      QString dump_file = getDumpFilePath(sf->GetFilename(), i);
      if (!dump_file.isEmpty())
        fitter.DumpOut( dump_file.toAscii().data());

      if (progress)
        SketchesWindow::sketches_window->GetProgressIndicator()->setProgress(++(*progress));
    }
  }
}

#endif  // NO_SURFACE_FIT

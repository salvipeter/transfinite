#include <cstring>

#include "transfinite-kapi.hh"

#include "surface-side-based.hh"
#include "surface-corner-based.hh"
#include "surface-generalized-coons.hh"
#include "surface-composite-ribbon.hh"

using Transfinite::DoubleVector;
using Transfinite::Point3D;
using Transfinite::PointVector;
using Transfinite::CurveVector;
using Transfinite::BSCurve;
using Transfinite::BSSurface;

struct BSplineCurve {
  BSplineCurve(const std::shared_ptr<BSCurve> &curve) : ptr_(curve) { }
  std::shared_ptr<BSCurve> ptr_;
};

struct BSplineSurface {
  BSplineSurface(const BSSurface &surface) : obj_(surface) { }
  BSSurface obj_;
};

long int transfinite_surface(BSplineCurve **curves, int curves_length,
                             char *surface_type,
                             char *split_or_trim,
                             double fit_tolerance,
                             BSplineSurface ***surfaces, int *surfaces_length) {
  CurveVector cv; cv.reserve(curves_length);
  for (int i = 0; i < curves_length; ++i)
    cv.push_back(curves[i]->ptr_);

  std::shared_ptr<Transfinite::Surface> surf;
  if (std::strcmp(surface_type, "SB") == 0)
    surf = std::make_shared<Transfinite::SurfaceSideBased>();
  else if (std::strcmp(surface_type, "CB") == 0)
    surf = std::make_shared<Transfinite::SurfaceCornerBased>();
  else if (std::strcmp(surface_type, "GC") == 0)
    surf = std::make_shared<Transfinite::SurfaceGeneralizedCoons>();
  else if (std::strcmp(surface_type, "CR") == 0)
    surf = std::make_shared<Transfinite::SurfaceCompositeRibbon>();
  else {
    return -1;                  // invalid surface type
  }

  surf->setCurves(cv);
  surf->setupLoop();            // TODO: this should validate the boundary curves [error code -2]
  surf->update();

  if (std::strcmp(split_or_trim, "split") == 0) {
    std::vector<BSSurface> result = surf->fitCentralSplit(fit_tolerance);
    *surfaces_length = result.size();
    *surfaces = new BSplineSurface*[*surfaces_length];
    for (int i = 0; i < *surfaces_length; ++i)
      (*surfaces)[i] = new BSplineSurface(result[i]);
  } else if (std::strcmp(split_or_trim, "trim") == 0) {
    BSSurface result = surf->fitTrimmed(fit_tolerance);
    *surfaces_length = 1;
    *surfaces = new BSplineSurface*[1];
    (*surfaces)[0] = new BSplineSurface(result);
  } else {
    return -3;                  // invalid fit mode
  }

  return 0;
}

long int create_B_SPLINE_CURVE(int degree,
                               double *knots, int knots_length,
                               double **coeffs, int coeffs_length,
                               BSplineCurve **result) {
  DoubleVector knots_v(knots, knots + knots_length);
  PointVector coeffs_v; coeffs_v.reserve(coeffs_length);
  for (int i = 0; i < coeffs_length; ++i)
    coeffs_v.push_back(Point3D(coeffs[i][0], coeffs[i][1], coeffs[i][2]));
  *result = new BSplineCurve(std::make_shared<BSCurve>(degree, knots_v, coeffs_v));
  return 0;
}

long int degree_B_SPLINE_CURVE(BSplineCurve *curve,
                               int *degree) {
  *degree = curve->ptr_->degree();
  return 0;
}

long int knots_B_SPLINE_CURVE(BSplineCurve *curve,
                              double **knots, int *knots_length) {
  DoubleVector kv = curve->ptr_->knotVector();
  *knots_length = kv.size();
  *knots = new double[*knots_length];
  for (int i = 0; i < *knots_length; ++i)
    (*knots)[i] = kv[i];
  return 0;
}

long int coeffs_B_SPLINE_CURVE(BSplineCurve *curve,
                               double ***coeffs, int *coeffs_length) {
  *coeffs_length = curve->ptr_->nrControlPoints();
  *coeffs = new double*[*coeffs_length];
  for (int i = 0; i < *coeffs_length; ++i) {
    const Point3D &p = curve->ptr_->controlPoint(i);
    (*coeffs)[i] = new double[3];
    (*coeffs)[i][0] = p[0];
    (*coeffs)[i][1] = p[1];
    (*coeffs)[i][2] = p[2];
  }
  return 0;
}

long int delete_B_SPLINE_CURVE(BSplineCurve *curve) {
  delete curve;
  return 0;
}

long int degrees_B_SPLINE_SURFACE(BSplineSurface *surface,
                                  int *deg_u, int *deg_v) {
  *deg_u = surface->obj_.deg_u_;
  *deg_v = surface->obj_.deg_v_;
  return 0;
}

long int knots_B_SPLINE_SURFACE(BSplineSurface *surface,
                                double **knots_u, int *knots_u_length,
                                double **knots_v, int *knots_v_length) {
  *knots_u_length = surface->obj_.knots_u_.size();
  *knots_u = new double[*knots_u_length];
  for (int i = 0; i < *knots_u_length; ++i)
    (*knots_u)[i] = surface->obj_.knots_u_[i];

  *knots_v_length = surface->obj_.knots_v_.size();
  *knots_v = new double[*knots_v_length];
  for (int i = 0; i < *knots_v_length; ++i)
    (*knots_v)[i] = surface->obj_.knots_v_[i];

  return 0;
}

long int n_coeffs_B_SPLINE_SURFACE(BSplineSurface *surface,
                                   int *n_coeffs_u, int *n_coeffs_v) {
  *n_coeffs_u = surface->obj_.cnet_.size();
  *n_coeffs_v = surface->obj_.cnet_[0].size();
  return 0;
}

long int coeffs_B_SPLINE_SURFACE(BSplineSurface *surface,
                                 double ***coeffs, int *coeffs_length) {
  size_t n = surface->obj_.cnet_.size(), m = surface->obj_.cnet_[0].size();
  *coeffs_length = n * m;
  *coeffs = new double*[*coeffs_length];
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < m; ++j) { // TODO: is this the good order for GeoDuck?
      size_t k = i * m + j;
      const Point3D &p = surface->obj_.cnet_[i][j];
      (*coeffs)[k] = new double[3];
      (*coeffs)[k][0] = p[0];
      (*coeffs)[k][1] = p[1];
      (*coeffs)[k][2] = p[2];
    }
  return 0;
}

long int trim_curves_B_SPLINE_SURFACE(BSplineSurface *surface,
                                      BSplineCurve ***curves, int *curves_length) {
  const CurveVector &cv = surface->obj_.param_curves_;
  *curves_length = cv.size();
  *curves = new BSplineCurve*[*curves_length];
  for (int i = 0; i < *curves_length; ++i)
    (*curves)[i] = new BSplineCurve(cv[i]);
  return 0;
}

long int delete_B_SPLINE_SURFACE(BSplineSurface *surface) {
  delete surface;
  return 0;
}

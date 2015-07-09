#pragma once

#include <ctime>

#include "geometry.hh"

#include "IgesExport.hh"

class IGES::IGESImpl {
  using NURBSCurv = BSplineCurve<Point<3, double>>;
  using NURBSSurf = BSplineSurface<Point<3, double>>;
public:
  IGESImpl(std::string filename)
    : filter_(IgesGlobal::create("Sketches", "Sketches", "1.0.2",
                                 "Product", "Product", "Author", "ShapEx",
                                 "protocol", std::time(nullptr), 0.1, filename.c_str(),
                                 IgesGlobal::IGU_MILLIMETERS, IgesGlobal::IGDSTD_NONE, 1.0, 5.0),
              IgesExport::IGVERSION_5_2) {}
  void writeSurface(const BSSurface &surface) {
    NURBSSurf s = convertToNURBS(surface);
    if (surface.curves_.empty()) {
      filter_ << s;
      return;
    }
    std::vector<NURBSCurv> cs;
    std::transform(surface.curves_.begin(), surface.curves_.end(), std::back_inserter(cs),
                   [](const std::shared_ptr<BSCurve> &c) { return convertToNURBS(*c); });

    TrimmedBSplineSurface tsurf;
    tsurf.surface = &s;
    // tsurf.domain_loop_sets.push_back(TrimmedBSplineSurface::DC_LoopType());
    // for (size_t i = 0; i < size; i++)
    //   tsurf.domain_loop_sets.back().push_back(&domain_trim_curves[i]);
    tsurf.space_loop_sets.push_back(TrimmedBSplineSurface::C_LoopType());
    for (const auto &c : cs)
      tsurf.space_loop_sets.back().push_back(&c);
    filter_ << tsurf;
  }
  void close() { filter_.close(); }
private:
  static NURBSCurv convertToNURBS(const BSCurve &curve) {
    NURBSCurv::KnotVectorType knots;
    DoubleVector knots_orig = curve.knotVector();
    knots.assign(knots_orig.begin(), knots_orig.end());
    size_t n = curve.nrControlPoints();
    NURBSCurv::ControlVectorType cv(n);
    for (size_t i = 0; i < n; ++i) {
      const Point3D &p = curve.controlPoint(i);
      cv[i] = Point<3, double>(p[0], p[1], p[2]);
    }
    return NURBSCurv(curve.degree(), knots, cv);
  }
  static NURBSSurf convertToNURBS(const BSSurface &surface) {
    NURBSSurf::KnotVectorType knots_u, knots_v;
    knots_u.assign(surface.knots_u_.begin(), surface.knots_u_.end());
    knots_v.assign(surface.knots_v_.begin(), surface.knots_v_.end());
    size_t n[] = { surface.cnet_.size(), surface.cnet_[0].size() };
    NURBSSurf::ControlNetType cnet(n);
    for (size_t i = 0; i < n[0]; ++i)
      for (size_t j = 0; j < n[1]; ++j) {
        const Point3D &p = surface.cnet_[i][j];
        cnet[i][j] = Point<3, double>(p[0], p[1], p[2]);
      }
    return NURBSSurf(surface.deg_u_, surface.deg_v_, knots_u, knots_v, cnet);
  }

  IgesExport filter_;
};

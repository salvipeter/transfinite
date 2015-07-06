#pragma once

#include <ctime>

#include "geometry.hh"

#include "IgesExport.hh"

class IGES::IGESImpl {
public:
  IGESImpl(std::string filename)
    : filter_(IgesGlobal::create("Sketches", "Sketches", "1.0.2",
                                 "Product", "Product", "Author", "ShapEx",
                                 "protocol", std::time(nullptr), 0.1, filename.c_str(),
                                 IgesGlobal::IGU_MILLIMETERS, IgesGlobal::IGDSTD_NONE, 1.0, 5.0),
              IgesExport::IGVERSION_5_2) {}
  void writeSurface(const BSSurface &surface) {
    BSplineSurface<Point<3, double>>::KnotVectorType knots_u, knots_v;
    knots_u.assign(surface.knots_u_.begin(), surface.knots_u_.end());
    knots_v.assign(surface.knots_v_.begin(), surface.knots_v_.end());
    size_t n[] = { surface.cnet_.size(), surface.cnet_[0].size() };
    BSplineSurface<Point<3, double>>::ControlNetType cnet(n);
    for (size_t i = 0; i < n[0]; ++i)
      for (size_t j = 0; j < n[1]; ++j) {
        const Point3D &p = surface.cnet_[i][j];
        cnet[i][j] = Point<3, double>(p[0], p[1], p[2]);
      }
    BSplineSurface<Point<3, double>> s(surface.deg_u_, surface.deg_v_,
                                       knots_u, knots_v, cnet);
    filter_ << s;
  }
  void writeTrimmedSurface(const BSSurface &surface, const CurveVector &curves) {
    // TODO
  }
  void close() { filter_.close(); }
private:
  IgesExport filter_;
};

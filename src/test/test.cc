#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "domain.hh"
#include "ribbon.hh"
#include "surface-side-based.hh"
#include "surface-corner-based.hh"
#include "surface-generalized-coons.hh"
#include "surface-composite-ribbon.hh"

CurveVector readLOP(std::string filename) {
  std::ifstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return CurveVector();
  }

  size_t n, deg, nk, nc;
  DoubleVector knots;
  PointVector cpts;
  CurveVector result;

  f >> n;
  result.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    f >> deg;
    f >> nk;
    knots.resize(nk);
    for (size_t j = 0; j < nk; ++j)
      f >> knots[j];
    f >> nc;
    cpts.resize(nc);
    for (size_t j = 0; j < nc; ++j)
      f >> cpts[j][0] >> cpts[j][1] >> cpts[j][2];
    result.push_back(std::make_shared<BSCurve>(deg, knots, cpts));
  }

  return result;
}

void surfaceTest(std::string filename, std::string type, size_t resolution,
                 std::shared_ptr<Transfinite::Surface> &&surf) {
  CurveVector cv = readLOP("../../models/" + filename + ".lop");

  surf->setCurves(cv);
  surf->setupLoop();             // should be called after curve pointers changed
  surf->update();                // should be called after curves changed
  surf->eval(resolution).writeOBJ("../../models/" + filename + "-" + type + ".obj");

  std::shared_ptr<const Transfinite::Domain> domain = surf->domain();
  const Point2DVector &v = domain->vertices();

  size_t n = cv.size();
  size_t res = 40;
  double step = 1.0e-4;

  // Positional test
  double max_pos_error = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= res; ++j) {
      double s = (double)j / (double)res;
      Point2D uv = domain->edgePoint(i, s);
      Point3D p = surf->eval(uv), q = cv[i]->eval(s);
      max_pos_error = std::max(max_pos_error, (p - q).norm());
    }
  }

  // Tangential test
  double max_tan_error = 0.0;
  VectorVector der;
  for (size_t i = 0; i < n; ++i) {
    Vector2D perp = v[i] - v[(i+n-1)%n];
    perp = Vector2D(perp[1], -perp[0]);
    if ((domain->center() - v[i]) * perp < 0)
      perp = -perp;
    perp.normalize();
    for (size_t j = 0; j <= res; ++j) {
      double s = (double)j / (double)res;
      Point2D uv = domain->edgePoint(i, s);
      Point2D uv2 = uv + perp * step;
      Point3D p = surf->eval(uv), q = surf->eval(uv2);
      cv[i]->eval(s, 1, der);
      Vector3D surf_normal = (der[1] ^ (q - p)).normalize();
      Vector3D normal = surf->ribbon(i)->normal(s);
      double angle = acos(surf_normal * normal);
      max_tan_error = std::max(max_tan_error, angle);
    }
  }

  std::cout << type << ":" << std::endl;
  std::cout << "  positional error: " << max_pos_error << std::endl;
  std::cout << "  tangential error: " << max_tan_error * 180.0 / M_PI << std::endl;

#ifndef NO_SURFACE_FIT

  std::vector<BSSurface> surfaces = surf->fitCentralSplit(0.01);
  IGES iges("../../models/" + filename + "-" + type + "-split.igs");
  for (const auto &s : surfaces)
    iges.writeSurface(s);
  iges.close();

  BSSurface trimmed = surf->fitTrimmed(0.01);
  IGES iges2("../../models/" + filename + "-" + type + "-trim.igs");
  iges2.writeSurface(trimmed);
  iges2.close();

#endif  // NO_SURFACE_FIT
}

int main(int argc, char **argv) {
#ifdef DEBUG
  std::cout << "Compiled in DEBUG mode" << std::endl;
#else  // !DEBUG
  std::cout << "Compiled in RELEASE mode" << std::endl;
#endif  // DEBUG

  size_t res = 15;
  std::string filename;
  if (argc < 2 || argc > 3) {
    std::cerr << "Usage: " << argv[0] << " model-name [resolution]" << std::endl;
    return 1;
  }
  filename = argv[1];
  if (argc == 3)
    res = atoi(argv[2]);

  surfaceTest(filename, "sb", res, std::make_shared<Transfinite::SurfaceSideBased>());
  surfaceTest(filename, "cb", res, std::make_shared<Transfinite::SurfaceCornerBased>());
  surfaceTest(filename, "gc", res, std::make_shared<Transfinite::SurfaceGeneralizedCoons>());
  surfaceTest(filename, "cr", res, std::make_shared<Transfinite::SurfaceCompositeRibbon>());

  return 0;
}

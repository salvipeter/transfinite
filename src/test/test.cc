#include <algorithm>
#include <iostream>

#include "domain.hh"
#include "ribbon.hh"
#include "surface-side-based.hh"
#include "surface-corner-based.hh"
#include "surface-generalized-bezier.hh"
#include "surface-generalized-coons.hh"
#include "surface-composite-ribbon.hh"
#include "surface-midpoint.hh"

#include "gb-fit.hh"
#include "io.hh"

void ribbonTest(std::string filename, size_t resolution, 
		double scaling, double ribbon_length) {
  CurveVector cv = readLOP("../../models/" + filename + ".lop");
  SurfaceSideBased surf;
  surf.setCurves(cv);
  surf.setupLoop();
  surf.update();

  // Fence
  TriMesh fence_mesh;
  size_t n = cv.size();
  size_t size = n * resolution * 2;
  PointVector pv; pv.reserve(size);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < resolution; ++j) {
      double u = (double)j / resolution;
      Point3D p = surf.ribbon(i)->curve()->eval(u);
      pv.push_back(p);
      p += surf.ribbon(i)->normal(u) * scaling;
      pv.push_back(p);
    }
  }
  fence_mesh.setPoints(pv);
  size_t index = 0;
  while (index < size - 2) {
    fence_mesh.addTriangle(index, index+1, index+2);
    ++index;
    fence_mesh.addTriangle(index, index+2, index+1);
    ++index;
  }
  fence_mesh.addTriangle(index, index+1, 0);
  fence_mesh.addTriangle(index+1, 1, 0);
  fence_mesh.writeOBJ("../../models/" + filename + "-fence.obj");

  // Ribbons
  TriMesh ribbon_mesh;
  pv.clear(); pv.reserve(n * (resolution + 1) * 2);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= resolution; ++j) {
      double u = (double)j / resolution;
      pv.push_back(surf.ribbon(i)->curve()->eval(u));
      pv.push_back(surf.ribbon(i)->eval(Point2D(u, ribbon_length)));
    }
  }
  ribbon_mesh.setPoints(pv);
  index = 0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < resolution; ++j) {
      ribbon_mesh.addTriangle(index, index+1, index+2);
      ++index;
      ribbon_mesh.addTriangle(index, index+2, index+1);
      ++index;
    }
    index += 2;
  }
  ribbon_mesh.writeOBJ("../../models/" + filename + "-ribbons.obj");  
}

void surfaceTest(std::string filename, std::string type, size_t resolution,
                 std::shared_ptr<Surface> &&surf) {
  CurveVector cv = readLOP("../../models/" + filename + ".lop");

  surf->setCurves(cv);
  surf->setupLoop();             // should be called after curve pointers changed
  surf->update();                // should be called after curves changed
  surf->eval(resolution).writeOBJ("../../models/" + filename + "-" + type + ".obj");

  std::shared_ptr<const Domain> domain = surf->domain();
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

  if (type == "mp") {
    Point3D midpoint(0,0,0);
    for (size_t i = 0; i < n; ++i)
      midpoint += surf->ribbon(i)->eval(Point2D(0.6, 0.4));
    midpoint /= n;
    dynamic_cast<SurfaceMidpoint *>(surf.get())->setMidpoint(midpoint);
    Point3D p = surf->eval(domain->center());
    std::cout << "  midpoint error: " << (p - midpoint).normSqr() << std::endl;
  }

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

void bezierTest(const std::string &filename) {
  SurfaceGeneralizedBezier surf = loadBezier("../../models/" + filename + ".gbp");

  // Generate mesh output
  TriMesh mesh = surf.eval(15);
  mesh.writeOBJ("../../models/bezier.obj");
  writeBezierControlPoints(surf, "../../models/bezier-cpts.obj");

  // Normal degree elevation
  SurfaceGeneralizedBezier sextic = elevateDegree(surf), sextic2, sextic3;
  writeBezierControlPoints(sextic, "../../models/bezier-elevated-cpts.obj");
  sextic.eval(15).writeOBJ("../../models/bezier-elevated.obj");
  SurfaceGeneralizedBezier elevated = surf;
  for (size_t i = 0; i < 30; ++i)
    elevated = elevateDegree(elevated);
  writeBezierControlPoints(elevated, "../../models/bezier-elevated-30-times.obj");
  for (size_t i = 0; i < 30; ++i)
    elevated = elevateDegree(elevated);
  writeBezierControlPoints(elevated, "../../models/bezier-elevated-60-times.obj");

  // Fit a sextic surface on the quintic mesh
  sextic2 = fitWithOriginal(sextic, mesh.points(), surf.domain()->parameters(15));
  sextic2.eval(15).writeOBJ("../../models/bezier-sextic.obj");
  writeBezierControlPoints(sextic2, "../../models/bezier-sextic-cpts.obj");

  // Fit again, now with projected parameterization
  sextic3 = sextic;
  for (size_t i = 0; i < 1; ++i)
    sextic3 = fitWithOriginal(sextic3, mesh.points(), parameterizePoints(sextic3, mesh.points()));
  sextic3.eval(15).writeOBJ("../../models/bezier-sextic-projected.obj");
  writeBezierControlPoints(sextic3, "../../models/bezier-sextic-projected-cpts.obj");

  // Fit again, now with smoothing
  sextic3 = sextic;
  for (size_t i = 0; i < 1; ++i)
    sextic3 = fitWithOriginal(sextic3, mesh.points(), parameterizePoints(sextic3, mesh.points()),
                              0.2);
  sextic3.eval(15).writeOBJ("../../models/bezier-sextic-projected-smooth.obj");
  writeBezierControlPoints(sextic3, "../../models/bezier-sextic-projected-smooth-cpts.obj");
}

void meshFitTest(const std::string &surfname, const std::string &meshname) {
  SurfaceGeneralizedBezier surf = loadBezier("../../models/" + surfname + ".gbp");
  TriMesh mesh = readOBJ("../../models/" + meshname + ".obj");

  surf.eval(15).writeOBJ("../../models/bezier.obj");
  writeBezierControlPoints(surf, "../../models/bezier-cpts.obj");

  SurfaceGeneralizedBezier fitted;
  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()));
  fitted.eval(15).writeOBJ("../../models/bezier-fitted.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-cpts.obj");
  saveBezier(fitted, "../../models/bezier-fitted.gbp");

  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()),
                           0.2);
  fitted.eval(15).writeOBJ("../../models/bezier-fitted-smooth02.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-smooth02-cpts.obj");
  saveBezier(fitted, "../../models/bezier-smooth02.gbp");

  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()),
                           0.5);
  fitted.eval(15).writeOBJ("../../models/bezier-fitted-smooth05.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-smooth05-cpts.obj");
  saveBezier(fitted, "../../models/bezier-smooth05.gbp");

  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()),
                           1.0);
  fitted.eval(15).writeOBJ("../../models/bezier-fitted-smooth10.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-smooth10-cpts.obj");
  saveBezier(fitted, "../../models/bezier-smooth10.gbp");
}

void classATest() {
#ifdef NO_SURFACE_FIT
  CurveVector cv = readLOP("../../models/pocket6sided.lop");

  CurveVector normal, class_a;
  for (const auto &c : cv) {
    VectorVector der;
    Point3D pa = c->eval(0.0, 1, der);
    Vector3D va = der[1];
    Point3D pb = c->eval(1.0, 1, der);
    Vector3D vb = der[1];
    PointVector pv = {pa, pa + va / 3.0, pb - vb / 3.0, pb};
    normal.push_back(std::make_shared<BSCurve>(pv));
    BCurve bc; bc.fitClassA(14, pa, va, pb, vb);
    class_a.push_back(std::make_shared<BSCurve>(bc.controlPoints()));
  }

  SurfaceCornerBased surf;
  surf.setCurves(normal);
  surf.setupLoop();
  surf.update();
  surf.eval(30).writeOBJ("../../models/bezier-normal.obj");
  surf.setCurves(class_a);
  surf.setupLoop();
  surf.update();
  surf.eval(30).writeOBJ("../../models/bezier-class-a.obj");
#endif  // NO_SURFACE_FIT
}

int main(int argc, char **argv) {
#ifdef DEBUG
  std::cout << "Compiled in DEBUG mode" << std::endl;
#else  // !DEBUG
  std::cout << "Compiled in RELEASE mode" << std::endl;
#endif  // DEBUG

  size_t res = 15;
  double scaling = 20.0;
  double ribbon_length = 0.25;

  if (argc < 2 || argc > 5) {
    std::cerr << "Usage:\n"
              << argv[0] << " model-name [resolution] [fence-scaling] [ribbon-length]" << std::endl
              << argv[0] << " bezier [model-name]" << std::endl
              << argv[0] << " mesh-fit [model-name] [mesh-name]" << std::endl;
    return 1;
  }

  std::string filename = argv[1];

  if (filename == "bezier") {
    if (argc == 2)
      bezierTest("cagd86");
    else
      bezierTest(argv[2]);
    return 0;
  } else if (filename == "class-a") {
    classATest();
    return 0;
  } else if (filename == "mesh-fit") {
    if (argc < 4) {
      std::cerr << "Not enough parameters!" << std::endl;
      return 1;
    }
    meshFitTest(argv[2], argv[3]);
    return 0;
  }

  if (argc >= 3)
    res = atoi(argv[2]);
  if (argc >= 4)
    scaling = strtod(argv[3], nullptr);
  if (argc >= 5)
    ribbon_length = strtod(argv[4], nullptr);

  ribbonTest(filename, res, scaling, ribbon_length);

  surfaceTest(filename, "sb", res, std::make_shared<SurfaceSideBased>());
  surfaceTest(filename, "cb", res, std::make_shared<SurfaceCornerBased>());
  surfaceTest(filename, "gc", res, std::make_shared<SurfaceGeneralizedCoons>());
  surfaceTest(filename, "cr", res, std::make_shared<SurfaceCompositeRibbon>());
  surfaceTest(filename, "mp", res, std::make_shared<SurfaceMidpoint>());

  return 0;
}

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "domain.hh"
#include "ribbon.hh"
#include "surface-side-based.hh"
#include "surface-corner-based.hh"
#include "surface-generalized-bezier.hh"
#include "surface-generalized-coons.hh"
#include "surface-composite-ribbon.hh"

using namespace Transfinite;

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

void bezierTest() {
  SurfaceGeneralizedBezier surf;
  surf.initNetwork(5, 5);
  surf.setCentralControlPoint(Point3D(62.6431, 9.65198, 22.5239));

  // Auto-generated data, contains redundancies
  surf.setControlPoint(0,0,0,Point3D(61.7919, 58, 116.732));
  surf.setControlPoint(0,1,0,Point3D(92.9346, 55.6215, 91.3775)); 
  surf.setControlPoint(0,2,0,Point3D(109.87, 58.9153, 54.4686)); 
  surf.setControlPoint(0,3,0,Point3D(123.086, 58.4892, 17.2603)); 
  surf.setControlPoint(0,4,0,Point3D(130.802, 57.9993, -20.5265)); 
  surf.setControlPoint(0,5,0,Point3D(131.793, 58, -60)); 
  surf.setControlPoint(0,0,1,Point3D(51.523, 47.9294, 117.163)); 
  surf.setControlPoint(0,1,1,Point3D(83.5958, 44.6549, 93.0348)); 
  surf.setControlPoint(0,2,1,Point3D(98.5855, 44.1235, 55.9962)); 
  surf.setControlPoint(0,3,1,Point3D(108.744, 39.0439, 18.0016)); 
  surf.setControlPoint(0,4,1,Point3D(114.053, 35.1734, -20.5281)); 
  surf.setControlPoint(0,5,1,Point3D(115.054, 35.1673, -60.0014)); 
  surf.setControlPoint(0,0,2,Point3D(40.1278, 40.5432, 117.313)); 
  surf.setControlPoint(0,1,2,Point3D(62.8623, 33.4245, 91.669)); 
  surf.setControlPoint(0,2,2,Point3D(74.2388, 25.4256, 50.9729)); 
  surf.setControlPoint(0,3,2,Point3D(84.1649, 17.9494, 8.6288)); 
  surf.setControlPoint(0,4,2,Point3D(94.8741, 12.6111, -28.4806)); 
  surf.setControlPoint(0,5,2,Point3D(98.3238, 12.2747, -61.3907)); 
  surf.setControlPoint(1,0,0,Point3D(131.793, 58, -60)); 
  surf.setControlPoint(1,1,0,Point3D(115.054, 35.1673, -60.0014)); 
  surf.setControlPoint(1,2,0,Point3D(98.3238, 12.2747, -61.3907)); 
  surf.setControlPoint(1,3,0,Point3D(88.4707, -13.7322, -60.2359)); 
  surf.setControlPoint(1,4,0,Point3D(86.0242, -43.0892, -56.2415)); 
  surf.setControlPoint(1,5,0,Point3D(72.1941, -67.4299, -59.9941)); 
  surf.setControlPoint(1,0,1,Point3D(130.802, 57.9993, -20.5265)); 
  surf.setControlPoint(1,1,1,Point3D(114.053, 35.1734, -20.5281)); 
  surf.setControlPoint(1,2,1,Point3D(94.8741, 12.6111, -28.4806)); 
  surf.setControlPoint(1,3,1,Point3D(81.7361, -12.943, -36.0907)); 
  surf.setControlPoint(1,4,1,Point3D(76.7726, -41.9252, -38.7019)); 
  surf.setControlPoint(1,5,1,Point3D(62.7974, -66.1704, -42.5395)); 
  surf.setControlPoint(1,0,2,Point3D(123.086, 58.4892, 17.2603)); 
  surf.setControlPoint(1,1,2,Point3D(108.744, 39.0439, 18.0016)); 
  surf.setControlPoint(1,2,2,Point3D(84.1649, 17.9494, 8.6288)); 
  surf.setControlPoint(1,3,2,Point3D(67.6886, -9.86071, -5.51556)); 
  surf.setControlPoint(1,4,2,Point3D(60.8794, -40.719, -20.0232)); 
  surf.setControlPoint(1,5,2,Point3D(51.0532, -65.35, -29.3928)); 
  surf.setControlPoint(2,0,0,Point3D(72.1941, -67.4299, -59.9941)); 
  surf.setControlPoint(2,1,0,Point3D(62.7974, -66.1704, -42.5395)); 
  surf.setControlPoint(2,2,0,Point3D(51.0532, -65.35, -29.3928)); 
  surf.setControlPoint(2,3,0,Point3D(36.6922, -64.8978, -20.2125)); 
  surf.setControlPoint(2,4,0,Point3D(19.6749, -64.7858, -15.0237)); 
  surf.setControlPoint(2,5,0,Point3D(0, -64.5101, -14.7862)); 
  surf.setControlPoint(2,0,1,Point3D(86.0242, -43.0892, -56.2415)); 
  surf.setControlPoint(2,1,1,Point3D(76.7726, -41.9252, -38.7019)); 
  surf.setControlPoint(2,2,1,Point3D(60.8794, -40.719, -20.0232)); 
  surf.setControlPoint(2,3,1,Point3D(40.9138, -39.6976, -3.51752)); 
  surf.setControlPoint(2,4,1,Point3D(19.675, -39.1309, 7.13632)); 
  surf.setControlPoint(2,5,1,Point3D(3.84774e-05, -38.8132, 7.32499)); 
  surf.setControlPoint(2,0,2,Point3D(88.4707, -13.7322, -60.2359)); 
  surf.setControlPoint(2,1,2,Point3D(81.7361, -12.943, -36.0907)); 
  surf.setControlPoint(2,2,2,Point3D(67.6886, -9.86071, -5.51556)); 
  surf.setControlPoint(2,3,2,Point3D(44.4584, -6.05418, 14.8337)); 
  surf.setControlPoint(2,4,2,Point3D(18.1019, -2.51218, 18.9975)); 
  surf.setControlPoint(2,5,2,Point3D(-0.00130326, -2.27912, 19.0098)); 
  surf.setControlPoint(3,0,0,Point3D(0, -64.5101, -14.7862)); 
  surf.setControlPoint(3,1,0,Point3D(3.84774e-05, -38.8132, 7.32499)); 
  surf.setControlPoint(3,2,0,Point3D(-0.00130326, -2.27912, 19.0098)); 
  surf.setControlPoint(3,3,0,Point3D(-0.000653163, 13.7061, 50.0875)); 
  surf.setControlPoint(3,4,0,Point3D(0.0070725, 16.102, 87.221)); 
  surf.setControlPoint(3,5,0,Point3D(-0.00806045, 31.2519, 117.651)); 
  surf.setControlPoint(3,0,1,Point3D(19.6749, -64.7858, -15.0237)); 
  surf.setControlPoint(3,1,1,Point3D(19.675, -39.1309, 7.13632)); 
  surf.setControlPoint(3,2,1,Point3D(18.1019, -2.51218, 18.9975)); 
  surf.setControlPoint(3,3,1,Point3D(16.0069, 13.6112, 50.2838)); 
  surf.setControlPoint(3,4,1,Point3D(14.443, 16.1255, 87.563)); 
  surf.setControlPoint(3,5,1,Point3D(14.4281, 31.3009, 117.98)); 
  surf.setControlPoint(3,0,2,Point3D(36.6922, -64.8978, -20.2125)); 
  surf.setControlPoint(3,1,2,Point3D(40.9138, -39.6976, -3.51752)); 
  surf.setControlPoint(3,2,2,Point3D(44.4584, -6.05418, 14.8337)); 
  surf.setControlPoint(3,3,2,Point3D(43.6807, 15.3926, 47.466)); 
  surf.setControlPoint(3,4,2,Point3D(37.5815, 23.1714, 89.3225)); 
  surf.setControlPoint(3,5,2,Point3D(27.7633, 34.9634, 117.604)); 
  surf.setControlPoint(4,0,0,Point3D(-0.00806045, 31.2519, 117.651)); 
  surf.setControlPoint(4,1,0,Point3D(14.4281, 31.3009, 117.98)); 
  surf.setControlPoint(4,2,0,Point3D(27.7633, 34.9634, 117.604)); 
  surf.setControlPoint(4,3,0,Point3D(40.1278, 40.5432, 117.313)); 
  surf.setControlPoint(4,4,0,Point3D(51.523, 47.9294, 117.163)); 
  surf.setControlPoint(4,5,0,Point3D(61.7919, 58, 116.732)); 
  surf.setControlPoint(4,0,1,Point3D(0.0070725, 16.102, 87.221)); 
  surf.setControlPoint(4,1,1,Point3D(14.443, 16.1255, 87.563)); 
  surf.setControlPoint(4,2,1,Point3D(37.5815, 23.1714, 89.3225)); 
  surf.setControlPoint(4,3,1,Point3D(62.8623, 33.4245, 91.669)); 
  surf.setControlPoint(4,4,1,Point3D(83.5958, 44.6549, 93.0348)); 
  surf.setControlPoint(4,5,1,Point3D(92.9346, 55.6215, 91.3775)); 
  surf.setControlPoint(4,0,2,Point3D(-0.000653163, 13.7061, 50.0875)); 
  surf.setControlPoint(4,1,2,Point3D(16.0069, 13.6112, 50.2838)); 
  surf.setControlPoint(4,2,2,Point3D(43.6807, 15.3926, 47.466)); 
  surf.setControlPoint(4,3,2,Point3D(74.2388, 25.4256, 50.9729)); 
  surf.setControlPoint(4,4,2,Point3D(98.5855, 44.1235, 55.9962)); 
  surf.setControlPoint(4,5,2,Point3D(109.87, 58.9153, 54.4686)); 

  surf.setupLoop();

  // Generate mesh output
  surf.eval(15).writeOBJ("../../models/bezier.obj");
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
  std::string filename;
  if (argc < 2 || argc > 5) {
    std::cerr << "Usage: " << argv[0]
	      << " model-name [resolution] [fence-scaling] [ribbon-length]" << std::endl;
    return 1;
  }
  filename = argv[1];
  if (argc >= 3)
    res = atoi(argv[2]);
  if (argc >= 4)
    scaling = strtod(argv[3], nullptr);
  if (argc >= 5)
    ribbon_length = strtod(argv[4], nullptr);

  if (filename == "bezier") {
    bezierTest();
    return 0;
  } else if (filename == "class-a") {
    classATest();
    return 0;
  }

  ribbonTest(filename, res, scaling, ribbon_length);

  surfaceTest(filename, "sb", res, std::make_shared<SurfaceSideBased>());
  surfaceTest(filename, "cb", res, std::make_shared<SurfaceCornerBased>());
  surfaceTest(filename, "gc", res, std::make_shared<SurfaceGeneralizedCoons>());
  surfaceTest(filename, "cr", res, std::make_shared<SurfaceCompositeRibbon>());

  return 0;
}

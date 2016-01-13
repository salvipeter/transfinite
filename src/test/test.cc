#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "../Eigen/LU"

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

void writeBezierControlPoints(const SurfaceGeneralizedBezier &surf, const std::string &filename) {
  // Slow but simple implementation creating a nice mesh
  size_t n = surf.domain()->vertices().size();
  size_t d = surf.degree();
  size_t l = surf.layers();
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;

  auto findControlPoint = [n,d,cp](size_t i, size_t j, size_t k) -> size_t {
    for (size_t c = 1, side = 0, col = 0, row = 0; c < cp; ++c, ++col) {
      if (col >= d - row) {
        if (++side >= n) {
          side = 0;
          ++row;
        }
        col = row;
      }
      size_t side_m = (side + n - 1) % n, side_p = (side + 1) % n;
      if ((i == side && j == col && k == row) ||
          (i == side_m && j == d - row && k == col) ||
          (i == side_p && j == row && k == d - col))
        return c + 1;
    }
    return 0;
  };

  std::ofstream f(filename);
  Point3D p = surf.centralControlPoint();
  f << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    p = surf.controlPoint(side, col, row);
    f << "v " << p[0] << " " << p[1] << " " << p[2] << std::endl;
  }

  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j <= d / 2; ++j)
      for (size_t k = 0; k < l - 1; ++k) {
        size_t a = findControlPoint(i, j, k);
        size_t b = findControlPoint(i, j + 1, k);
        size_t c = findControlPoint(i, j + 1, k + 1);
        size_t d = findControlPoint(i, j, k + 1);
        f << "f " << a << " " << b << " " << c << " " << d << std::endl;
      }
  if (d % 2 == 0)
    for (size_t i = 0; i < n; ++i) {
      size_t im = (i + n - 1) % n;
      size_t a = findControlPoint(i, l - 1, l - 1);
      size_t b = findControlPoint(i, l, l - 1);
      size_t c = 1;
      size_t d = findControlPoint(im, l, l - 1);
      f << "f " << a << " " << b << " " << c << " " << d << std::endl;
    }
  else
    for (size_t i = 0; i < n; ++i) {
      size_t a = findControlPoint(i, l - 1, l - 1);
      size_t b = findControlPoint(i, l, l - 1);
      size_t c = 1;
      f << "f " << a << " " << b << " " << c << std::endl;
    }

  f.close();
}

SurfaceGeneralizedBezier elevateDegree(const SurfaceGeneralizedBezier &surf) {
  size_t n = surf.domain()->vertices().size();
  size_t d = surf.degree() + 1;
  SurfaceGeneralizedBezier result;
  result.initNetwork(n, d);
  size_t l = result.layers();

  // Set the corner control points
  for (size_t i = 0; i < n; ++i)
    result.setControlPoint(i, 0, 0, surf.controlPoint(i, 0, 0));

  // Set the boundary control points
  for (size_t j = 1; j < d; ++j) {
    double eta = (double)j / d;
    for (size_t i = 0; i < n; ++i)
      result.setControlPoint(i, j, 0,
                             surf.controlPoint(i, j - 1, 0) * eta +
                             surf.controlPoint(i, j, 0) * (1 - eta));
  }

  // Set all other control points (except the center)
  // When d is odd, the topmost layer of the original patch
  // is simulated by the control points of the adjacent boundary
  for (size_t j = 1; j <= d / 2; ++j) {
    double eta = (double)j / d;
    for (size_t k = 1; k < l; ++k) {
      double theta = (double)k / d;
      for (size_t i = 0; i < n; ++i) {
        Point3D p1, p2, p3, p4;
        p1 = surf.controlPoint(i, j - 1, k - 1);
        p2 = surf.controlPoint(i, j, k - 1);
        if (d % 2 == 0 || k < d / 2) {
          p3 = surf.controlPoint(i, j - 1, k);
          p4 = surf.controlPoint(i, j, k);
        } else {
          size_t im = (i + n - 1) % n;
          p3 = surf.controlPoint(im, d - k - 1, j - 1);
          if (j == d / 2)
            p4 = surf.centralControlPoint();
          else
            p4 = surf.controlPoint(im, d - k - 1, j);
        }
        result.setControlPoint(i, j, k,
                               p1 * eta * theta + p2 * (1 - eta) * theta +
                               p3 * eta * (1 - theta) + p4 * (1 - eta) * (1 - theta));
      }
    }
  }

  // Set center as the mass center of the innermost control points
  Point3D center(0, 0, 0);
  for (size_t i = 0; i < n; ++i)
    center += result.controlPoint(i, l, l - 1);
  center /= n;
  result.setCentralControlPoint(center);

  result.setupLoop();
  return result;
}

SurfaceGeneralizedBezier loadBezier(const std::string &filename) {
  std::ifstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return SurfaceGeneralizedBezier();
  }

  size_t n, d;
  f >> n >> d;
  size_t l = (d + 1) / 2;
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;          // # of control points

  SurfaceGeneralizedBezier surf;
  surf.initNetwork(5, 5);

  Point3D p;
  f >> p[0] >> p[1] >> p[2];
  surf.setCentralControlPoint(p);

  for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    f >> p[0] >> p[1] >> p[2];
    surf.setControlPoint(side, col, row, p);
  }
  f.close();

  surf.setupLoop();

  return surf;
}

SurfaceGeneralizedBezier fitWithOriginal(const SurfaceGeneralizedBezier &original,
                                         const PointVector &points,
                                         const Point2DVector &params) {
  size_t n = original.domain()->vertices().size();
  size_t d = original.degree();
  size_t l = original.layers();
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;          // # of control points
  size_t bcp = n * (d - 1) * 2; // # of boundary control points
  size_t mcp = cp - bcp;        // # of movable control points
  size_t m = points.size();     // # of samples

  // Initialize patch with the fixed boundaries
  SurfaceGeneralizedBezier surf;
  surf.initNetwork(n, d);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < d - 1; ++j)
      for (size_t k = 0; k < 2; ++k)
        surf.setControlPoint(i, j, k, original.controlPoint(i, j, k));
  surf.setupLoop();

  Eigen::MatrixXd A(m, mcp);
  Eigen::MatrixXd b(m, 3);
  b.setZero();

  // Fill the matrices
  for (size_t j = 0; j < m; ++j) {
    double blend_sum = 0.0;
    for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
      if (col >= d - row) {
        if (++side >= n) {
          side = 0;
          ++row;
        }
        col = row;
      }
      double blend = surf.weight(side, col, row, params[j]);
      if (col < l)
        blend += surf.weight((side + n - 1) % n, d - row, col, params[j]);
      if (col > d - l)
        blend += surf.weight((side + 1) % n, row, d - col, params[j]);
      if (i > bcp)
        A(j, i - bcp) = blend;
      else {
        Point3D p = surf.controlPoint(side, col, row);
        b(j, 0) -= p[0] * blend; b(j, 1) -= p[1] * blend; b(j, 2) -= p[2] * blend;
      }
      blend_sum += blend;
    }
    A(j, 0) = 1.0 - blend_sum;
    b(j, 0) += points[j][0]; b(j, 1) += points[j][1]; b(j, 2) += points[j][2];
  }

  // LSQ Fit
  Eigen::MatrixXd x = A.fullPivLu().solve(b);

  // Fill control points
  surf.setCentralControlPoint(Point3D(x(0, 0), x(0, 1), x(0, 2)));
  for (size_t i = bcp + 1, side = 0, col = 2, row = 2; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    surf.setControlPoint(side, col, row, Point3D(x(i - bcp, 0), x(i - bcp, 1), x(i - bcp, 2)));
  }

  return surf;
}

void bezierTest(const std::string &filename) {
  SurfaceGeneralizedBezier surf = loadBezier("../../models/" + filename + ".gbp");

  // Generate mesh output
  TriMesh mesh = surf.eval(15);
  mesh.writeOBJ("../../models/bezier.obj");
  writeBezierControlPoints(surf, "../../models/bezier-cpts.obj");

  // Normal degree elevation
  SurfaceGeneralizedBezier sextic = elevateDegree(surf);
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
  surf = fitWithOriginal(sextic, mesh.points(), surf.domain()->parameters(15));
  surf.eval(15).writeOBJ("../../models/bezier-sextic.obj");
  writeBezierControlPoints(surf, "../../models/bezier-sextic-cpts.obj");
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
              << argv[0] << " bezier [model-name]" << std::endl;
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

  return 0;
}

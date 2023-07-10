#include <fstream>
#include <iostream>
#include <sstream>

#include "io.hh"

#include "domain.hh"
#include "utilities.hh"

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
    result.push_back(std::make_shared<BSplineCurve>(BSCurve(deg, knots, cpts)));
  }

  return result;
}

TriMesh readOBJ(const std::string &filename) {
  TriMesh result;
  std::ifstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return result;
  }
  bool points_set = false;
  std::string line;
  std::istringstream ss;
  Point3D p;
  TriMesh::Triangle t;
  PointVector pv;
  while (!f.eof()) {
    std::getline(f, line);
    if (line.empty())
      continue;
    switch (line[0]) {
    case 'v':
      ss.str(line);
      ss.seekg(2); // skip the first two characters
      ss >> p[0] >> p[1] >> p[2];
      pv.push_back(p);
      break;
    case 'f':
      if (!points_set) {
        result.setPoints(pv);
        points_set = true;
      }
      ss.str(line);
      ss.seekg(2); // skip the first two characters
      ss >> t[0] >> t[1] >> t[2];
      result.addTriangle(t[0] - 1, t[1] - 1, t[2] - 1);
      break;
    default:
      break;
    }
  }
  return result;
}

void writePCP(const std::string &filename, const Point2DVector &uvs, const PointVector &points) {
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }

  size_t n = uvs.size();
  f << n << std::endl;
  for (size_t i = 0; i < n; ++i) {
    const Point2D &uv = uvs[i];
    const Point3D &p = points[i];
    f << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << uv[0] << ' ' << uv[1] << std::endl;
  }

  f.close();
}

void loadBezier(const std::string &filename, SurfaceGeneralizedBezier *surf) {
  std::ifstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }

  size_t n, d;
  f >> n >> d;
  size_t l = (d + 1) / 2;
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;          // # of control points

  surf->initNetwork(n, d);

  Point3D p;
  f >> p[0] >> p[1] >> p[2];
  surf->setCentralControlPoint(p);

  for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    f >> p[0] >> p[1] >> p[2];
    surf->setControlPoint(side, col, row, p);
  }
  f.close();

  surf->setupLoop();
}

void saveBezier(const SurfaceGeneralizedBezier &surf, const std::string &filename) {
  std::ofstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return;
  }

  size_t n = surf.domain()->vertices().size();
  size_t d = surf.degree();
  f << n << ' ' << d << std::endl;
  size_t l = (d + 1) / 2;
  size_t cp = 1 + d / 2;
  cp = n * cp * l + 1;          // # of control points

  Point3D p = surf.centralControlPoint();
  f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;

  for (size_t i = 1, side = 0, col = 0, row = 0; i < cp; ++i, ++col) {
    if (col >= d - row) {
      if (++side >= n) {
        side = 0;
        ++row;
      }
      col = row;
    }
    p = surf.controlPoint(side, col, row);
    f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f.close();
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

SurfaceSPatch loadSPatch(const std::string &filename) {
  std::ifstream f(filename);
  f.exceptions(std::ios::failbit | std::ios::badbit);
  size_t n, d;
  f >> n >> d;
  size_t n_cp = binomial(n + d - 1, d);
  SurfaceSPatch surf;
  surf.initNetwork(n, d);
  SurfaceSPatch::Index index(n);
  for (size_t i = 0; i < n_cp; ++i) {
    for (size_t j = 0; j < n; ++j)
      f >> index[j];
    Point3D p;
    for (size_t j = 0; j < 3; ++j)
      f >> p[j];
    surf.setControlPoint(index, p);
  }
  surf.setupLoop();
  return surf;
}

std::vector<SurfaceSuperD> loadSuperDModel(const std::string &filename) {
  std::vector<SurfaceSuperD> result;
  std::ifstream f(filename);
  if (!f.is_open()) {
    std::cerr << "Unable to open file: " << filename << std::endl;
    return {};
  }
  size_t patches;
  f >> patches;
  for (size_t i = 0; i < patches; ++i) {
    SurfaceSuperD surf;
    size_t n;
    f >> n;
    surf.initNetwork(n);
    double x, y, z;
    for (size_t j = 0; j < n; ++j) {
      f >> x >> y >> z;
      surf.setFaceControlPoint(j, { x, y, z });
    }
    for (size_t j = 0; j < n; ++j) {
      f >> x >> y >> z;
      surf.setEdgeControlPoint(j, { x, y, z });
    }
    f >> x >> y >> z;
    surf.setVertexControlPoint({ x, y, z });
    surf.setupLoop();
    surf.updateRibbons();
    result.push_back(surf);
  }
  return result;
}

#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <sstream>

#include "domain.hh"
#include "ribbon.hh"
#include "surface-biharmonic.hh"
#include "surface-c0coons.hh"
#include "surface-composite-ribbon.hh"
#include "surface-corner-based.hh"
#include "surface-elastic.hh"
#include "surface-generalized-bezier-corner.hh"
#include "surface-generalized-bezier.hh"
#include "surface-generalized-coons.hh"
#include "surface-harmonic.hh"
#include "surface-hybrid.hh"
#include "surface-midpoint-coons.hh"
#include "surface-midpoint.hh"
#include "surface-nsided.hh"
#include "surface-polar.hh"
#include "surface-side-based.hh"
#include "surface-spatch.hh"
#include "surface-superd.hh"

#include "bezier.hh"
#include "gb-fit.hh"
#include "io.hh"

// Parameters
std::string filename;
double fullness = 0.5;
size_t resolution = 15;
double ribbon_length = 0.25;
double scaling = 20.0;

void ribbonTest(const std::shared_ptr<Surface> &surf) {
  // Fence
  TriMesh fence_mesh;
  size_t n = surf->domain()->size();
  size_t size = n * resolution * 2;
  PointVector pv; pv.reserve(size);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < resolution; ++j) {
      double u = (double)j / resolution;
      Point3D p = surf->ribbon(i)->curve()->eval(u);
      pv.push_back(p);
      p += surf->ribbon(i)->normal(u) * scaling;
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
      pv.push_back(surf->ribbon(i)->curve()->eval(u));
      pv.push_back(surf->ribbon(i)->eval(Point2D(u, ribbon_length)));
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

void showDeviations(const std::shared_ptr<Surface> &surf) {
  size_t n = surf->domain()->size();
  size_t res = 40;
  double step = 1.0e-4;
  std::shared_ptr<const Domain> domain = surf->domain();

  // Special handling of Generalized Bezier surfaces
  CurveVector inner_curves;
  bool is_bezier =
    bool(std::dynamic_pointer_cast<SurfaceGeneralizedBezier>(surf)) &&
    !bool(std::dynamic_pointer_cast<SurfaceHybrid>(surf));
  if (is_bezier) {
    auto bs = std::dynamic_pointer_cast<SurfaceGeneralizedBezier>(surf);
    auto degree = bs->degree();
    DoubleVector knots;
    knots.insert(knots.end(), degree + 1, 0.0);
    knots.insert(knots.end(), degree + 1, 1.0);
    PointVector cpts(degree + 1);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j <= degree; ++j)
        cpts[j] = bs->controlPoint(i, j, 1);
      inner_curves.push_back(std::make_shared<BSplineCurve>(BSCurve(degree, knots, cpts)));
    }
  }

  // Positional test
  double max_pos_error = 0.0;
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j <= res; ++j) {
      double s = (double)j / (double)res;
      Point2D uv = domain->edgePoint(i, s);
      Point3D p = surf->eval(uv), q = surf->ribbon(i)->curve()->eval(s);
      max_pos_error = std::max(max_pos_error, (p - q).norm());
    }
  }

  // Tangential test
  double max_tan_error = 0.0;
  VectorVector der;
  const Point2DVector &v = domain->vertices();
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
      surf->ribbon(i)->curve()->eval(s, 1, der);
      Vector3D surf_normal = (der[1] ^ (q - p)).normalize();
      Vector3D normal;
      if (is_bezier) {
        auto cross = inner_curves[i]->eval(s) - surf->ribbon(i)->curve()->eval(s);
        normal = (der[1] ^ cross).normalize();
      } else {
        normal = surf->ribbon(i)->normal(s);
      }
      double angle = acos(surf_normal * normal);
      max_tan_error = std::max(max_tan_error, angle);
    }
  }

  std::cout << "  positional error: " << max_pos_error << std::endl;
  std::cout << "  tangential error: " << max_tan_error * 180.0 / M_PI << std::endl;
}

void surfaceTest(std::string type, std::shared_ptr<Surface> &&surf) {
  CurveVector cv = readLOP("../../models/" + filename + ".lop");
  if (cv.empty())
    return;

  std::chrono::steady_clock::time_point begin, end;

  surf->setCurves(cv);
  surf->setupLoop();             // should be called after curve pointers changed
  surf->update();                // should be called after curves changed

  ribbonTest(surf);

  begin = std::chrono::steady_clock::now();
  surf->eval(resolution).writeOBJ("../../models/" + filename + "-" + type + ".obj");
  end = std::chrono::steady_clock::now();

  std::cout << type << ":" << std::endl;
  std::cout << "  evaluation time : "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;
  showDeviations(surf);

  if (type == "mp" || type == "mc") {
    size_t n = surf->domain()->size();
    Point3D midpoint(0,0,0);
    for (size_t i = 0; i < n; ++i)
      midpoint += surf->ribbon(i)->eval(Point2D(0.6, 0.4));
    midpoint /= n;
    dynamic_cast<SurfaceMidpoint *>(surf.get())->setMidpoint(midpoint);
    Point3D p = surf->eval(surf->domain()->center());
    std::cout << "  midpoint error: " << (p - midpoint).normSqr() << std::endl;
  }
}

void bezierTest() {
  auto surf = std::make_shared<SurfaceGeneralizedBezier>();
  loadBezier("../../models/" + filename + ".gbp", surf.get());

  std::chrono::steady_clock::time_point begin, end;
  begin = std::chrono::steady_clock::now();
  surf->eval(100).writeOBJ("../../models/" + filename + "-GB.obj");
  end = std::chrono::steady_clock::now();
  std::cout << "  evaluation time : "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;
  showDeviations(surf);
  
  // Generate mesh output
  TriMesh mesh = surf->eval(resolution);
  mesh.writeOBJ("../../models/bezier.obj");
  writeBezierControlPoints(*surf, "../../models/bezier-cpts.obj");

  // Normal degree elevation
  SurfaceGeneralizedBezier sextic = elevateDegree(*surf), sextic2, sextic3;
  writeBezierControlPoints(sextic, "../../models/bezier-elevated-cpts.obj");
  sextic.eval(resolution).writeOBJ("../../models/bezier-elevated.obj");
  SurfaceGeneralizedBezier elevated = *surf;
  for (size_t i = 0; i < 30; ++i)
    elevated = elevateDegree(elevated);
  writeBezierControlPoints(elevated, "../../models/bezier-elevated-30-times.obj");
  for (size_t i = 0; i < 30; ++i)
    elevated = elevateDegree(elevated);
  writeBezierControlPoints(elevated, "../../models/bezier-elevated-60-times.obj");

  // Fit a sextic surface on the quintic mesh
  sextic2 = fitWithOriginal(sextic, mesh.points(), surf->domain()->parameters(resolution));
  sextic2.eval(resolution).writeOBJ("../../models/bezier-sextic.obj");
  writeBezierControlPoints(sextic2, "../../models/bezier-sextic-cpts.obj");

  // Fit again, now with projected parameterization
  sextic3 = sextic;
  for (size_t i = 0; i < 1; ++i)
    sextic3 = fitWithOriginal(sextic3, mesh.points(), parameterizePoints(sextic3, mesh.points()));
  sextic3.eval(resolution).writeOBJ("../../models/bezier-sextic-projected.obj");
  writeBezierControlPoints(sextic3, "../../models/bezier-sextic-projected-cpts.obj");

  // Fit again, now with smoothing
  sextic3 = sextic;
  for (size_t i = 0; i < 1; ++i)
    sextic3 = fitWithOriginal(sextic3, mesh.points(), parameterizePoints(sextic3, mesh.points()),
                              0.2);
  sextic3.eval(resolution).writeOBJ("../../models/bezier-sextic-projected-smooth.obj");
  writeBezierControlPoints(sextic3, "../../models/bezier-sextic-projected-smooth-cpts.obj");
}

void cornerBezierTest() {
  auto surf = std::make_shared<SurfaceGeneralizedBezierCorner>();
  loadBezier("../../models/" + filename + ".gbp", surf.get());

  std::chrono::steady_clock::time_point begin, end;
  begin = std::chrono::steady_clock::now();
  surf->eval(100).writeOBJ("../../models/" + filename + "-GBC.obj");
  end = std::chrono::steady_clock::now();
  std::cout << "  evaluation time : "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;
  showDeviations(surf);
  
  // Generate mesh output
  TriMesh mesh = surf->eval(resolution);
  mesh.writeOBJ("../../models/bezier.obj");
  writeBezierControlPoints(*surf, "../../models/bezier-cpts.obj");
}

void hybridTest() {
  auto surf = std::make_shared<SurfaceHybrid>();
  loadBezier("../../models/" + filename + ".gbp", surf.get());
  surf->update();

  ribbonTest(surf);

  std::chrono::steady_clock::time_point begin, end;
  begin = std::chrono::steady_clock::now();
  surf->eval(resolution).writeOBJ("../../models/" + filename + "-HB.obj");
  end = std::chrono::steady_clock::now();
  std::cout << "  evaluation time : "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;
  showDeviations(surf);
  writeBezierControlPoints(*surf, "../../models/hybrid-cpts.obj");
}

void spatchTest() {
  auto surface = loadSPatch("../../models/" + filename + ".sp");

  std::chrono::steady_clock::time_point begin, end;
  begin = std::chrono::steady_clock::now();
  surface.eval(resolution).writeOBJ("../../models/" + filename + "-SP.obj");
  end = std::chrono::steady_clock::now();
  std::cout << "  evaluation time : "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;
}

void superDTest() {
  std::vector<SurfaceSuperD> surfaces = loadSuperDModel("../../models/" + filename + ".sdm");

  std::chrono::steady_clock::time_point begin, end;
  begin = std::chrono::steady_clock::now();
  for (size_t i = 0; i < surfaces.size(); ++i) {
    std::stringstream s;
    s << "../../models/" << filename << "-SD-" << i << ".obj";
    surfaces[i].eval(resolution).writeOBJ(s.str());
  }
  end = std::chrono::steady_clock::now();
  std::cout << "  evaluation time : "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count()
            << "ms" << std::endl;
}

void cloudTest() {
  CurveVector cv = readLOP("../../models/" + filename + ".lop");
  if (cv.empty())
    return;

  std::shared_ptr<Surface> surf = std::make_shared<SurfaceGeneralizedCoons>();
  surf->setCurves(cv);
  surf->setupLoop();
  surf->update();
  Point2DVector uvs = surf->domain()->parameters(resolution);
  PointVector points; points.reserve(uvs.size());
  std::transform(uvs.begin(), uvs.end(), std::back_inserter(points),
                 [surf](const Point2D &uv) { return surf->eval(uv); });

  writePCP("../../models/" + filename + ".pcp", uvs, points);
}

void meshFitTest(const std::string &surfname, const std::string &meshname) {
  SurfaceGeneralizedBezier surf;
  loadBezier("../../models/" + surfname + ".gbp", &surf);
  TriMesh mesh = readOBJ("../../models/" + meshname + ".obj");

  surf.eval(resolution).writeOBJ("../../models/bezier.obj");
  writeBezierControlPoints(surf, "../../models/bezier-cpts.obj");

  SurfaceGeneralizedBezier fitted;
  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()));
  fitted.eval(resolution).writeOBJ("../../models/bezier-fitted.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-cpts.obj");
  saveBezier(fitted, "../../models/bezier-fitted.gbp");

  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()),
                           0.2);
  fitted.eval(resolution).writeOBJ("../../models/bezier-fitted-smooth02.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-smooth02-cpts.obj");
  saveBezier(fitted, "../../models/bezier-smooth02.gbp");

  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()),
                           0.5);
  fitted.eval(resolution).writeOBJ("../../models/bezier-fitted-smooth05.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-smooth05-cpts.obj");
  saveBezier(fitted, "../../models/bezier-smooth05.gbp");

  fitted = fitWithOriginal(surf, mesh.points(), parameterizePoints(surf, mesh.points()),
                           1.0);
  fitted.eval(resolution).writeOBJ("../../models/bezier-fitted-smooth10.obj");
  writeBezierControlPoints(fitted, "../../models/bezier-fitted-smooth10-cpts.obj");
  saveBezier(fitted, "../../models/bezier-smooth10.gbp");
}

void printStatistics(const DoubleVector &data) {
  size_t n = data.size();
  double mean = std::accumulate(data.begin(), data.end(), 0.0) / n;
  DoubleVector dev; dev.reserve(n);
  std::transform(data.begin(), data.end(), std::back_inserter(dev),
                 [mean](double x) { return pow(x - mean, 2); });
  double stddev = sqrt(std::accumulate(dev.begin(), dev.end(), 0.0) / (n - 1));
  DoubleVector sorted = data; std::sort(sorted.begin(), sorted.end());
  double quartiles[5];
  quartiles[0] = sorted.front(); quartiles[4] = sorted.back();
  for (size_t i = 1; i < 4; ++i) {
    if (n % 4 == 0)
      quartiles[i] = sorted[n * i / 4];
    else if (n % 4 == 1)
      quartiles[i] = sorted[n * i / 4] * 0.75 + sorted[n * i / 4 + 1] * 0.25;
    else if (n % 4 == 2)
      quartiles[i] = sorted[n * i / 4] * 0.5 + sorted[n * i / 4 + 1] * 0.5;
    else
      quartiles[i] = sorted[n * i / 4] * 0.25 + sorted[n * i / 4 + 1] * 0.75;
  }
  std::cout << "Statistics:" << std::endl
            << "Mean:\t" << mean << "\nStd.Dev.:\t" << stddev << std::endl
            << "Min:\t" << quartiles[0] << "\nMax:\t" << quartiles[4] << std::endl
            << "Quartiles:\t" << quartiles[1] << ", " << quartiles[2]
            << ", " << quartiles[3] << std::endl
            << "Inter-quartile histogram:" << std::endl;
  size_t bins = 25, height = 15;
  double bin_width = (quartiles[3] - quartiles[1]) / bins;
  std::vector<size_t> histogram(bins, 0);
  size_t current_bin = 0;
  double current_bin_min = quartiles[1];
  double current_bin_max = quartiles[1] + bin_width;
  for (double x : sorted) {
    if (x < current_bin_min)
      continue;
    if (x > current_bin_max) {
      ++current_bin;
      current_bin_min += bin_width;
      current_bin_max += bin_width;
    }
    if (current_bin >= bins)
      break;
    ++histogram[current_bin];
  }
  size_t max = *std::max_element(histogram.begin(), histogram.end());
  std::transform(histogram.begin(), histogram.end(), histogram.begin(),
                 [height,max](size_t k) { return k * height / max; });
  for (size_t i = height; i > 0; --i) {
    for (size_t j = 0; j < bins; ++j)
      if (histogram[j] >= i)
        std::cout << " xx";
      else
        std::cout << "   ";
    std::cout << std::endl;
  }
  for (size_t j = 0; j < bins; ++j)
    std::cout << "L__";
  std::cout << "I" << std::endl;
}

DoubleVector deviationFromMesh(const Surface &surf, const TriMesh &mesh) {
  size_t resolution = 30;
  TriMesh surf_mesh = surf.eval(resolution);
  PointVector vertices = surf_mesh.points();
  DoubleVector result;
  for (const auto &p : mesh.points()) {
    TriMesh::Triangle tri = surf_mesh.closestTriangle(p);
    const Point3D &a = vertices[tri[0]], &b = vertices[tri[1]], &c = vertices[tri[2]];
    Vector3D n = ((b - a) ^ (c - a)).normalize();
    Point3D q = p + n * ((a - p) * n);
    result.push_back((p - q).norm());
  }
  return result;
}

void deviationTest(const std::string &surfname, const std::string &meshname) {
  // SurfaceGeneralizedBezier surf;
  // loadBezier("../../models/" + surfname + ".gbp", &surf);
  SurfaceMidpointCoons surf;
  {
    CurveVector cv = readLOP("../../models/" + surfname + ".lop");
    if (cv.empty())
      return;
    surf.setCurves(cv);
    surf.setupLoop();             // should be called after curve pointers changed
    surf.update();                // should be called after curves changed
  }

  TriMesh mesh = readOBJ("../../models/" + meshname + ".obj");
  DoubleVector deviations = deviationFromMesh(surf, mesh);
  printStatistics(deviations);
}

void classATest() {
  CurveVector cv = readLOP("../../models/pocket6sided.lop");

  CurveVector normal, class_a;
  for (const auto &c : cv) {
    VectorVector der;
    Point3D pa = c->eval(0.0, 1, der);
    Vector3D va = der[1];
    Point3D pb = c->eval(1.0, 1, der);
    Vector3D vb = der[1];
    PointVector pv = {pa, pa + va / 3.0, pb - vb / 3.0, pb};
    normal.push_back(std::make_shared<BSplineCurve>(BSCurve(pv)));
    BCurve bc; bc.fitClassA(14, pa, va, pb, vb);
    class_a.push_back(std::make_shared<BSplineCurve>(BSCurve(bc.controlPoints())));
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
}

int main(int argc, char **argv) {
#ifdef DEBUG
  std::cout << "Compiled in DEBUG mode" << std::endl;
#else  // !DEBUG
  std::cout << "Compiled in RELEASE mode" << std::endl;
#endif  // DEBUG

  if (argc < 2 || argc > 5) {
    std::cerr << "Usage:\n"
              << argv[0] << " model-name [resolution] [fence-scaling] [ribbon-length]" << std::endl
              << argv[0] << " bezier [model-name]" << std::endl
              << argv[0] << " bezier-corner [model-name]" << std::endl
              << argv[0] << " hybrid [model-name]" << std::endl
              << argv[0] << " cloud [model-name]" << std::endl
              << argv[0] << " class-a" << std::endl
              << argv[0] << " mesh-fit [model-name] [mesh-name]" << std::endl
              << argv[0] << " deviation [model-name] [mesh-name]" << std::endl
              << argv[0] << " spatch [model-name]" << std::endl
              << argv[0] << " superd [model-name]" << std::endl;
    return 1;
  }

  std::string type = argv[1];
  if (argc == 2)
    filename = "cagd86";
  else
    filename = argv[2];

  if (type == "bezier")
    bezierTest();
  else if (type == "bezier-corner")
    cornerBezierTest();
  else if (type == "hybrid")
    hybridTest();
  else if (type == "cloud")
    cloudTest();
  else if (type == "class-a")
    classATest();
  else if (type == "mesh-fit" || type == "deviation") {
    if (argc < 4) {
      std::cerr << "Not enough parameters!" << std::endl;
      return 1;
    }
    if (type == "mesh-fit")
      meshFitTest(filename, argv[3]);
    else
      deviationTest(filename, argv[3]);
  } else if (type == "spatch")
    spatchTest();
  else if (type == "superd") {
    if (argc == 2)
      filename = "trebol";
    if (argc > 3)
      fullness = strtod(argv[3], nullptr);
    superDTest();
  } else {

    // First argument treated as filename
    filename = type;
    if (argc >= 3)
      resolution = atoi(argv[2]);
    if (argc >= 4)
      scaling = strtod(argv[3], nullptr);
    if (argc >= 5)
      ribbon_length = strtod(argv[4], nullptr);

    // surfaceTest("sb", std::make_shared<SurfaceSideBased>());
    // surfaceTest("cb", std::make_shared<SurfaceCornerBased>());
    // surfaceTest("gc", std::make_shared<SurfaceGeneralizedCoons>());
    // surfaceTest("cr", std::make_shared<SurfaceCompositeRibbon>());
    // surfaceTest("mp", std::make_shared<SurfaceMidpoint>());
    surfaceTest("mc", std::make_shared<SurfaceMidpointCoons>());
    // surfaceTest("pp", std::make_shared<SurfacePolar>());
    // surfaceTest("ns", std::make_shared<SurfaceNSided>());
    // surfaceTest("cc", std::make_shared<SurfaceC0Coons>());
    // surfaceTest("ep", std::make_shared<SurfaceElastic>());
    // surfaceTest("hp", std::make_shared<SurfaceHarmonic>());
    // surfaceTest("bp", std::make_shared<SurfaceBiharmonic>());
  }
}

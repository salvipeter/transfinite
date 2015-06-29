#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "rmf.hh"
#include "domain-regular.hh"
#include "domain-circular.hh"
#include "parameterization-bilinear.hh"
#include "parameterization-barycentric.hh"
#include "ribbon-compatible.hh"
#include "surface-side-based.hh"

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

void surfaceTest(std::shared_ptr<Surface> &&surf, const CurveVector &cv,
                 std::string filename, size_t resolution) {
  surf->setCurves(cv);
  surf->setupLoop();             // should be called after curve pointers changed
  surf->update();                // should be called after curves changed
  surf->eval(resolution).writeOBJ(filename);
}

int main(int argc, char **argv) {
  Vector2D u(1, 3), v(4, 5);
  u += v;
  Vector2D w = u + v;
  std::cout << w[0] << ", " << w[1] << std::endl;

  BSCurve c(2, {0,0,0,1,2,4,4,4}, {{0,0,0},{1,1,0},{2,2,0},{3,1,0},{4,0,0}});
  VectorVector der;
  Point3D p = c.eval(0.5, 2, der);
  std::cout << p[0] << ", " << p[1] << std::endl;
  std::cout << der[1][0] << ", " << der[1][1] << std::endl;
  std::cout << der[2][0] << ", " << der[2][1] << std::endl;
  c.reverse();
  c.normalize();
  std::ofstream f("/tmp/test");
  for (size_t i = 0; i <= 100; ++i) {
    double u = (double)i / 100.0;
    Point3D p = c.eval(u);
    f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f.close();
  std::cout << std::setprecision(15) << c.arcLength(0.0, 0.3) << ", " <<
    c.arcLength(0.5, 0.7) << ", " << c.arcLength(0.8, 1.0) << std::endl;
  
  // RMF test
  DoubleVector knots =
    {0, 0, 0, 0, 0.133053, 0.232419, 0.338451, 0.43994, 0.541278, 0.667394, 1, 1, 1, 1};
  PointVector cpts =
    {{0.155279, -2.14286, 0}, {0.0532173, -1.80779, 0}, {0.190782, -1.07383, 1},
     {1.09142, -0.69642, 1}, {1.81621, -0.426967, 2}, {2.81035, -0.467704, 1},
     {3.07214, 0.640326, 1}, {5.1987, -0.432608, 0}, {4.46645, -2.13998, 1},
     {3.25645, -1.96703, -1}};
  std::shared_ptr<BSCurve> c2 = std::make_shared<BSCurve>(3, knots, cpts);
  Vector3D start(0.9041749728349117, 0.27541001827971723, -0.32652249590212407);
  Vector3D end(0.6838937155616026, -0.5642165240201893, -0.46254632182941535);
  
  RMF rmf;
  rmf.setCurve(c2);
  rmf.setStart(start);
  rmf.setEnd(end);
  rmf.update();
  size_t res = 57;
  std::ofstream f2("/tmp/rmftest");
  for (size_t i = 0; i <= res; ++i) {
    double u = (double)i / (double)res;
    Point3D p = c2->eval(u);
    Vector3D n = rmf.eval(u);
    f2 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    f2 << p[0] + n[0] << ' ' << p[1] + n[1] << ' ' << p[2] + n[2] << std::endl;
    f2 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f2.close();
  c2->reverse();
  c2->normalize();
  rmf.setStart(end);
  rmf.setEnd(start);
  rmf.update();
  std::ofstream f3("/tmp/rmftest2");
  for (size_t i = 0; i <= res; ++i) {
    double u = (double)i / (double)res;
    Point3D p = c2->eval(u);
    Vector3D n = rmf.eval(u);
    f3 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    f3 << p[0] + n[0] << ' ' << p[1] + n[1] << ' ' << p[2] + n[2] << std::endl;
    f3 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f3.close();

  // Domain test
  const size_t mesh_res = 15;
  std::shared_ptr<Domain> d = std::make_shared<DomainCircular>();
  std::shared_ptr<BSCurve> cs = std::make_shared<BSCurve>(c);
  d->setSides({cs,cs,c2,cs,c2});
  d->update();
  std::shared_ptr<Parameterization> par = std::make_shared<ParameterizationBarycentric>();
  par->setDomain(d);
  par->update();
  const Point2DVector &params = d->parameters(mesh_res);
  TriMesh mesh = d->meshTopology(mesh_res);
  PointVector points(params.size());
  std::transform(params.begin(), params.end(), points.begin(),
                 [par](const Point2D &p) {
                   return Point3D(p[0], p[1], par->mapToRibbon(1, p)[0]);
                 });
  mesh.setPoints(points);
  mesh.writeOBJ("/tmp/test-s.obj");
  std::transform(params.begin(), params.end(), points.begin(),
                 [par](const Point2D &p) {
                   return Point3D(p[0], p[1], par->mapToRibbon(1, p)[1]);
                 });
  mesh.setPoints(points);
  mesh.writeOBJ("/tmp/test-d.obj");

  // Ribbon test
  CurveVector cv = readLOP("../../models/cagd86.lop");
  size_t n = cv.size();
  double max_d = 0.25;
  size_t rib_res = 15;
  TriMesh mesh2;
  mesh2.resizePoints((rib_res + 1) * (rib_res + 1) * n);
  PointVector pv; pv.reserve((rib_res + 1) * (rib_res + 1) * n);
  size_t index = 0;
  std::vector<std::shared_ptr<Ribbon>> ribbons;
  for (size_t side = 0; side < n; ++side) {
    ribbons.push_back(std::make_shared<RibbonCompatible>());
    ribbons.back()->setCurve(cv[side]);
  }
  for (size_t side = 0; side < n; ++side) {
    ribbons[side]->setNeighbors(ribbons[(side+n-1)%n], ribbons[(side+1)%n]);
    ribbons[side]->update();

    for (size_t si = 0; si <= rib_res; ++si) {
      double s = (double)si / (double)rib_res;
      for (size_t di = 0; di <= rib_res; ++di) {
        double d = (double)di / (double)rib_res;
        d *= max_d;
        pv.push_back(ribbons[side]->eval(Point2D(s, d)));

        if (si != 0 && di != 0) {
          mesh2.addTriangle(index, index - rib_res - 1, index - 1);
          mesh2.addTriangle(index - rib_res - 1, index - rib_res - 2, index - 1);
        }

        ++index;
      }
    }
    mesh2.setPoints(pv);
    mesh2.writeOBJ("/tmp/ribtest.obj");
  }

  // Surface test
  surfaceTest(std::make_shared<SurfaceSideBased>(), cv, "/tmp/surf-sb.obj", 15);
  // surfaceTest(std::make_shared<SurfaceCornerBased>(), cv, "/tmp/surf-cb.obj", 15);
  // surfaceTest(std::make_shared<SurfaceGeneralizedCoons>(), cv, "/tmp/surf-gc.obj", 15);
  // surfaceTest(std::make_shared<SurfaceCompositeRibbon>(), cv, "/tmp/surf-cr.obj", 15);

  return 0;
}

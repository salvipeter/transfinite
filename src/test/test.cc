#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "rmf.hh"
#include "domain-regular.hh"
#include "domain-circular.hh"
#include "parameterization-bilinear.hh"
#include "parameterization-barycentric.hh"

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
  BSCurve c2(3, {0, 0, 0, 0, 0.133053, 0.232419, 0.338451, 0.43994, 0.541278, 0.667394, 1, 1, 1, 1},
             {{0.155279, -2.14286, 0}, {0.0532173, -1.80779, 0}, {0.190782, -1.07383, 1},
              {1.09142, -0.69642, 1}, {1.81621, -0.426967, 2}, {2.81035, -0.467704, 1},
              {3.07214, 0.640326, 1}, {5.1987, -0.432608, 0}, {4.46645, -2.13998, 1},
              {3.25645, -1.96703, -1}});
  Vector3D start(0.9041749728349117, 0.27541001827971723, -0.32652249590212407);
  Vector3D end(0.6838937155616026, -0.5642165240201893, -0.46254632182941535);
  RMF rmf(c2, start, end);
  size_t res = 57;
  std::ofstream f2("/tmp/rmftest");
  for (size_t i = 0; i <= res; ++i) {
    double u = (double)i / (double)res;
    Point3D p = c2.eval(u);
    Vector3D n = rmf.eval(u);
    f2 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    f2 << p[0] + n[0] << ' ' << p[1] + n[1] << ' ' << p[2] + n[2] << std::endl;
    f2 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f2.close();
  c2.reverse();
  c2.normalize();
  RMF rmf2(c2, end, start);
  std::ofstream f3("/tmp/rmftest2");
  for (size_t i = 0; i <= res; ++i) {
    double u = (double)i / (double)res;
    Point3D p = c2.eval(u);
    Vector3D n = rmf2.eval(u);
    f3 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
    f3 << p[0] + n[0] << ' ' << p[1] + n[1] << ' ' << p[2] + n[2] << std::endl;
    f3 << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f3.close();

  // Domain test
  const size_t mesh_res = 15;
  std::shared_ptr<Domain> d = std::make_shared<DomainCircular>();
  d->setSides({c,c,c2,c,c2});
  std::shared_ptr<Parameterization> par = std::make_shared<ParameterizationBarycentric>();
  par->setDomain(d);
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

  return 0;
}

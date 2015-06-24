#include <fstream>
#include <iostream>

#include "geometry.hh"

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
  std::ofstream f("/tmp/test2");
  for (size_t i = 0; i <= 100; ++i) {
    double u = (double)i / 100.0;
    Point3D p = c.eval(u);
    f << p[0] << ' ' << p[1] << ' ' << p[2] << std::endl;
  }
  f.close();
  
  return 0;
}

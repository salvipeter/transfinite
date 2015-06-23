#include <iostream>

#include "geometry.hh"

int main(int argc, char **argv) {
  Vector2D u(1, 3), v(4, 5);
  u += v;
  Vector2D w = u + v;
  std::cout << w[0] << ", " << w[1] << std::endl;

  BSCurve c(3, {0,0,0,0,1,1,1,1}, {{0,0,0},{1,1,0},{2,1,0},{3,0,0}});
  VectorVector der;
  Point3D p = c.eval(0.5, 2, der);
  std::cout << p[0] << ", " << p[1] << std::endl;
  std::cout << der[1][0] << ", " << der[1][1] << std::endl;
  std::cout << der[2][0] << ", " << der[2][1] << std::endl;
  
  return 0;
}

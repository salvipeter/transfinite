#include <iostream>

#include "geometry.hh"

int main(int argc, char **argv) {
  Vector2D u(1, 3), v(4, 5);
  u += v;
  Vector2D w = u + v;
  std::cout << w[0] << ", " << w[1] << std::endl;
  return 0;
}

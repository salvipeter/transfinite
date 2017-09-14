#include <cmath>

double
blendHermite(double x) {
  double x2 = x * x;
  return 2.0 * x * x2 - 3.0 * x2 + 1.0;
}

double
hermite(int i, double t) {
  switch(i) {
  case 0: return std::pow(1 - t, 3) + 3.0 * std::pow(1 - t, 2) * t;
  case 1: return std::pow(1 - t, 2) * t;
  case 2: return (1 - t) * std::pow(t, 2);
  case 3: return 3.0 * (1 - t) * std::pow(t, 2) + std::pow(t, 3);
  default: ;
  }
  return -1.0;                  // should not come here
}

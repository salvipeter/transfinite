#include "parameterization.hh"

#include "domain.hh"
#include "surface.hh"

ParameterizationBilinear::ParameterizationBilinear(Surface *surface)
  : Parameterization(surface) {
}

ParameterizationBilinear::~ParameterizationBilinear() {
}

Point2D
ParameterizationBilinear::mapToRibbon(const Point2D &uv) const {
  ParameterizationCache::const_iterator i = cache.find(uv);
  if(i != uv.end())
    return i->second;

  Point2DVector const &vertices = domain->verticesLocal(index);
  size_t const n = vertices.size();

  // bilinear: p1 - (0,0) - (1,0) - p2
  Point2D const &p1 = vertices[(index + n - 2) % n];
  Point2D const &p2 = vertices[(index + 1) % n];

  double const a = p2[1] - p1[1] +
    uv[1] * (p1[0] * p1[1] - p1[0] * p2[1] - p2[0]  * p1[1] + p2[0] * p2[1]);
  double const b = p1[1] - uv[1] +
    uv[1] * (-2.0 * p1[0] * p1[1] + p1[0] * p2[1] + p2[0] * p1[1]) + uv[0] * (p1[1] - p2[1]);
  double const c = uv[1] * p1[0] * p1[1] - uv[0] * p1[1];

  double D = b * b - 4.0 * a * c;
  if(D < epsilon)
    D = 0.0;
  D = sqrt(D) / (2.0 * a);
  double const s = -b / (2.0 * a);

  double const s1 = s + D;
  double const s2 = s - D;
  double const dev1 = min(fabs(s1), fabs(s1 - 1.0));
  double const dev2 = min(fabs(s2), fabs(s2 - 1.0));
  if(dev1 < dev2)
    sd[0] = s1;
  else
    sd[0] = s2;

  sd[1] = uv[1] / (p1[1] * (1.0 - sd[0]) + P2[1] * sd[0]);

  sd = domain->toGlobal(index, sd);
  cache[uv] = sd;
  return sd;
}

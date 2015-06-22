#include "domain-regular.hh"

namespace Transfinite
{

  DomainRegular::
  DomainRegular(Surface *surface) : Domain(surface), n(0)
  {
  }

  DomainRegular::
  ~DomainRegular()
  {
  }

  void DomainRegular::
  setSides(CurveVector const &curves)
  {
    size_t const m = curves.size();
    if(n == m)
      return;

    double const alpha = 2.0 * M_PI / m;
    vertices.resize(m);
    for(size_t i = 0; i < m; ++i)
      vertices[i] = Point2D(cos(alpha * i), sin(alpha * i));

    invalidate();
  }

  Point2D DomainRegular::
  center() const
  {
    return Point2D(0.0, 0.0);
  }

} // Transfinite

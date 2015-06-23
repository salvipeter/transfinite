#pragma once

#include "domain.hh"

namespace Transfinite
{

  class DomainRegular : public Domain
  {
  public:
    DomainRegular(Surface *surface);
    virtual ~DomainRegular();
    virtual void setSides(CurveVector const &curves);
    virtual Point2D center() const;
  };

} // Transfinite

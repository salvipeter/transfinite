#pragma once

#include "domain.hh"

namespace Transfinite {

class DomainRegular : public Domain {
public:
  virtual ~DomainRegular();
  virtual bool update();
  virtual void computeCenter();
};

} // namespace Transfinite

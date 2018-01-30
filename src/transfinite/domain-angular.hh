#pragma once

#include "domain-circular.hh"

namespace Transfinite {

class DomainAngular : public DomainCircular {
public:
  virtual ~DomainAngular();
  virtual bool update() override;
};

} // namespace Transfinite

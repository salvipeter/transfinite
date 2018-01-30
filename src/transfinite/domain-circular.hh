#pragma once

#include "domain.hh"

namespace Transfinite {

class DomainCircular : public Domain {
public:
  virtual ~DomainCircular();
  virtual bool update() override;
};

} // namespace Transfinite

#pragma once

#include "domain.hh"

namespace Transfinite {

class DomainRegular : public Domain {
public:
  virtual ~DomainRegular();
  virtual bool update() override;
  virtual void computeCenter() override;
};

} // namespace Transfinite

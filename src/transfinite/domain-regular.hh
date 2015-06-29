#pragma once

#include "domain.hh"

class DomainRegular : public Domain {
public:
  virtual ~DomainRegular();
  virtual bool update();
  virtual void computeCenter();
};

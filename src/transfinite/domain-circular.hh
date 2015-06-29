#pragma once

#include "domain.hh"

class DomainCircular : public Domain {
public:
  virtual ~DomainCircular();
  virtual bool update();
  virtual void computeCenter();
};

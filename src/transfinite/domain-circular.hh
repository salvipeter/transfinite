#pragma once

#include "domain.hh"

class DomainCircular : public Domain {
public:
  virtual ~DomainCircular();
  virtual void setSides(const CurveVector &curves);
  virtual void computeCenter();
};

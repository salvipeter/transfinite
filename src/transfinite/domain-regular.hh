#pragma once

#include "domain.hh"

class DomainRegular : public Domain {
public:
  virtual ~DomainRegular();
  virtual void setSides(const CurveVector &curves);
  virtual void computeCenter();
};

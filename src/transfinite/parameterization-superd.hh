#pragma once

#include "parameterization-parallel.hh"

namespace Transfinite {

class ParameterizationSuperD : public ParameterizationParallel {
public:
  virtual ~ParameterizationSuperD();
  virtual void updateMultipliers() override;
};

} // namespace Transfinite

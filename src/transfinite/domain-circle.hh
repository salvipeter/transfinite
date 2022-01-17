#pragma once

#include "domain.hh"

namespace Transfinite {

class DomainCircle : public Domain {
public:
  virtual ~DomainCircle();
  virtual bool update() override;
  virtual size_t meshSize(size_t resolution) const override;
  virtual TriMesh meshTopology(size_t resolution) const override;
  virtual bool onEdge(size_t resolution, size_t index) const override;
  virtual Point2D edgePoint(size_t i, double s) const override;
protected:
  virtual void computeCenter() override;
  virtual Point2DVector parametersImpl(size_t resolution) const override;
};

} // namespace Transfinite

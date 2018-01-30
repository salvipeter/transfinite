#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceMidpoint : public Surface {
public:
  SurfaceMidpoint();
  SurfaceMidpoint(const SurfaceMidpoint &) = default;
  virtual ~SurfaceMidpoint();
  SurfaceMidpoint &operator=(const SurfaceMidpoint &) = default;
  virtual void update(size_t i) override;
  virtual void update() override;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;
  void setMidpoint(const Point3D &p);
  void unsetMidpoint();

protected:
  virtual double deficiency(const Point2D &p) const;
  virtual std::shared_ptr<Ribbon> newRibbon() const override;

  Point3D central_cp_;

private:
  void updateCentralControlPoint();

  bool midpoint_set_;
  Point3D midpoint_;
};

} // namespace Transfinite

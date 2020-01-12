#pragma once

#include "surface.hh"

namespace Transfinite {

class SurfaceSuperD : public Surface {
public:
  using QuarticCurve = std::array<Point3D, 5>;
  using QuarticSurface = std::array<QuarticCurve, 5>;
  SurfaceSuperD();
  SurfaceSuperD(const SurfaceSuperD &) = default;
  virtual ~SurfaceSuperD();
  SurfaceSuperD &operator=(const SurfaceSuperD &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;
  void initNetwork(size_t n);
  virtual void setupLoop() override;
  void updateRibbons();
  Point3D vertexControlPoint() const;
  void setVertexControlPoint(const Point3D &p);
  Point3D faceControlPoint(size_t i) const;
  void setFaceControlPoint(size_t i, const Point3D &p);
  Point3D edgeControlPoint(size_t i) const;
  void setEdgeControlPoint(size_t i, const Point3D &p);
  double fullness() const;
  void setFullness(double f);

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
  // Functions generating the "real" ribbon:
  QuarticCurve generateQuartic(const Point3D &a, const Point3D &b, const Point3D &c) const;
  QuarticCurve generateBase(size_t i) const;
  QuarticCurve generateMid(size_t i) const;
  QuarticCurve generateOpp(size_t i) const;
  QuarticSurface generateRibbon(size_t i) const;

private:
  double fullness_;
  Point3D cp_v_;
  PointVector cp_f_, cp_e_;
  std::vector<QuarticSurface> quartic_ribbons_;
};

} // namespace Transfinite
